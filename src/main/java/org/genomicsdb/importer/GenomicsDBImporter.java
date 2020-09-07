/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.genomicsdb.importer;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import org.genomicsdb.Constants;
import org.genomicsdb.GenomicsDBLibLoader;
import org.genomicsdb.exception.GenomicsDBException;
import org.genomicsdb.importer.extensions.CallSetMapExtensions;
import org.genomicsdb.importer.extensions.JsonFileExtensions;
import org.genomicsdb.importer.extensions.VidMapExtensions;
import org.genomicsdb.importer.model.ChromosomeInterval;
import org.genomicsdb.importer.model.SampleInfo;
import org.genomicsdb.model.*;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.genomicsdb.GenomicsDBUtils.createTileDBWorkspace;
import static org.genomicsdb.GenomicsDBUtils.listGenomicsDBFragments;
import static org.genomicsdb.GenomicsDBUtils.deleteFile;
import static org.genomicsdb.GenomicsDBUtils.writeToFile;
import static org.genomicsdb.GenomicsDBUtils.getMaxValidRowIndex;

/**
 * Java wrapper for vcf2genomicsdb - imports VCFs into GenomicsDB.
 * All vid information is assumed to be set correctly by the user (JSON files)
 */
public class GenomicsDBImporter extends GenomicsDBImporterJni implements JsonFileExtensions, CallSetMapExtensions,
        VidMapExtensions {
    static {
        try {
            boolean loaded = GenomicsDBLibLoader.loadLibrary();
            if (!loaded) throw new GenomicsDBException("Could not load genomicsdb native library");
        } catch (UnsatisfiedLinkError ule) {
            throw new GenomicsDBException("Could not load genomicsdb native library", ule);
        }
    }

    private ImportConfig config;
    private String mLoaderJSONFile = null;
    private int mRank = 0;
    //For buffered streams
    private boolean mContainsBufferStreams = false;
    private long mGenomicsDBImporterObjectHandle = 0;
    private ArrayList<GenomicsDBImporterStreamWrapper> mBufferStreamWrapperVector = new ArrayList<GenomicsDBImporterStreamWrapper>();
    private boolean mIsLoaderSetupDone = false;
    private long[] mExhaustedBufferStreamIdentifiers = null;
    private long mNumExhaustedBufferStreams = 0;
    //Done flag - useful only for buffered streams
    private boolean mDone = false;
    //JSON object that specifies callset/sample name to row_idx mapping in the buffer
    private JSONObject mCallsetMappingJSON = new JSONObject();

    /**
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     */
    public GenomicsDBImporter(final String loaderJSONFile) {
        initialize(loaderJSONFile, 0);
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public GenomicsDBImporter(final String loaderJSONFile, final int rank) {
        initialize(loaderJSONFile, rank);
    }

    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param importConfig      Top level import configuration object
     * @param sampleToReaderMap Feature Readers objects corresponding to input GVCF files
     * @param rank              Rank of object - corresponds to the partition index in the loader
     * @throws IOException when load into TileDB array fails
     */
    private GenomicsDBImporter(final ImportConfig importConfig,
                               final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                               final int rank) throws IOException {
        File importJSONFile = dumpTemporaryLoaderJSONFile(importConfig.getImportConfiguration(), "");

        initialize(importJSONFile.getAbsolutePath(), rank);
   
        //jniCopyVidMap(mGenomicsDBImporterObjectHandle, vidMapPB.toByteArray());
        //jniCopyCallsetMap(mGenomicsDBImporterObjectHandle, callsetMappingPB.toByteArray());

        GenomicsDBImportConfiguration.Partition partition = importConfig.getImportConfiguration().getColumnPartitions(rank);
        String chromosomeName = "";
        int chromosomeStart = 0;
        int chromosomeEnd = 0;
        if (partition.getBegin().hasContigPosition()) {
            chromosomeName = partition.getBegin().getContigPosition().getContig();
            chromosomeStart = (int) partition.getBegin().getContigPosition().getPosition();
            chromosomeEnd = (int) partition.getEnd().getContigPosition().getPosition();
        }
        
        if(sampleToReaderMap != null)
            for (final Map.Entry<String, FeatureReader<VariantContext>> currEntry : sampleToReaderMap.entrySet()) {

                String sampleName = currEntry.getKey();

                FeatureReader<VariantContext> featureReader = currEntry.getValue();

                Iterator<VariantContext> iterator;
                if (chromosomeName != null && !chromosomeName.isEmpty()) {
                    iterator = featureReader.query(chromosomeName, chromosomeStart, chromosomeEnd);
                }
                else {
                    // we may have multiple chromosomes per partition...
                    try {
                        iterator = GenomicsDBImporter.columnPartitionIterator(
                                featureReader, 
                                importJSONFile.getAbsolutePath(), rank);
                    }
                    catch (ParseException ex) {
                        throw new RuntimeException("Error parsing import json file:"+ex);
                    }
                }

                addSortedVariantContextIterator(
                        getStreamNameFromSampleName(sampleName),
                        (VCFHeader) featureReader.getHeader(),
                        iterator,
                        importConfig.getImportConfiguration().getSizePerColumnPartition()/sampleToReaderMap.size(),
                        importConfig.isPassAsVcf() ? VariantContextWriterBuilder.OutputType.VCF_STREAM
                                : VariantContextWriterBuilder.OutputType.BCF_STREAM,
                        null);
            }
    }

    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for running Chromosome intervals imports in parallel.
     *
     * @param config Parallel import configuration
     * @throws FileNotFoundException when files could not be read/written
     * @throws com.googlecode.protobuf.format.JsonFormat.ParseException when existing callset jsons are invalid. incremental case
     */
    public GenomicsDBImporter(final ImportConfig config) throws FileNotFoundException, 
            com.googlecode.protobuf.format.JsonFormat.ParseException {
        this.config = config;
        this.config.setImportConfiguration(addExplicitValuesToImportConfiguration(config));
        long lbRowIdx = this.config.getImportConfiguration().hasLbCallsetRowIdx()
            ? this.config.getImportConfiguration().getLbCallsetRowIdx()
            : 0L;
        // if lower bound row index is zero, and this is incremental import
        // then let's try to read from metadata 
        String workspace = this.config.getImportConfiguration().getColumnPartitions(0).getWorkspace();
        if (lbRowIdx == 0 && this.config.isIncrementalImport()) {
            GenomicsDBImportConfiguration.Partition partition = this.config.getImportConfiguration().getColumnPartitions(0); 
            String array;
            if (partition.hasGenerateArrayNameFromPartitionBounds()) {
                if (partition.getBegin().hasContigPosition()) {
                    final String chromosomeName = partition.getBegin().getContigPosition().getContig();
                    final int chromosomeStart = (int) partition.getBegin().getContigPosition().getPosition();
                    final int chromosomeEnd = (int) partition.getEnd().getContigPosition().getPosition();
                    array = String.format(Constants.CHROMOSOME_INTERVAL_FOLDER, chromosomeName, chromosomeStart, chromosomeEnd);
                }
                else {
                    // if array name isn't using contig position, then just named using
                    // partition rank. Let's pick "0" to get max row index
                    array = "0";
                }
            } else {
                array = partition.getArrayName();
            }
            lbRowIdx = getMaxValidRowIndex(workspace, array) + 1;
            this.config.setImportConfiguration(this.config.getImportConfiguration().toBuilder()
                    .setLbCallsetRowIdx(lbRowIdx)
                    .build());
        }
        //This sorts the list sampleNames if !useSamplesInOrder
        //Why you should use this? If you are writing multiple partitions in different machines,
        //you must have consistent ordering of samples across partitions. If file order is different
        //in different processes, then set useSamplesInOrder to false and let the sort in
        //generateSortedCallSetMap ensure consistent ordering across samples
        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB = null;
        if(this.config.sampleToReaderMapCreator() != null) //caller will create the feature readers
            callsetMappingPB =
                this.generateSortedCallSetMap(new ArrayList<>(this.config.getSampleNameToVcfPath().keySet()),
                        this.config.isUseSamplesInOrder(), lbRowIdx);
        else   //else let GenomicsDB C++ modules read the files directly
            callsetMappingPB =
                this.generateSortedCallSetMapFromNameToPathMap(this.config.getSampleNameToVcfPath(),
                        this.config.isUseSamplesInOrder(), lbRowIdx);

        String outputCallsetmapJsonFilePath = this.config.getOutputCallsetmapJsonFile();
        if (this.config.isIncrementalImport()) {
            if (outputCallsetmapJsonFilePath == null || outputCallsetmapJsonFilePath.isEmpty()) {
                throw new GenomicsDBException("Incremental import must specify callset file name");
            }
            // write out fragment names to file to help with backup/recovery
            String[] fragments = listGenomicsDBFragments(workspace);
            if(writeToFile(outputCallsetmapJsonFilePath+".fragmentlist", String.join("\n",fragments))!=0) {
                System.err.println("Warning: Could not write original fragment list for backup");
            }
            callsetMappingPB = this.mergeCallsetsForIncrementalImport(outputCallsetmapJsonFilePath,
                                       this.config.getSampleNameToVcfPath(),
                                       callsetMappingPB);
        }
        //Vid map
        String vidmapOutputFilepath = this.config.getOutputVidmapJsonFile();
        GenomicsDBVidMapProto.VidMappingPB vidMapPB;
        if (config.isIncrementalImport()) {
            if (vidmapOutputFilepath == null || vidmapOutputFilepath.isEmpty()) {
                throw new GenomicsDBException("Incremental import must specify vid file name");
            }
            vidMapPB = generateVidMapFromFile(vidmapOutputFilepath);
        } else {
            vidMapPB = generateVidMapFromMergedHeader(this.config.getMergedHeader());
        }

        //Write out callset map if needed
        if (outputCallsetmapJsonFilePath != null && !outputCallsetmapJsonFilePath.isEmpty())
            this.writeCallsetMapJSONFile(outputCallsetmapJsonFilePath, callsetMappingPB);

        //Write out vidmap if needed, don't write for incremental import
        if (!config.isIncrementalImport() && vidmapOutputFilepath != null && !vidmapOutputFilepath.isEmpty())
            this.writeVidMapJSONFile(vidmapOutputFilepath, vidMapPB);

        //Write out merged header if needed, don't write for incremental import
        String vcfHeaderOutputFilepath = this.config.getOutputVcfHeaderFile();
        if (!config.isIncrementalImport() && vcfHeaderOutputFilepath != null && !vcfHeaderOutputFilepath.isEmpty())
            this.writeVcfHeaderFile(vcfHeaderOutputFilepath, this.config.getMergedHeader());

        //Set callset map and vid map in the top level config object
        this.config.setImportConfiguration(this.config.getImportConfiguration().toBuilder()
                .setCallsetMapping(callsetMappingPB)
                .setVidMapping(vidMapPB)
                .build());

        //Create workspace folder to avoid issues with concurrency
	if (createTileDBWorkspace(workspace, false) < 0) {
	    throw new IllegalStateException(String.format("Cannot create '%s' workspace.", workspace));
        }
    }

    /**
     * Function to return the vid mapping protobuf object. 
     *
     * @return protobuf object for vid mapping
     */
    public GenomicsDBVidMapProto.VidMappingPB getProtobufVidMapping() {
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder = 
                this.config.getImportConfiguration().toBuilder();
        if (importConfigurationBuilder.hasVidMapping()) {
            return importConfigurationBuilder.getVidMapping();
        }
        else {
            return GenomicsDBVidMapProto.VidMappingPB.newBuilder().build();
        }
    }

    /**
     * Function to update vid mapping protobuf object in the top level config object.
     * Used in cases where the VCF header doesn't contain accurate information about
     * how to parse fields. For instance, allele specific annotations 
     *
     * @param vidMapPB vid mapping protobuf object to use as new
     */
    public void updateProtobufVidMapping(GenomicsDBVidMapProto.VidMappingPB vidMapPB) {
        //Write out vidmap if needed, don't write for incremental import
        String vidmapOutputFilepath = this.config.getOutputVidmapJsonFile();
        if (!config.isIncrementalImport() && vidmapOutputFilepath != null && !vidmapOutputFilepath.isEmpty())
            this.writeVidMapJSONFile(vidmapOutputFilepath, vidMapPB);

        this.config.setImportConfiguration(this.config.getImportConfiguration().toBuilder()
                .setVidMapping(vidMapPB)
                .build());
    }

    /**
     * Utility function that returns a list of ChromosomeInterval objects for
     * the column partition specified by the loader JSON file and rank/partition index
     *
     * @param loaderJSONFile path to loader JSON file
     * @param partitionIdx   rank/partition index
     * @return list of ChromosomeInterval objects for the specified partition
     * @throws ParseException when there is a bug in the JNI interface and a faulty JSON is returned
     */
    private static ArrayList<ChromosomeInterval> getChromosomeIntervalsForColumnPartition(
            final String loaderJSONFile, final int partitionIdx) throws ParseException {
        final String chromosomeIntervalsJSONString = jniGetChromosomeIntervalsForColumnPartition(loaderJSONFile, partitionIdx);
    /* JSON format
      {
        "contigs": [
           { "chr1": [ 100, 200] },
           { "chr2": [ 500, 600] }
        ]
      }
    */
        ArrayList<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        JSONParser parser = new JSONParser();
        JSONObject topObj = (JSONObject) (parser.parse(chromosomeIntervalsJSONString));
        assert topObj.containsKey("contigs");
        JSONArray listOfDictionaries = (JSONArray) (topObj.get("contigs"));
        for (Object currDictObj : listOfDictionaries) {
            JSONObject currDict = (JSONObject) currDictObj;
            assert currDict.size() == 1; //1 entry
            for (Object currEntryObj : currDict.entrySet()) {
                Map.Entry<?, ?> currEntry = (Map.Entry<?, ?>)currEntryObj;
                JSONArray currValue = (JSONArray)currEntry.getValue();
                assert currValue.size() == 2;
                chromosomeIntervals.add(new ChromosomeInterval((String)currEntry.getKey(), (Long) (currValue.get(0)),
                        (Long) (currValue.get(1))));
            }
        }
        return chromosomeIntervals;
    }

    /**
     * Utility function that returns a MultiChromosomeIterator given an FeatureReader
     * that will iterate over the VariantContext objects provided by the reader belonging
     * to the column partition specified by the loader JSON file and rank/partition index
     *
     * @param reader         AbstractFeatureReader over VariantContext objects
     * @param loaderJSONFile path to loader JSON file
     * @param partitionIdx   rank/partition index
     * @return MultiChromosomeIterator that iterates over VariantContext objects in the reader
     * belonging to the specified column partition
     * @throws IOException    when the reader's query method throws an IOException
     * @throws ParseException when there is a bug in the JNI interface and a faulty JSON is returned
     */
    public static MultiChromosomeIterator columnPartitionIterator(
            FeatureReader<VariantContext> reader,
            final String loaderJSONFile,
            final int partitionIdx) throws ParseException, IOException {
        return new MultiChromosomeIterator(reader, GenomicsDBImporter.getChromosomeIntervalsForColumnPartition(
                loaderJSONFile, partitionIdx));
    }

    private GenomicsDBImportConfiguration.ImportConfiguration addExplicitValuesToImportConfiguration(ImportConfig config) {
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder =
                config.getImportConfiguration().toBuilder();
        importConfigurationBuilder.setSegmentSize(config.getImportConfiguration().getSegmentSize())
                .setFailIfUpdating(config.getImportConfiguration().getFailIfUpdating())
                .setEnableSharedPosixfsOptimizations(config.getImportConfiguration().getEnableSharedPosixfsOptimizations())
                //TODO: making the following attributes explicit since the C++ layer is not working with the
                // protobuf object and it's defaults
                .setTreatDeletionsAsIntervals(true).setCompressTiledbArray(true).setNumCellsPerTile(1000)
                .setRowBasedPartitioning(false).setProduceTiledbArray(true).build();
        return importConfigurationBuilder.build();
    }

    /**
     * Initialize variables
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     */
    private void initialize(String loaderJSONFile, int rank) {
        mLoaderJSONFile = loaderJSONFile;
        mRank = rank;
        //Initialize C++ module
        mGenomicsDBImporterObjectHandle = jniInitializeGenomicsDBImporterObject(mLoaderJSONFile, mRank);
        if (mGenomicsDBImporterObjectHandle == 0) throw new GenomicsDBException(
                "Could not initialize GenomicsDBImporter object");
    }

    /**
     * Add a sorted VC iterator as the data source - caller must:
     * 1. Call setupGenomicsDBImporter() once all iterators are added
     * 2. Call doSingleImport()
     * 3. Done!
     *
     * @param streamName        Name of the stream being added - must be unique with respect
     *                          to this GenomicsDBImporter object
     * @param vcfHeader         VCF header for the stream
     * @param vcIterator        Iterator over VariantContext objects
     * @param bufferCapacity    Capacity of the stream buffer in bytes
     * @param streamType        BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object
     *                          which contains row index and globally unique name
     *                          can be set to null, which implies that the mapping is
     *                          stored in a callsets JSON file
     * @return returns the stream index
     */
    public int addSortedVariantContextIterator(
            final String streamName,
            final VCFHeader vcfHeader,
            Iterator<VariantContext> vcIterator,
            final long bufferCapacity,
            final VariantContextWriterBuilder.OutputType streamType,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException {
        return addBufferStream(streamName, vcfHeader, bufferCapacity, streamType, vcIterator, sampleIndexToInfo);
    }

    /**
     * Add a buffer stream or VC iterator - internal function
     *
     * @param streamName        Name of the stream being added - must be unique with respect to this
     *                          GenomicsDBImporter object
     * @param vcfHeader         VCF header for the stream
     * @param bufferCapacity    Capacity of the stream buffer in bytes
     * @param streamType        BCF_STREAM or VCF_STREAM
     * @param vcIterator        Iterator over VariantContext objects - can be null
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which
     *                          contains row index and globally unique name can be set to null,
     *                          which implies that the mapping is stored in a callsets JSON file
     * @return returns the stream index
     */
    @SuppressWarnings("unchecked")
    public int addBufferStream(final String streamName,
                               final VCFHeader vcfHeader,
                               final long bufferCapacity,
                               final VariantContextWriterBuilder.OutputType streamType,
                               Iterator<VariantContext> vcIterator,
                               final Map<Integer, SampleInfo> sampleIndexToInfo)
            throws GenomicsDBException {
        if (mIsLoaderSetupDone) throw new GenomicsDBException(
                "Cannot add buffer streams after setupGenomicsDBImporter() is called");
        //First time a buffer is added
        if (!mContainsBufferStreams)
            mContainsBufferStreams = true;
        mBufferStreamWrapperVector.add(new GenomicsDBImporterStreamWrapper(vcfHeader, bufferCapacity, streamType, vcIterator));
        int currIdx = mBufferStreamWrapperVector.size() - 1;
        SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(currIdx).getStream();
        jniAddBufferStream(mGenomicsDBImporterObjectHandle, streamName, streamType == VariantContextWriterBuilder.OutputType.BCF_STREAM,
                bufferCapacity, currStream.getBuffer(), currStream.getNumValidBytes());
        if (sampleIndexToInfo != null) {
            for (Map.Entry<Integer, SampleInfo> currEntry : sampleIndexToInfo.entrySet()) {
                JSONObject sampleJSON = new JSONObject();
                sampleJSON.put("row_idx", currEntry.getValue().getRowIdx());
                sampleJSON.put("stream_name", streamName);
                sampleJSON.put("idx_in_file", currEntry.getKey());
                mCallsetMappingJSON.put(currEntry.getValue().getName(), sampleJSON);
            }
        }
        return currIdx;
    }

    /**
     * Setup the importer after all the buffer streams are added, but before any
     * data is inserted into any stream
     * No more buffer streams can be added once setupGenomicsDBImporter() is called
     *
     * @throws IOException throws IOException if modified callsets JSON cannot be written
     */
    @SuppressWarnings("unchecked")
    public void setupGenomicsDBImporter() throws IOException {
        if (mIsLoaderSetupDone) return;
        //Callset mapping JSON - convert to string
        JSONObject topCallsetJSON = new JSONObject();
        topCallsetJSON.put("callsets", mCallsetMappingJSON);
        StringWriter stringWriter = new StringWriter();
        topCallsetJSON.writeJSONString(stringWriter);
        //Call native setupGenomicsDBImporter()
        long mMaxBufferStreamIdentifiers = jniSetupGenomicsDBLoader(mGenomicsDBImporterObjectHandle,
                stringWriter.toString());
        //Why 2* - each identifier is a pair<buffer_stream_idx, partition_idx>
        //Why +1 - the last element will contain the number of exhausted stream identifiers
        //when doSingleImport() is called
        mExhaustedBufferStreamIdentifiers = new long[2 * ((int) mMaxBufferStreamIdentifiers) + 1];
        //Set all streams to empty
        //Add all streams to mExhaustedBufferStreamIdentifiers - this way when doSingleImport is
        //called the first time all streams' data are written
        for (int i = 0, idx = 0; i < mBufferStreamWrapperVector.size(); ++i, idx += 2) {
            SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(i).getStream();
            currStream.setNumValidBytes(0);
            mExhaustedBufferStreamIdentifiers[idx] = i;
            mExhaustedBufferStreamIdentifiers[idx + 1] = 0;
        }
        //Set number of exhausted buffer streams - all streams are exhausted the first time
        mNumExhaustedBufferStreams = mBufferStreamWrapperVector.size();
        mExhaustedBufferStreamIdentifiers[(int) (2 * mMaxBufferStreamIdentifiers)] = mNumExhaustedBufferStreams;
        mIsLoaderSetupDone = true;
    }

    /**
     * Write VariantContext object to stream - may fail if the buffer is full
     * It's the caller's responsibility keep track of the VC object that's not written
     *
     * @param vc        VariantContext object
     * @param streamIdx index of the stream returned by the addBufferStream() call
     * @return true if the vc object was written successfully, false otherwise
     */
    public boolean add(VariantContext vc, final int streamIdx) throws GenomicsDBException, RuntimeIOException {
        if (streamIdx < 0 || streamIdx >= mBufferStreamWrapperVector.size()) throw new GenomicsDBException(
                "Invalid stream idx " + Integer.toString(streamIdx) + " must be between [0-"
                        + Long.toString(mBufferStreamWrapperVector.size() - 1) + "]");
        if (!mIsLoaderSetupDone) throw new GenomicsDBException("Cannot add VariantContext objects to streams before " +
                "calling setupGenomicsDBImporter()");
        GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(streamIdx);
        currWrapper.getVcWriter().add(vc);
        SilentByteBufferStream currStream = currWrapper.getStream();
        //at least one record already existed in the buffer
        if (currStream.overflow() && currStream.getMarker() > 0) {
            //Set num valid bytes to marker - marker points to location
            //after the last valid serialized vc in the buffer
            currStream.setNumValidBytes(currStream.getMarker());
            return false;
        }
        //the first record to be added to the buffer is too large, resize buffer
        while (currStream.overflow()) {
            currStream.resize(2 * currStream.size() + 1);
            currStream.setNumValidBytes(0);
            currStream.setOverflow(false);
            currWrapper.getVcWriter().add(vc);
        }
        //Update marker - marker points to location after the last valid serialized vc in the buffer
        currStream.setMarker(currStream.getNumValidBytes());
        return true;
    }

    /**
     * Only to be used in cases where iterator of VariantContext are not used. The data is written to buffers
     * directly after which this function is called. See TestBufferStreamGenomicsDBImporter.java for an example
     *
     * @return true if the import process is done
     * @throws IOException if the import fails
     */
    public boolean doSingleImport() throws IOException {
        if (mDone) return true;
        if (!mIsLoaderSetupDone) setupGenomicsDBImporter();
        boolean allExhaustedStreamsHaveIterators = true;
        while (!mDone && allExhaustedStreamsHaveIterators) {
            //Write data from buffer streams exhausted in the previous round into GenomicsDB
            for (int i = 0, idx = 0; i < mNumExhaustedBufferStreams; ++i, idx += 2) {
                int bufferStreamIdx = (int) mExhaustedBufferStreamIdentifiers[idx];
                GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(bufferStreamIdx);
                //If iterator is provided, get data from iterator
                if (currWrapper.hasIterator()) {
                    while (currWrapper.getCurrentVC() != null) {
                        boolean added = add(currWrapper.getCurrentVC(), bufferStreamIdx);
                        if (added) currWrapper.next();
                        else break; //buffer full
                    }
                }
                SilentByteBufferStream currStream = currWrapper.getStream();
                jniWriteDataToBufferStream(mGenomicsDBImporterObjectHandle, bufferStreamIdx, 0, currStream.getBuffer(),
                        currStream.getNumValidBytes());
            }
            mDone = jniImportBatch(mGenomicsDBImporterObjectHandle, mExhaustedBufferStreamIdentifiers);
            mNumExhaustedBufferStreams = mExhaustedBufferStreamIdentifiers[mExhaustedBufferStreamIdentifiers.length - 1];
            //Reset markers, numValidBytesInBuffer and overflow flag for the exhausted streams
            for (long i = 0, idx = 0; i < mNumExhaustedBufferStreams; ++i, idx += 2) {
                int bufferStreamIdx = (int) mExhaustedBufferStreamIdentifiers[(int) idx];
                GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(bufferStreamIdx);
                if (!currWrapper.hasIterator()) allExhaustedStreamsHaveIterators = false;
                SilentByteBufferStream currStream = currWrapper.getStream();
                currStream.setOverflow(false);
                currStream.setMarker(0);
                currStream.setNumValidBytes(0);
            }
            if (mDone) {
                mGenomicsDBImporterObjectHandle = 0;
                mContainsBufferStreams = false;
                mIsLoaderSetupDone = false;
            }
        }

        return mDone;
    }

    /**
     * Import multiple chromosome interval
     *
     * @throws InterruptedException when there is an exception in any of the threads in the stream
     */
    public void executeImport() throws InterruptedException {
        executeImport(0);
    }

    /**
     * Import multiple chromosome interval
     * @param numThreads number of threads used to import partitions
     */
    public void executeImport(final int numThreads) {
        final int batchSize = this.config.getBatchSize();
        final int sampleCount = this.config.getSampleNameToVcfPath().size();
        final int updatedBatchSize = (batchSize <= 0) ? sampleCount : Math.min(batchSize, sampleCount);
        final int numberPartitions = this.config.getImportConfiguration().getColumnPartitionsList().size();
        final ExecutorService executor = numThreads == 0 ? ForkJoinPool.commonPool() : Executors.newFixedThreadPool(numThreads);
        final boolean performConsolidation = this.config.getImportConfiguration().getConsolidateTiledbArrayAfterLoad();

        //Set size_per_column_partition once
        this.config.setImportConfiguration(this.config.getImportConfiguration().toBuilder()
                .setSizePerColumnPartition(this.config.getImportConfiguration().getSizePerColumnPartition()
                        * updatedBatchSize).build());

        //Iterate over sorted sample list in batches
        iterateOverSamplesInBatches(sampleCount, updatedBatchSize, numberPartitions, executor, performConsolidation);
        executor.shutdown();
        deleteBackupFiles();
        mDone = true;
    }

    private void deleteBackupFiles() {
        if (this.config.isIncrementalImport()) {
            String outputCallsetmapJsonFilePath = this.config.getOutputCallsetmapJsonFile();
            if (outputCallsetmapJsonFilePath == null || outputCallsetmapJsonFilePath.isEmpty()) {
                throw new GenomicsDBException("Incremental import must specify callset file name");
            }
            if(deleteFile(outputCallsetmapJsonFilePath+".inc.backup")!=0) {
                System.err.println("Warning: Could not delete backup callset file"); 
            }
            if(deleteFile(outputCallsetmapJsonFilePath+".fragmentlist")!=0) {
                System.err.println("Warning: Could not delete fragment list backup file"); 
            }
        }
    }

    private void iterateOverSamplesInBatches(final int sampleCount, final int updatedBatchSize, final int numberPartitions,
                                             final ExecutorService executor, final boolean performConsolidation) {
        final boolean failIfUpdating = this.config.getImportConfiguration().getFailIfUpdating();
        final int totalBatchCount = (sampleCount + updatedBatchSize - 1)/updatedBatchSize; //ceil
        BatchCompletionCallbackFunctionArgument callbackFunctionArgument = new BatchCompletionCallbackFunctionArgument(0, totalBatchCount);
        int origLbRowIdx = (this.config.getImportConfiguration().hasLbCallsetRowIdx())
                           ? ((int)this.config.getImportConfiguration().getLbCallsetRowIdx())
                           : (0);
        for (int i = origLbRowIdx, batchCount = 1; i < sampleCount+origLbRowIdx; i += updatedBatchSize, ++batchCount) {
            final int index = i;
            IntStream.range(0, numberPartitions).forEach(rank -> updateConfigPartitionsAndLbUb(this.config, index, rank));
            GenomicsDBImportConfiguration.ImportConfiguration updatedConfig =
                    this.config.getImportConfiguration().toBuilder().setConsolidateTiledbArrayAfterLoad(
                            ((i + updatedBatchSize) >= (sampleCount+origLbRowIdx)) && performConsolidation)
                    .setFailIfUpdating(failIfUpdating && i == 0) //fail if updating should be set to true iff user sets it and this is the first batch
                    .build();
            this.config.setImportConfiguration(updatedConfig);
            // sending index-origLbRowIdx because iterativeOverChromosomeIntervals uses that to index sampleMap
            // and in incremental import case that will only have the incremental callset
            // whereas updateConfigPartitionsAndLbUb sets config that will index into merged callsetPB
            List<CompletableFuture<Boolean>> futures = iterateOverChromosomeIntervals(updatedBatchSize, numberPartitions, executor, index-origLbRowIdx);
            List<Boolean> result = futures.stream().map(CompletableFuture::join).collect(Collectors.toList());

            if (result.contains(false)) {
                executor.shutdown();
                throw new IllegalStateException("There was an unhandled exception during chromosome interval import.");
            }
            if(this.config.getFunctionToCallOnBatchCompletion() != null) {
                callbackFunctionArgument.batchCount = batchCount;
                this.config.getFunctionToCallOnBatchCompletion().apply(callbackFunctionArgument);
            }
        }
    }

    private List<CompletableFuture<Boolean>> iterateOverChromosomeIntervals(final int updatedBatchSize, final int numberPartitions,
                                                                            final ExecutorService executor, final int index) {
        return IntStream.range(0, numberPartitions).mapToObj(rank ->
                CompletableFuture.supplyAsync(() -> {
                    try {
                        final Map<String, FeatureReader<VariantContext>> sampleToReaderMap =
                                (this.config.sampleToReaderMapCreator() != null)
                                    ? this.config.sampleToReaderMapCreator().apply(
                                        this.config.getSampleNameToVcfPath(), updatedBatchSize, index)
                                    : null;
                        GenomicsDBImporter importer = new GenomicsDBImporter(this.config, sampleToReaderMap, rank);
                        Boolean result = importer.doSingleImport();
                        if(sampleToReaderMap != null) {
                            sampleToReaderMap.values().forEach(v -> {
                                try {
                                    v.close();
                                } catch (IOException e) {
                                    throw new IllegalStateException(e);
                                }
                            });
                        }
                        return result;
                    } catch (IOException ex) {
                        throw new IllegalStateException("There was an unhandled exception during chromosome interval import.", ex);
                    }
                }, executor)
        ).collect(Collectors.toList());
    }

    private void updateConfigPartitionsAndLbUb(ImportConfig importConfig, final int index, final int rank) {
        GenomicsDBImportConfiguration.Partition partition = importConfig.getImportConfiguration().getColumnPartitions(rank);
        if (partition.hasGenerateArrayNameFromPartitionBounds()) {
            String arrayName = Integer.toString(rank);
            if (partition.getBegin().hasContigPosition()) {
                final String chromosomeName = partition.getBegin().getContigPosition().getContig();
                final int chromosomeStart = (int) partition.getBegin().getContigPosition().getPosition();
                final int chromosomeEnd = (int) partition.getEnd().getContigPosition().getPosition();

                arrayName = String.format(Constants.CHROMOSOME_INTERVAL_FOLDER, chromosomeName, chromosomeStart, chromosomeEnd);
            }
            partition = partition.toBuilder().setArrayName(arrayName).build();
        }

        GenomicsDBImportConfiguration.ImportConfiguration importConfiguration = importConfig
                .getImportConfiguration().toBuilder()
                .setLbCallsetRowIdx((long) index)
                .setUbCallsetRowIdx((long) (index + importConfig.getBatchSize() - 1))
                .setColumnPartitions(rank, partition).build();
        importConfig.setImportConfiguration(importConfiguration);
    }

    /**
     * @return get number of buffer streams for which new data must be supplied
     */
    public long getNumExhaustedBufferStreams() {
        return mNumExhaustedBufferStreams;
    }

    /**
     * Get buffer stream index of i-th exhausted stream
     * There are mNumExhaustedBufferStreams and the caller must provide data for streams
     * with indexes getExhaustedBufferStreamIndex(0), getExhaustedBufferStreamIndex(1),...,
     * getExhaustedBufferStreamIndex(mNumExhaustedBufferStreams-1)
     *
     * @param i i-th exhausted buffer stream
     * @return the buffer stream index of the i-th exhausted stream
     */
    public int getExhaustedBufferStreamIndex(final long i) {
        assert i < mNumExhaustedBufferStreams && i >= 0;
        //why 2* - exhausted buffer stream identifier is a pair<stream_idx, partition_idx>
        return (int) mExhaustedBufferStreamIdentifiers[2 * ((int) i)];
    }

    /**
     * Is the import process completed
     *
     * @return true if complete, false otherwise
     */
    public boolean isDone() {
        return mDone;
    }

    /**
     * Utility function that returns a MultiChromosomeIterator given an FeatureReader
     * that will iterate over the VariantContext objects provided by the reader belonging
     * to the column partition specified by this object's loader JSON file and rank/partition index
     *
     * @param reader   AbstractFeatureReader over VariantContext objects
     * @return MultiChromosomeIterator that iterates over VariantContext objects in the reader
     * belonging to the specified column partition
     * @throws IOException    when the reader's query method throws an IOException
     * @throws ParseException when there is a bug in the JNI interface and a faulty JSON is returned
     */
    public MultiChromosomeIterator columnPartitionIterator(
            FeatureReader<VariantContext> reader) throws ParseException, IOException {
        return GenomicsDBImporter.columnPartitionIterator(reader, mLoaderJSONFile, mRank);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     *
     */
    public void write() throws GenomicsDBException {
        write(mLoaderJSONFile, mRank);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     *
     * @param rank     Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public void write(final int rank) throws GenomicsDBException {
        write(mLoaderJSONFile, rank);
    }

    /**
     * Write to TileDB/GenomicsDB
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public void write(final String loaderJSONFile, final int rank)
            throws GenomicsDBException {
        mDone = false;
        if (loaderJSONFile == null) throw new GenomicsDBException("Loader JSON file not specified");
        if (mContainsBufferStreams)
            throw new GenomicsDBException("Cannot call write() functions if buffer streams are added");
        int status = jniGenomicsDBImporter(loaderJSONFile, rank);
        if (status != 0) throw new GenomicsDBException("GenomicsDBImporter write failed for loader JSON: "
                + loaderJSONFile + " rank: " + rank);
        mDone = true;
    }

    /*
     * Create and return Partition protobuf object using start and workspace attributes
     */
    private GenomicsDBImportConfiguration.Partition.Builder initializePartitionPB(final long start,
            final String workspace) {
        GenomicsDBImportConfiguration.Partition.Builder partitionBuilder =
                GenomicsDBImportConfiguration.Partition.newBuilder();
        Coordinates.GenomicsDBColumn.Builder columnBuilder =
                Coordinates.GenomicsDBColumn.newBuilder().setTiledbColumn(start);
        partitionBuilder.setBegin(columnBuilder.build());
        partitionBuilder.setWorkspace(workspace);
        partitionBuilder.setGenerateArrayNameFromPartitionBounds(true);
        return partitionBuilder;
    }

    /*
     * Returns Vid contigs that match loader column partitions
     */
    private List<GenomicsDBVidMapProto.Chromosome> getContigsFromPartitions(
            GenomicsDBImportConfiguration.ImportConfiguration.Builder builder) {
        assert builder.hasVidMapping() : "ImportConfiguration must specify VidMapping";
        GenomicsDBVidMapProto.VidMappingPB vidPB = builder.getVidMapping();
        assert builder.getColumnPartitions(0).getBegin().hasContigPosition() : 
                "Partitions must be specified using ContigPositions";

        Set<String> contigNames = builder.getColumnPartitionsList().stream()
                .map(x -> x.getBegin().getContigPosition().getContig())
                .collect(Collectors.toSet());

        return vidPB.getContigsList().stream()
                .filter(x -> contigNames.contains(x.getName()))
                .collect(Collectors.toList());
    }

    /*
     * Returns Vid protobuf object that only has contigs specified in contigList
     */
    private GenomicsDBVidMapProto.VidMappingPB regenerateVidUsingContigList(
            GenomicsDBImportConfiguration.ImportConfiguration.Builder builder,
             final List<GenomicsDBVidMapProto.Chromosome> contigList) {
        assert builder.hasVidMapping() : "ImportConfiguration must specify VidMapping";
        GenomicsDBVidMapProto.VidMappingPB.Builder vidBuilder = 
                builder.getVidMapping().toBuilder();

        // only need to regenerate vid if interval list (read: map)
        // doesn't include all the contigs
        if (vidBuilder.getContigsCount() != contigList.size()) {
            vidBuilder.clearContigs();
            long columnOffset = 0;
            for (GenomicsDBVidMapProto.Chromosome entry : contigList) {
                GenomicsDBVidMapProto.Chromosome.Builder chr = GenomicsDBVidMapProto.Chromosome.newBuilder();
                chr.setName(entry.getName()).setLength(entry.getLength()).setTiledbColumnOffset(columnOffset);
                vidBuilder.addContigs(chr.build());
                columnOffset += entry.getLength();
            }
            GenomicsDBVidMapProto.VidMappingPB vidPB = vidBuilder.build();
            builder.setVidMapping(vidPB);

            // if we wrote out a vid file earlier, amend it
            String vidmapOutputFilepath = this.config.getOutputVidmapJsonFile();
            if (!config.isIncrementalImport() && vidmapOutputFilepath != null && !vidmapOutputFilepath.isEmpty())
                this.writeVidMapJSONFile(vidmapOutputFilepath, vidPB);

            return vidPB;
        }
        return vidBuilder.build();
    }

    /**
     * Coalesce contigs into fewer GenomicsDB partitions
     *
     * @param partitions Approximate number of partitions to coalesce into
     * @throws GenomicsDBException Coalescing contigs into partitions fails
     */
    public void coalesceContigsIntoNumPartitions(final int partitions)
            throws GenomicsDBException {
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder =
                this.config.getImportConfiguration().toBuilder();
        if (importConfigurationBuilder.hasVidMapping()) {
            GenomicsDBVidMapProto.VidMappingPB vidPB = importConfigurationBuilder.getVidMapping();
            // we assume here that contigs in vid are sorted in ascending tiledb column offset order
            // this is true for vids generated by GenomicsDBImporter...I suppose it could be false
            // for hand written vids...but we'll throw an exception later if not true
            final String workspace = importConfigurationBuilder.getColumnPartitions(0).getWorkspace();
            List<GenomicsDBVidMapProto.Chromosome> contigList = getContigsFromPartitions(importConfigurationBuilder);
            vidPB = regenerateVidUsingContigList(importConfigurationBuilder, contigList);
            final int numContigs = vidPB.getContigsCount();
            final long genomeLength = vidPB.getContigs(numContigs-1).getTiledbColumnOffset() +
                    vidPB.getContigs(numContigs-1).getLength() + 1;
            importConfigurationBuilder.clearColumnPartitions();

            // initialize begin for first partition to 0
            GenomicsDBImportConfiguration.Partition.Builder partitionBuilder =
                    initializePartitionPB(0, workspace);
            long partitionSum = 0;
            long previousEnd = -1;
            final long goalSize = genomeLength / partitions;
            for (GenomicsDBVidMapProto.Chromosome contig: vidPB.getContigsList()) {
                if (previousEnd + 1 != contig.getTiledbColumnOffset()) {
                    throw new GenomicsDBException("GenomicsDBImporter expects contigs in " +
                            "vid to be contiguous and be sorted by tiledb_column_offset. " + 
                            "Contig "+contig.getName()+" has tiledb_column_offset "+contig.getTiledbColumnOffset()+
                            " while previous contig ended at "+previousEnd);
                }
                partitionSum += contig.getLength();
                previousEnd = contig.getTiledbColumnOffset() + contig.getLength() - 1;
                if (partitionSum >= goalSize) {
                    Coordinates.GenomicsDBColumn.Builder endBuilder =
                        Coordinates.GenomicsDBColumn.newBuilder().setTiledbColumn(previousEnd);
                    partitionBuilder.setEnd(endBuilder.build());
                    importConfigurationBuilder.addColumnPartitions(partitionBuilder.build());
                    partitionBuilder = initializePartitionPB(previousEnd+1, workspace);
                    partitionSum = 0;
                }
            }

            // finish up the final partition
            Coordinates.GenomicsDBColumn.Builder endBuilder =
                    Coordinates.GenomicsDBColumn.newBuilder().setTiledbColumn(genomeLength-1);
            partitionBuilder.setEnd(endBuilder.build());
            importConfigurationBuilder.addColumnPartitions(partitionBuilder.build());
            this.config.setImportConfiguration(importConfigurationBuilder.build());
        }
        else {
            throw new GenomicsDBException("GenomicsDBImporter needs vid mapping to coalesce contigs into partitions");
        }
    }
}
