/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
 * Copyright (c) 2018-2019 Omics Data Automation, Inc.
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

package org.genomicsdb.reader;

import htsjdk.tribble.*;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBVidMapProto;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.importer.extensions.JsonFileExtensions;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toList;
import static org.genomicsdb.Constants.CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX;
import static org.genomicsdb.GenomicsDBUtils.listGenomicsDBArrays;
import static org.genomicsdb.GenomicsDBUtils.getArrayColumnBounds;
import static com.googlecode.protobuf.format.JsonFormat.ParseException;

/**
 * A reader for GenomicsDB that implements {@link htsjdk.tribble.FeatureReader}
 * Currently, the reader only return {@link htsjdk.variant.variantcontext.VariantContext}
 */
public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T>, JsonFileExtensions {
    private String loaderJSONFile;
    private GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration;
    private FeatureCodec<T, SOURCE> codec;
    private FeatureCodecHeader featureCodecHeader;
    private List<String> sequenceNames;
    private Map<String, Coordinates.ContigInterval> intervalsPerArray = new HashMap<>();
    /**
     * Constructor
     *
     * @param exportConfiguration query parameters
     * @param codec               FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                            and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param loaderJSONFile      GenomicsDB loader JSON configuration file
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration,
                                   final FeatureCodec<T, SOURCE> codec,
                                   final Optional<String> loaderJSONFile) throws IOException {
        GenomicsDBExportConfiguration.ExportConfiguration.Builder builder = 
                exportConfiguration.toBuilder();
        if(exportConfiguration.getQueryColumnRangesCount() == 0
                && exportConfiguration.getQueryContigIntervalsCount() == 0
                && exportConfiguration.getQueryRowRangesCount() == 0
	        && exportConfiguration.getQuerySampleNamesCount() == 0) {
            builder.setScanFull(true);
        }
        this.exportConfiguration = builder.build();
        this.codec = codec;
        this.loaderJSONFile = loaderJSONFile.orElse("");
        List<String> chromosomeIntervalArrays = this.exportConfiguration.hasArrayName() ? new ArrayList<String>() {{
            add(exportConfiguration.getArrayName());
        }} : getArrayListFromWorkspace(exportConfiguration.getWorkspace(), Optional.empty());
        if (chromosomeIntervalArrays == null || chromosomeIntervalArrays.size() < 1)
            throw new IllegalStateException("There is no genome data stored in the database");
        generateHeadersForQuery(chromosomeIntervalArrays.get(0));
    }

    /**
     * Return the VCF header of the combined gVCF stream
     *
     * @return the VCF header of the combined gVCF stream
     */
    public Object getHeader() {
        return this.featureCodecHeader.getHeaderValue();
    }

    /**
     * Return the list of contigs in the combined VCF header
     *
     * @return list of strings of the contig names
     */
    public List<String> getSequenceNames() {
        return this.sequenceNames;
    }

    public void close() throws IOException {
    }

    /**
     * Return an iterator over {@link htsjdk.variant.variantcontext.VariantContext}
     * objects for the specified TileDB array and query configuration
     *
     * @return iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects
     */
    public CloseableTribbleIterator<T> iterator() throws IOException {
        if (this.exportConfiguration.hasArrayName()) {
            return new GenomicsDBFeatureIterator<>(this.loaderJSONFile,
                    this.exportConfiguration, Optional.empty(), 
                    this.featureCodecHeader, this.codec);
        }
        else {
            List<String> chromosomeIntervalArraysPaths = 
                    resolveChromosomeArrayFolderList(Optional.empty());
            return new GenomicsDBFeatureIterator<>(this.loaderJSONFile,
                    this.exportConfiguration, Optional.of(chromosomeIntervalArraysPaths),
                    this.featureCodecHeader, this.codec);
        }
    }

    /**
     * Return an iterator over {@link htsjdk.variant.variantcontext.VariantContext}
     * objects for the specified TileDB array and queried position
     *
     * @param chr   contig name
     * @param start start position (1-based)
     * @param end   end position, inclusive (1-based)
     * @return iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects
     */
    public CloseableTribbleIterator<T> query(final String chr, final int start, final int end) throws IOException {
        if (this.exportConfiguration.hasArrayName()) {
            return new GenomicsDBFeatureIterator<>(this.loaderJSONFile,
                    this.exportConfiguration, Optional.empty(), 
                    this.featureCodecHeader, this.codec, chr, OptionalInt.of(start), OptionalInt.of(end));
        }
        else {
            Optional<Coordinates.ContigInterval> contigInterval = Optional.of(
                    Coordinates.ContigInterval.newBuilder().setContig(chr).setBegin(start).setEnd(end).build());
            List<String> chromosomeIntervalArraysPaths = 
                    resolveChromosomeArrayFolderList(contigInterval);
            return new GenomicsDBFeatureIterator<>(this.loaderJSONFile,
                    this.exportConfiguration, Optional.of(chromosomeIntervalArraysPaths),
                    this.featureCodecHeader, this.codec, chr, OptionalInt.of(start), OptionalInt.of(end));
        }
    }


    private List<String> resolveChromosomeArrayFolderList(final Optional<Coordinates.ContigInterval> chromosome) {
        List<String> chromosomeIntervalArraysNames = getArrayListFromWorkspace(exportConfiguration.getWorkspace(), chromosome);
        chromosomeIntervalArraysNames.sort(new ChrArrayFolderComparator());
        return chromosomeIntervalArraysNames;
    }


    private boolean checkIfContigOverlapsArray(String workspace, String array, 
            Coordinates.ContigInterval contigInterval) {
        String[] ref = array.split(CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX);
        if (ref.length == 3) {
            return contigInterval.getContig().equals(ref[0]) && (contigInterval.getBegin() <= Integer.parseInt(ref[2])
                    && contigInterval.getEnd() >= Integer.parseInt(ref[1]));
        }
        else {
            try {
                GenomicsDBVidMapProto.VidMappingPB vidPB = this.exportConfiguration.hasVidMapping() ?
                        this.exportConfiguration.getVidMapping() :
                        generateVidMapFromFile(this.exportConfiguration.getVidMappingFile());
                long[] bounds = getArrayColumnBounds(workspace, array);
                // here we only check if contig starts within column bounds
                // this works because currently the contig coalescing respects contig boundaries
                // that is, we are guaranteed that entire contigs reside in the same array
                // if that changes we need to get smarter here
                return checkVidIfContigStartsWithinColumnBounds(vidPB, bounds, contigInterval);
            } catch (ParseException e) {
                throw new RuntimeException("Error parsing vid from file");
            }
        }
    }


    private void generateHeadersForQuery(final String randomExistingArrayName) throws IOException {
        GenomicsDBExportConfiguration.ExportConfiguration.Builder fullExportConfigurationBuilder =
                GenomicsDBExportConfiguration.ExportConfiguration.newBuilder(this.exportConfiguration)
                        .setArrayName(randomExistingArrayName);
        if(this.exportConfiguration.getQueryColumnRangesCount() == 0
	        && this.exportConfiguration.getQueryContigIntervalsCount() == 0
                && this.exportConfiguration.getQueryRowRangesCount() == 0
	        && this.exportConfiguration.getQuerySampleNamesCount() == 0) {
            fullExportConfigurationBuilder.setScanFull(true);
        }
        GenomicsDBQueryStream gdbStream = new GenomicsDBQueryStream(
                this.loaderJSONFile, fullExportConfigurationBuilder.build(),
                this.codec instanceof BCF2Codec, true);
        SOURCE source = this.codec.makeSourceFromStream(gdbStream);
        this.featureCodecHeader = this.codec.readHeader(source);
        //Store sequence names
        VCFHeader vcfHeader = (VCFHeader) (this.featureCodecHeader.getHeaderValue());
        this.sequenceNames = new ArrayList<>(vcfHeader.getContigLines().size());
        for (final VCFContigHeaderLine contigHeaderLine : vcfHeader.getContigLines())
            this.sequenceNames.add(contigHeaderLine.getID());
        gdbStream.close();
    }

    private List<String> getArrayListFromWorkspace(final String workspace_str, Optional<Coordinates.ContigInterval> chromosome) {
	List<String> folders = Arrays.asList(listGenomicsDBArrays(workspace_str));
        return chromosome.map(contigInterval -> folders.stream().filter(name -> 
            checkIfContigOverlapsArray(workspace_str, name, contigInterval)
        ).collect(toList())).orElse(folders);
    }
}
