/*
 * The MIT License (MIT)
 * Copyright (c) 2018-2023 Omics Data Automation, Inc.
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



package org.genomicsdb.importer.extensions;

import com.google.protobuf.util.JsonFormat;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.genomicsdb.model.GenomicsDBCallsetsMapProto;
import org.genomicsdb.GenomicsDBUtils;
import org.genomicsdb.exception.GenomicsDBException;

import java.io.IOException;
import java.net.URI;
import java.util.*;

public interface CallSetMapExtensions {
    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     *
     * @param sampleNames               list of sample names
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(final List<String> sampleNames,
                                                                                 final boolean useSamplesInOrderProvided) {
        return generateSortedCallSetMap(sampleNames, useSamplesInOrderProvided, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes. This method is synchronized
     * to block multiple invocations (if by any chance) disturb the order in
     * which TileDB row indexes are generated
     *
     * @param sampleNames               list of sample names
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(List<String> sampleNames,
                                                                                 final boolean useSamplesInOrderProvided,
                                                                                 final long lbRowIdx) {

        if (!useSamplesInOrderProvided)
            Collections.sort(sampleNames);
        LinkedHashMap<String, String> sampleNameToStreamName = new LinkedHashMap<String, String>();
        for(final String sample : sampleNames)
            sampleNameToStreamName.put(sample, getStreamNameFromSampleName(sample));
        return generateSortedCallSetMap(sampleNameToStreamName, lbRowIdx);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes. This method is synchronized
     * to block multiple invocations (if by any chance) disturb the order in
     * which TileDB row indexes are generated
     *
     * @param inputSampleNameToPath    map from sample name to VCF/BCF file path
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMapFromNameToPathMap(
            final Map<String, URI> inputSampleNameToPath,
            final boolean useSamplesInOrderProvided,
            final long lbRowIdx) {

        LinkedHashMap<String, String> sampleNameToStreamName = new LinkedHashMap<String, String>();
        for(Map.Entry<String, URI> currEntry : inputSampleNameToPath.entrySet())
            sampleNameToStreamName.put(currEntry.getKey(), currEntry.getValue().toString());
        return generateSortedCallSetMap(sampleNameToStreamName, useSamplesInOrderProvided, lbRowIdx);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes. This method is synchronized
     * to block multiple invocations (if by any chance) disturb the order in
     * which TileDB row indexes are generated
     *
     * @param inputSampleNameToStreamName    map from sample name to VCF/BCF file path
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(final Map<String, String> inputSampleNameToStreamName,
                                                                                 final boolean useSamplesInOrderProvided,
                                                                                 final long lbRowIdx) {

        TreeMap<String, String> sortedMap = new TreeMap<String, String>();
        if (!useSamplesInOrderProvided)
            for(Map.Entry<String, String> currEntry : inputSampleNameToStreamName.entrySet())
                sortedMap.put(currEntry.getKey(), currEntry.getValue());
        LinkedHashMap<String, String> sampleNameToStreamName = new LinkedHashMap<String, String>();
        Iterator<Map.Entry<String, String>> iter = useSamplesInOrderProvided
            ? inputSampleNameToStreamName.entrySet().iterator()
            : sortedMap.entrySet().iterator();
        while(iter.hasNext()) {
            Map.Entry<String, String> currEntry = iter.next();
            sampleNameToStreamName.put(currEntry.getKey(), currEntry.getValue());
        }
        return generateSortedCallSetMap(sampleNameToStreamName, lbRowIdx);
    }

    /**
     * Creates CallSets Protobuf structure for given map
     *
     * @param sampleNameToStreamName    map from sample name to VCF/BCF file path
     *                                  use the order they appear in
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final LinkedHashMap<String, String> sampleNameToStreamName,
            final long lbRowIdx) {

        GenomicsDBCallsetsMapProto.CallsetMappingPB.Builder callsetMapBuilder =
                GenomicsDBCallsetsMapProto.CallsetMappingPB.newBuilder();

        long tileDBRowIndex = lbRowIdx;

        for (Map.Entry<String,String> currEntry : sampleNameToStreamName.entrySet()) {
            GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.Builder idMapBuilder =
                    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.newBuilder();

            idMapBuilder.setSampleName(currEntry.getKey()).setRowIdx(tileDBRowIndex++).setIdxInFile(0)
                    .setStreamName(currEntry.getValue());

            GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sampleIDToTileDBIDMap =
                    idMapBuilder.build();

            callsetMapBuilder.addCallsets(sampleIDToTileDBIDMap);
        }

        return callsetMapBuilder.build();
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     * <p>
     * Assume one sample per input GVCF file
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param validateSampleToReaderMap If True, check i) whether sample names are consistent
     *                                  with headers and ii) feature readers are valid
     *                                  in sampleToReaderMap
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
            final boolean validateSampleToReaderMap,
            final boolean useSamplesInOrderProvided) {
        return generateSortedCallSetMap(sampleToReaderMap, useSamplesInOrderProvided, validateSampleToReaderMap, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     * <p>
     * Assume one sample per input GVCF file
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param validateSampleMap         Check i) whether sample names are consistent
     *                                  with headers and ii) feature readers are valid
     *                                  in sampleToReaderMap
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
            final boolean useSamplesInOrderProvided,
            final boolean validateSampleMap,
            final long lbRowIdx) {
        if (!validateSampleMap) return generateSortedCallSetMap(new ArrayList<>(sampleToReaderMap.keySet()),
                useSamplesInOrderProvided, lbRowIdx);

        List<String> listOfSampleNames = new ArrayList<>(sampleToReaderMap.size());
        for (Map.Entry<String, FeatureReader<VariantContext>> mapObject : sampleToReaderMap.entrySet()) {
            String sampleName = mapObject.getKey();
            FeatureReader<VariantContext> featureReader = mapObject.getValue();

            if (featureReader == null)
                throw new IllegalArgumentException("Null FeatureReader found for sample: " + sampleName);

            VCFHeader header = (VCFHeader) featureReader.getHeader();
            List<String> sampleNamesInHeader = header.getSampleNamesInOrder();
            if (sampleNamesInHeader.size() > 1) {
                StringBuilder messageBuilder = new StringBuilder("Multiple samples ");
                for (String name : sampleNamesInHeader)
                    messageBuilder.append(name).append(" ");
                messageBuilder.append(" appear in header for ").append(sampleName);
                throw new IllegalArgumentException(messageBuilder.toString());
            } else {
                if (!sampleName.equals(sampleNamesInHeader.get(0))) {
                    System.err.println("Sample " + header + " does not match with " + sampleNamesInHeader.get(0) +
                            " in header");
                }
            }
            listOfSampleNames.add(sampleName);
        }
        return generateSortedCallSetMap(listOfSampleNames, useSamplesInOrderProvided, lbRowIdx);
    }

    /**
     * Merges incremental import's callsets with existing callsets.
     * @param callsetMapJSONFilePath path to existing callset map file
     * @param inputSampleNameToPath    map from sample name to VCF/BCF file path
     * @param newCallsetMapPB callset mapping protobuf for new callsets
     * @return merged callsets for callsets to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB mergeCallsetsForIncrementalImport(
            final String callsetMapJSONFilePath,
            final Map<String, URI> inputSampleNameToPath,
            final GenomicsDBCallsetsMapProto.CallsetMappingPB newCallsetMapPB)
            throws GenomicsDBException {
        String existingCallsetsJSON = GenomicsDBUtils.readEntireFile(callsetMapJSONFilePath);
        checkDuplicateCallsetsForIncrementalImport(existingCallsetsJSON, inputSampleNameToPath);
        GenomicsDBCallsetsMapProto.CallsetMappingPB.Builder callsetMapBuilder = GenomicsDBCallsetsMapProto.CallsetMappingPB.newBuilder();
        try {
            JsonFormat.parser().merge(existingCallsetsJSON, callsetMapBuilder);
            callsetMapBuilder.mergeFrom(newCallsetMapPB);
        } catch (IOException e) {
          throw new GenomicsDBException(String.format("Could not merge incremental import's callsets with existing callsets : %s", e.getMessage()));
        }
        return callsetMapBuilder.build();
    }

    /**
     * Checks for duplicates between existing and incremental callsets
     * @param existingCallsetsJSON json string with existing callset
     * @param inputSampleNameToPath    map from sample name to VCF/BCF file path
     */
    default void checkDuplicateCallsetsForIncrementalImport(
            final String existingCallsetsJSON,
            final Map<String, URI> inputSampleNameToPath) 
            throws GenomicsDBException {
        // create PB from existing json string to iterate through it
        // maybe better to use some json library here since we just need to 
        // iterate through the existing callset
        GenomicsDBCallsetsMapProto.CallsetMappingPB.Builder callsetMapBuilder =
                GenomicsDBCallsetsMapProto.CallsetMappingPB.newBuilder();
        try {
          JsonFormat.parser().merge(existingCallsetsJSON, callsetMapBuilder);
        } catch (IOException e) {
          throw new GenomicsDBException(String.format("Could not check for duplicates between existing and incremental callsets : %s", e.getMessage()));
        }
        GenomicsDBCallsetsMapProto.CallsetMappingPB existingCallsetsPB = callsetMapBuilder.build();
        int callsetsCount = existingCallsetsPB.getCallsetsCount();
        for (int i=0; i<callsetsCount; i++) {
            String sampleName = existingCallsetsPB.getCallsets(i).getSampleName();
            URI value = inputSampleNameToPath.get(sampleName);
            if (value != null) {
                throw new GenomicsDBException("Duplicate sample name found: "+sampleName+". Sample "+
                        "was originally in "+value);
            }
        }
    }

    /**
     * Returns stream name given a sample name
     * @param sampleName sample name
     * @return stream name
     */
    default String getStreamNameFromSampleName(final String sampleName) {
        return sampleName + "_stream";
    }
}
