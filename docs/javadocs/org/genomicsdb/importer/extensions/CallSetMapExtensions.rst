.. java:import:: htsjdk.tribble FeatureReader

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: org.genomicsdb.model GenomicsDBCallsetsMapProto

.. java:import:: org.genomicsdb GenomicsDBUtils

.. java:import:: org.genomicsdb.exception GenomicsDBException

.. java:import:: java.net URI

CallSetMapExtensions
====================

.. java:package:: org.genomicsdb.importer.extensions
   :noindex:

.. java:type:: public interface CallSetMapExtensions

Methods
-------
checkDuplicateCallsetsForIncrementalImport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  void checkDuplicateCallsetsForIncrementalImport(String existingCallsetsJSON, Map<String, URI> inputSampleNameToPath) throws ParseException, GenomicsDBException
   :outertype: CallSetMapExtensions

   Checks for duplicates between existing and incremental callsets

   :param existingCallsetsJSON: json string with existing callset
   :param inputSampleNameToPath: map from sample name to VCF/BCF file path
   :throws ParseException: when there is an error parsing callset jsons
   :throws GenomicsDBException: when duplicate samples are found between existing and new callsets

generateSortedCallSetMap
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(List<String> sampleNames, boolean useSamplesInOrderProvided)
   :outertype: CallSetMapExtensions

   Creates a sorted list of callsets and generates unique TileDB row indices for them. Sorted to maintain order between distributed share-nothing load processes.

   :param sampleNames: list of sample names
   :param useSamplesInOrderProvided: If True, do not sort the samples, use the order they appear in
   :return: Mappings of callset (sample) names to TileDB rows

generateSortedCallSetMap
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(List<String> sampleNames, boolean useSamplesInOrderProvided, long lbRowIdx)
   :outertype: CallSetMapExtensions

   Creates a sorted list of callsets and generates unique TileDB row indices for them. Sorted to maintain order between distributed share-nothing load processes. This method is synchronized to block multiple invocations (if by any chance) disturb the order in which TileDB row indexes are generated

   :param sampleNames: list of sample names
   :param useSamplesInOrderProvided: If True, do not sort the samples, use the order they appear in
   :param lbRowIdx: Smallest row idx which should be imported by this object
   :return: Mappings of callset (sample) names to TileDB rows

generateSortedCallSetMap
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(Map<String, String> inputSampleNameToStreamName, boolean useSamplesInOrderProvided, long lbRowIdx)
   :outertype: CallSetMapExtensions

   Creates a sorted list of callsets and generates unique TileDB row indices for them. Sorted to maintain order between distributed share-nothing load processes. This method is synchronized to block multiple invocations (if by any chance) disturb the order in which TileDB row indexes are generated

   :param inputSampleNameToStreamName: map from sample name to VCF/BCF file path
   :param useSamplesInOrderProvided: If True, do not sort the samples, use the order they appear in
   :param lbRowIdx: Smallest row idx which should be imported by this object
   :return: Mappings of callset (sample) names to TileDB rows

generateSortedCallSetMap
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(LinkedHashMap<String, String> sampleNameToStreamName, long lbRowIdx)
   :outertype: CallSetMapExtensions

   Creates CallSets Protobuf structure for given map

   :param sampleNameToStreamName: map from sample name to VCF/BCF file path use the order they appear in
   :param lbRowIdx: Smallest row idx which should be imported by this object
   :return: Mappings of callset (sample) names to TileDB rows

generateSortedCallSetMap
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(Map<String, FeatureReader<VariantContext>> sampleToReaderMap, boolean validateSampleToReaderMap, boolean useSamplesInOrderProvided)
   :outertype: CallSetMapExtensions

   Creates a sorted list of callsets and generates unique TileDB row indices for them. Sorted to maintain order between distributed share-nothing load processes.

   Assume one sample per input GVCF file

   :param sampleToReaderMap: Variant Readers objects of the input GVCF files
   :param useSamplesInOrderProvided: If True, do not sort the samples, use the order they appear in
   :param validateSampleToReaderMap: If True, check i) whether sample names are consistent with headers and ii) feature readers are valid in sampleToReaderMap
   :return: Mappings of callset (sample) names to TileDB rows

generateSortedCallSetMap
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(Map<String, FeatureReader<VariantContext>> sampleToReaderMap, boolean useSamplesInOrderProvided, boolean validateSampleMap, long lbRowIdx)
   :outertype: CallSetMapExtensions

   Creates a sorted list of callsets and generates unique TileDB row indices for them. Sorted to maintain order between distributed share-nothing load processes.

   Assume one sample per input GVCF file

   :param sampleToReaderMap: Variant Readers objects of the input GVCF files
   :param useSamplesInOrderProvided: If True, do not sort the samples, use the order they appear in
   :param validateSampleMap: Check i) whether sample names are consistent with headers and ii) feature readers are valid in sampleToReaderMap
   :param lbRowIdx: Smallest row idx which should be imported by this object
   :return: Mappings of callset (sample) names to TileDB rows

generateSortedCallSetMapFromNameToPathMap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMapFromNameToPathMap(Map<String, URI> inputSampleNameToPath, boolean useSamplesInOrderProvided, long lbRowIdx)
   :outertype: CallSetMapExtensions

   Creates a sorted list of callsets and generates unique TileDB row indices for them. Sorted to maintain order between distributed share-nothing load processes. This method is synchronized to block multiple invocations (if by any chance) disturb the order in which TileDB row indexes are generated

   :param inputSampleNameToPath: map from sample name to VCF/BCF file path
   :param useSamplesInOrderProvided: If True, do not sort the samples, use the order they appear in
   :param lbRowIdx: Smallest row idx which should be imported by this object
   :return: Mappings of callset (sample) names to TileDB rows

getStreamNameFromSampleName
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  String getStreamNameFromSampleName(String sampleName)
   :outertype: CallSetMapExtensions

   Returns stream name given a sample name

   :param sampleName: sample name
   :return: stream name

mergeCallsetsForIncrementalImport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBCallsetsMapProto.CallsetMappingPB mergeCallsetsForIncrementalImport(String callsetMapJSONFilePath, Map<String, URI> inputSampleNameToPath, GenomicsDBCallsetsMapProto.CallsetMappingPB newCallsetMapPB) throws ParseException
   :outertype: CallSetMapExtensions

   Merges incremental import's callsets with existing callsets. Also create a copy of original callset file to aid with recovery if the incremental import goes awry

   :param callsetMapJSONFilePath: path to existing callset map file
   :param inputSampleNameToPath: map from sample name to VCF/BCF file path
   :param newCallsetMapPB: callset mapping protobuf for new callsets
   :throws ParseException: when there is an error parsing callset jsons
   :return: merged callsets for callsets to TileDB rows

