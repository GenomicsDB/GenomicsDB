.. java:import:: org.genomicsdb.importer Constants

.. java:import:: org.genomicsdb.model GenomicsDBVidMapProto

.. java:import:: org.genomicsdb GenomicsDBUtils

.. java:import:: java.util ArrayList

.. java:import:: java.util List

.. java:import:: java.util Set

.. java:import:: java.util Map

.. java:import:: java.util HashMap

VidMapExtensions
================

.. java:package:: org.genomicsdb.importer.extensions
   :noindex:

.. java:type:: public interface VidMapExtensions

Fields
------
vcfHeaderLineFilterIdx
^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static int vcfHeaderLineFilterIdx
   :outertype: VidMapExtensions

vcfHeaderLineFormatIdx
^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static int vcfHeaderLineFormatIdx
   :outertype: VidMapExtensions

vcfHeaderLineInfoIdx
^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static int vcfHeaderLineInfoIdx
   :outertype: VidMapExtensions

vcfHeaderLineTypeString
^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static String[] vcfHeaderLineTypeString
   :outertype: VidMapExtensions

Methods
-------
checkForDuplicateVcfFieldNames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static DuplicateVcfFieldNamesCheckResult checkForDuplicateVcfFieldNames(String vcfFieldName, int vcfHeaderLineTypeIdx, String vcfType, String vcfLengthDescriptor, List<GenomicsDBVidMapProto.GenomicsDBFieldInfo> infoFields, Map<String, List<DuplicateFieldInfoTrackClass>> fieldNameToDuplicateCheckInfo)
   :outertype: VidMapExtensions

generateVidMapFromFile
^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBVidMapProto.VidMappingPB generateVidMapFromFile(String vidFile) throws ParseException
   :outertype: VidMapExtensions

   Generate the ProtoBuf data structure for vid mapping from the existing vid file

   :param vidFile: file with existing vid info
   :throws ParseException: when there is an error parsing existing vid json
   :return: a vid map containing all field names, lengths and types from the merged GVCF header. for incremental import case

generateVidMapFromMergedHeader
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  GenomicsDBVidMapProto.VidMappingPB generateVidMapFromMergedHeader(Set<VCFHeaderLine> mergedHeader)
   :outertype: VidMapExtensions

   Generate the ProtoBuf data structure for vid mapping Also remember, contigs are 1-based which means TileDB column offsets should start from 1

   :param mergedHeader: Header from all input GVCFs are merged to create one vid map for all
   :return: a vid map containing all field names, lengths and types from the merged GVCF header

getVcfHeaderLength
^^^^^^^^^^^^^^^^^^

.. java:method::  String getVcfHeaderLength(VCFHeaderLine headerLine)
   :outertype: VidMapExtensions

   Maps the "Number" from INFO or FORMAT fields in VCF header to GenomicsDB lengths. If unbounded, the field becomes a variable length attribute (discussed in more detail in TileDB tutorial tiledb.org). If "A", "R" or "G", the field also becomes a variable length attribute, but the order and length depends on order/number of alleles or genotypes. Integers remain integers

   :param headerLine: Info or Format header line from VCF
   :return: VAR, A, R, G, or integer values from VCF

