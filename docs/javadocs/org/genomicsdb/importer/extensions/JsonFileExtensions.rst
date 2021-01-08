.. java:import:: com.googlecode.protobuf.format JsonFormat

.. java:import:: htsjdk.variant.variantcontext.writer VariantContextWriter

.. java:import:: htsjdk.variant.variantcontext.writer VariantContextWriterBuilder

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: htsjdk.variant.vcf VCFHeaderLine

.. java:import:: org.genomicsdb.model GenomicsDBImportConfiguration

.. java:import:: org.genomicsdb GenomicsDBUtils

.. java:import:: org.genomicsdb.model GenomicsDBCallsetsMapProto

.. java:import:: org.genomicsdb.model GenomicsDBVidMapProto

.. java:import:: java.util Set

JsonFileExtensions
==================

.. java:package:: org.genomicsdb.importer.extensions
   :noindex:

.. java:type:: public interface JsonFileExtensions

Methods
-------
dumpTemporaryLoaderJSONFile
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  File dumpTemporaryLoaderJSONFile(GenomicsDBImportConfiguration.ImportConfiguration importConfiguration, String filename) throws IOException
   :outertype: JsonFileExtensions

   Create a JSON file from the import configuration

   :param importConfiguration: The configuration object
   :param filename: File to dump the loader JSON to
   :throws FileNotFoundException: PrintWriter throws this exception when file not found
   :return: New file (with the specified name) written to local storage

writeCallsetMapJSONFile
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  void writeCallsetMapJSONFile(String outputCallsetMapJSONFilePath, GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB)
   :outertype: JsonFileExtensions

   Writes a JSON file from a vidmap protobuf object. This method is not called implicitly by GenomicsDB constructor. It needs to be called explicitly if the user wants these objects to be written. Called explicitly from GATK-4 GenomicsDBImport tool

   :param outputCallsetMapJSONFilePath: Full path of file to be written
   :param callsetMappingPB: Protobuf callset map object

writeVcfHeaderFile
^^^^^^^^^^^^^^^^^^

.. java:method::  void writeVcfHeaderFile(String outputVcfHeaderFilePath, Set<VCFHeaderLine> headerLines)
   :outertype: JsonFileExtensions

   Writes a VCF Header file from a set of VCF header lines. This method is not called implicitly by GenomicsDB constructor. It needs to be called explicitly if the user wants these objects to be written. Called explicitly from GATK-4 GenomicsDBImport tool

   :param outputVcfHeaderFilePath: Full path of file to be written
   :param headerLines: Set of header lines

writeVidMapJSONFile
^^^^^^^^^^^^^^^^^^^

.. java:method::  void writeVidMapJSONFile(String outputVidMapJSONFilePath, GenomicsDBVidMapProto.VidMappingPB vidMappingPB)
   :outertype: JsonFileExtensions

   Writes a JSON file from a set of VCF header lines. This method is not called implicitly by GenomicsDB constructor. It needs to be called explicitly if the user wants these objects to be written. Called explicitly from GATK-4 GenomicsDBImport tool

   :param outputVidMapJSONFilePath: Full path of file to be written
   :param vidMappingPB: ProtoBuf data structure for vid mapping

