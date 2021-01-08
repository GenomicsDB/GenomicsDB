.. java:import:: htsjdk.tribble FeatureReader

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.vcf VCFHeaderLine

.. java:import:: java.util.function Function

.. java:import:: java.lang Integer

.. java:import:: java.util.stream IntStream

.. java:import:: java.util.stream Collectors

.. java:import:: java.net URI

ImportConfig
============

.. java:package:: org.genomicsdb.model
   :noindex:

.. java:type:: public class ImportConfig

   This implementation extends what is in GenomicsDBImportConfiguration. Add extra data that is needed for parallel import.

Constructors
------------
ImportConfig
^^^^^^^^^^^^

.. java:constructor:: public ImportConfig(GenomicsDBImportConfiguration.ImportConfiguration importConfiguration, boolean validateSampleToReaderMap, boolean passAsVcf, int batchSize, Set<VCFHeaderLine> mergedHeader, Map<String, URI> sampleNameToVcfPath, Func<Map<String, URI>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator, boolean incrementalImport)
   :outertype: ImportConfig

   Main ImportConfig constructor

   :param importConfiguration: GenomicsDBImportConfiguration protobuf object
   :param validateSampleToReaderMap: Flag for validating sample to reader map
   :param passAsVcf: Flag for indicating that a VCF is being passed
   :param batchSize: Batch size
   :param mergedHeader: Required header
   :param sampleNameToVcfPath: Sample name to VCF path map
   :param sampleToReaderMapCreator: Function used for creating sampleToReaderMap
   :param incrementalImport: Flag for indicating incremental import

ImportConfig
^^^^^^^^^^^^

.. java:constructor:: protected ImportConfig()
   :outertype: ImportConfig

Methods
-------
getBatchSize
^^^^^^^^^^^^

.. java:method:: public int getBatchSize()
   :outertype: ImportConfig

getFunctionToCallOnBatchCompletion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public Function<BatchCompletionCallbackFunctionArgument, Void> getFunctionToCallOnBatchCompletion()
   :outertype: ImportConfig

getImportConfiguration
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBImportConfiguration.ImportConfiguration getImportConfiguration()
   :outertype: ImportConfig

getMergedHeader
^^^^^^^^^^^^^^^

.. java:method:: public Set<VCFHeaderLine> getMergedHeader()
   :outertype: ImportConfig

getOutputCallsetmapJsonFile
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public String getOutputCallsetmapJsonFile()
   :outertype: ImportConfig

getOutputVcfHeaderFile
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public String getOutputVcfHeaderFile()
   :outertype: ImportConfig

getOutputVidmapJsonFile
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public String getOutputVidmapJsonFile()
   :outertype: ImportConfig

getSampleNameToVcfPath
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public Map<String, URI> getSampleNameToVcfPath()
   :outertype: ImportConfig

isIncrementalImport
^^^^^^^^^^^^^^^^^^^

.. java:method:: public boolean isIncrementalImport()
   :outertype: ImportConfig

isPassAsVcf
^^^^^^^^^^^

.. java:method:: public boolean isPassAsVcf()
   :outertype: ImportConfig

isUseSamplesInOrder
^^^^^^^^^^^^^^^^^^^

.. java:method:: public boolean isUseSamplesInOrder()
   :outertype: ImportConfig

isValidateSampleToReaderMap
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public boolean isValidateSampleToReaderMap()
   :outertype: ImportConfig

sampleToReaderMapCreator
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public Func<Map<String, URI>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator()
   :outertype: ImportConfig

setBatchSize
^^^^^^^^^^^^

.. java:method:: public void setBatchSize(int batchSize)
   :outertype: ImportConfig

setFunctionToCallOnBatchCompletion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setFunctionToCallOnBatchCompletion(Function<BatchCompletionCallbackFunctionArgument, Void> functionToCallOnBatchCompletion)
   :outertype: ImportConfig

setImportConfiguration
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setImportConfiguration(GenomicsDBImportConfiguration.ImportConfiguration importConfiguration)
   :outertype: ImportConfig

setIncrementalImport
^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setIncrementalImport(boolean incrementalImport)
   :outertype: ImportConfig

setMergedHeader
^^^^^^^^^^^^^^^

.. java:method:: public void setMergedHeader(Set<VCFHeaderLine> mergedHeader)
   :outertype: ImportConfig

setOutputCallsetmapJsonFile
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setOutputCallsetmapJsonFile(String outputCallsetMapJsonFile)
   :outertype: ImportConfig

setOutputVcfHeaderFile
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setOutputVcfHeaderFile(String outputVcfHeaderFile)
   :outertype: ImportConfig

setOutputVidmapJsonFile
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setOutputVidmapJsonFile(String outputVidMapJsonFile)
   :outertype: ImportConfig

setPassAsVcf
^^^^^^^^^^^^

.. java:method:: public void setPassAsVcf(boolean passAsVcf)
   :outertype: ImportConfig

setSampleNameToVcfPath
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setSampleNameToVcfPath(Map<String, URI> sampleNameToVcfPath)
   :outertype: ImportConfig

setSampleToReaderMapCreator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setSampleToReaderMapCreator(Func<Map<String, URI>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator)
   :outertype: ImportConfig

setUseSamplesInOrder
^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setUseSamplesInOrder(boolean useSamplesInOrder)
   :outertype: ImportConfig

setValidateSampleToReaderMap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setValidateSampleToReaderMap(boolean validateSampleToReaderMap)
   :outertype: ImportConfig

validateChromosomeIntervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  void validateChromosomeIntervals()
   :outertype: ImportConfig

