.. java:import:: org.genomicsdb GenomicsDBLibLoader

.. java:import:: org.genomicsdb.exception GenomicsDBException

GenomicsDBQuery
===============

.. java:package:: org.genomicsdb.reader
   :noindex:

.. java:type:: public class GenomicsDBQuery

Fields
------
defaultSegmentSize
^^^^^^^^^^^^^^^^^^

.. java:field:: public static final long defaultSegmentSize
   :outertype: GenomicsDBQuery

Methods
-------
connect
^^^^^^^

.. java:method:: public long connect(String workspace, String vidMappingFile, String callsetMappingFile, String referenceGenome, List<String> attributes) throws GenomicsDBException
   :outertype: GenomicsDBQuery

connect
^^^^^^^

.. java:method:: public long connect(String workspace, String vidMappingFile, String callsetMappingFile, String referenceGenome, List<String> attributes, long segmentSize) throws GenomicsDBException
   :outertype: GenomicsDBQuery

connectJSON
^^^^^^^^^^^

.. java:method:: public long connectJSON(String queryJSONFile)
   :outertype: GenomicsDBQuery

connectJSON
^^^^^^^^^^^

.. java:method:: public long connectJSON(String queryJSONFile, String loaderJSONFile)
   :outertype: GenomicsDBQuery

disconnect
^^^^^^^^^^

.. java:method:: public void disconnect(long handle)
   :outertype: GenomicsDBQuery

generateVCF
^^^^^^^^^^^

.. java:method:: public void generateVCF(long handle, String arrayName, List<Pair> columnRanges, List<Pair> rowRanges, String outputFilename, String outputFormat)
   :outertype: GenomicsDBQuery

generateVCF
^^^^^^^^^^^

.. java:method:: public void generateVCF(long handle, String arrayName, List<Pair> columnRanges, List<Pair> rowRanges, String outputFilename, String outputFormat, boolean overwrite)
   :outertype: GenomicsDBQuery

generateVCF
^^^^^^^^^^^

.. java:method:: public void generateVCF(long handle, String outputFilename, String outputFormat)
   :outertype: GenomicsDBQuery

generateVCF
^^^^^^^^^^^

.. java:method:: public void generateVCF(long handle, String outputFilename, String outputFormat, boolean overwrite)
   :outertype: GenomicsDBQuery

queryVariantCalls
^^^^^^^^^^^^^^^^^

.. java:method:: public List<Interval> queryVariantCalls(long handle, String arrayName)
   :outertype: GenomicsDBQuery

queryVariantCalls
^^^^^^^^^^^^^^^^^

.. java:method:: public List<Interval> queryVariantCalls(long handle, String arrayName, List<Pair> columnRanges)
   :outertype: GenomicsDBQuery

queryVariantCalls
^^^^^^^^^^^^^^^^^

.. java:method:: public List<Interval> queryVariantCalls(long handle, String arrayName, List<Pair> columnRanges, List<Pair> rowRanges)
   :outertype: GenomicsDBQuery

version
^^^^^^^

.. java:method:: public String version()
   :outertype: GenomicsDBQuery

