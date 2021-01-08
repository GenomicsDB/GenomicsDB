.. java:import:: htsjdk.tribble FeatureReader

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.vcf VCFHeaderLine

.. java:import:: java.util.function Function

.. java:import:: java.lang Integer

.. java:import:: java.util.stream IntStream

.. java:import:: java.util.stream Collectors

.. java:import:: java.net URI

ImportConfig.SortGenomicsDBPartition
====================================

.. java:package:: org.genomicsdb.model
   :noindex:

.. java:type::  class SortGenomicsDBPartition implements Comparator<Integer>
   :outertype: ImportConfig

Constructors
------------
SortGenomicsDBPartition
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public SortGenomicsDBPartition(List<GenomicsDBImportConfiguration.Partition> partitionList)
   :outertype: ImportConfig.SortGenomicsDBPartition

Methods
-------
compare
^^^^^^^

.. java:method:: public int compare(Integer lIdx, Integer rIdx)
   :outertype: ImportConfig.SortGenomicsDBPartition

