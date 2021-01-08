.. java:import:: htsjdk.tribble CloseableTribbleIterator

.. java:import:: htsjdk.tribble Feature

.. java:import:: htsjdk.tribble FeatureCodec

.. java:import:: htsjdk.tribble FeatureCodecHeader

.. java:import:: htsjdk.variant.bcf2 BCF2Codec

.. java:import:: org.genomicsdb.model GenomicsDBExportConfiguration

.. java:import:: java.io IOException

.. java:import:: java.util Iterator

.. java:import:: java.util List

.. java:import:: java.util Optional

.. java:import:: java.util OptionalInt

.. java:import:: java.util.stream Collectors

.. java:import:: java.util Collections

.. java:import:: java.lang RuntimeException

GenomicsDBFeatureIterator
=========================

.. java:package:: org.genomicsdb.reader
   :noindex:

.. java:type:: public class GenomicsDBFeatureIterator<T extends Feature, SOURCE> implements CloseableTribbleIterator<T>

   Iterator over \ :java:ref:`htsjdk.variant.variantcontext.VariantContext`\  objects. Uses \ :java:ref:`GenomicsDBQueryStream`\  to obtain combined gVCF records (as BCF2) from TileDB/GenomicsDB

Constructors
------------
GenomicsDBFeatureIterator
^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor::  GenomicsDBFeatureIterator(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, Optional<List<String>> arrayNames, FeatureCodecHeader featureCodecHeader, FeatureCodec<T, SOURCE> codec) throws IOException
   :outertype: GenomicsDBFeatureIterator

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf object
   :param arrayNames: List of array names
   :param featureCodecHeader: htsjdk Feature codec header
   :param codec: FeatureCodec, currently only \ :java:ref:`htsjdk.variant.bcf2.BCF2Codec`\  and \ :java:ref:`htsjdk.variant.vcf.VCFCodec`\  are tested
   :throws IOException: when data cannot be read from the stream

GenomicsDBFeatureIterator
^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor::  GenomicsDBFeatureIterator(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, Optional<List<String>> arrayNames, FeatureCodecHeader featureCodecHeader, FeatureCodec<T, SOURCE> codec, String chr, OptionalInt start, OptionalInt end) throws IOException
   :outertype: GenomicsDBFeatureIterator

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param arrayNames: List of array names
   :param featureCodecHeader: htsjdk Feature codec header
   :param codec: FeatureCodec, currently only \ :java:ref:`htsjdk.variant.bcf2.BCF2Codec`\  and \ :java:ref:`htsjdk.variant.vcf.VCFCodec`\  are tested
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :throws IOException: when data cannot be read from the stream

Methods
-------
close
^^^^^

.. java:method:: @Override public void close()
   :outertype: GenomicsDBFeatureIterator

hasNext
^^^^^^^

.. java:method:: @Override public boolean hasNext()
   :outertype: GenomicsDBFeatureIterator

iterator
^^^^^^^^

.. java:method:: @Override public Iterator<T> iterator()
   :outertype: GenomicsDBFeatureIterator

next
^^^^

.. java:method:: @Override public T next()
   :outertype: GenomicsDBFeatureIterator

remove
^^^^^^

.. java:method:: @Override public void remove()
   :outertype: GenomicsDBFeatureIterator

