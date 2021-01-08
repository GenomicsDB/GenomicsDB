.. java:import:: htsjdk.tribble CloseableTribbleIterator

.. java:import:: htsjdk.tribble FeatureCodec

.. java:import:: htsjdk.tribble.readers PositionalBufferedStream

.. java:import:: htsjdk.variant.bcf2 BCF2Codec

.. java:import:: htsjdk.variant.variantcontext Allele

.. java:import:: htsjdk.variant.variantcontext Genotype

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: org.apache.spark.sql.catalyst InternalRow

.. java:import:: org.apache.spark.sql.catalyst.util ArrayData

.. java:import:: org.apache.spark.sql.sources.v2.reader InputPartitionReader

.. java:import:: org.apache.spark.sql.types StructField

.. java:import:: org.apache.spark.unsafe.types UTF8String

.. java:import:: org.genomicsdb.reader GenomicsDBFeatureReader

.. java:import:: org.genomicsdb.model Coordinates

.. java:import:: org.genomicsdb.model GenomicsDBExportConfiguration

.. java:import:: org.json.simple.parser ParseException

.. java:import:: scala.collection JavaConverters

.. java:import:: java.io IOException

.. java:import:: java.lang.reflect Array

.. java:import:: java.util ArrayList

.. java:import:: java.util List

.. java:import:: java.util Map

.. java:import:: java.util Optional

GenomicsDBInputPartitionReader
==============================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBInputPartitionReader implements InputPartitionReader<InternalRow>

Constructors
------------
GenomicsDBInputPartitionReader
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputPartitionReader()
   :outertype: GenomicsDBInputPartitionReader

GenomicsDBInputPartitionReader
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputPartitionReader(GenomicsDBInputPartition iPartition)
   :outertype: GenomicsDBInputPartitionReader

Methods
-------
close
^^^^^

.. java:method:: public void close()
   :outertype: GenomicsDBInputPartitionReader

get
^^^

.. java:method:: public InternalRow get()
   :outertype: GenomicsDBInputPartitionReader

next
^^^^

.. java:method:: public boolean next()
   :outertype: GenomicsDBInputPartitionReader

