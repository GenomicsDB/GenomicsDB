.. java:import:: org.genomicsdb.reader GenomicsDBFeatureReader

.. java:import:: htsjdk.tribble CloseableTribbleIterator

.. java:import:: htsjdk.tribble Feature

.. java:import:: org.apache.hadoop.classification InterfaceAudience

.. java:import:: org.apache.hadoop.classification InterfaceStability

.. java:import:: org.apache.hadoop.mapreduce InputSplit

.. java:import:: org.apache.hadoop.mapreduce RecordReader

.. java:import:: org.apache.hadoop.mapreduce TaskAttemptContext

.. java:import:: org.apache.spark.sql.sources.v2.reader InputPartition

.. java:import:: java.io IOException

GenomicsDBRecordReader
======================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: @InterfaceAudience.Public @InterfaceStability.Stable public class GenomicsDBRecordReader<VCONTEXT extends Feature, SOURCE> extends RecordReader<String, VCONTEXT>

Constructors
------------
GenomicsDBRecordReader
^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor::  GenomicsDBRecordReader(GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader)
   :outertype: GenomicsDBRecordReader

Methods
-------
close
^^^^^

.. java:method:: public void close() throws IOException
   :outertype: GenomicsDBRecordReader

getCurrentKey
^^^^^^^^^^^^^

.. java:method:: @Override public String getCurrentKey() throws IOException, InterruptedException
   :outertype: GenomicsDBRecordReader

getCurrentValue
^^^^^^^^^^^^^^^

.. java:method:: public VCONTEXT getCurrentValue() throws IOException, InterruptedException
   :outertype: GenomicsDBRecordReader

getProgress
^^^^^^^^^^^

.. java:method:: public float getProgress() throws IOException, InterruptedException
   :outertype: GenomicsDBRecordReader

hasNext
^^^^^^^

.. java:method:: public Boolean hasNext()
   :outertype: GenomicsDBRecordReader

initialize
^^^^^^^^^^

.. java:method:: public void initialize(InputSplit inputSplit, TaskAttemptContext taskAttemptContext) throws IOException, InterruptedException
   :outertype: GenomicsDBRecordReader

nextKeyValue
^^^^^^^^^^^^

.. java:method:: public boolean nextKeyValue() throws IOException, InterruptedException
   :outertype: GenomicsDBRecordReader

