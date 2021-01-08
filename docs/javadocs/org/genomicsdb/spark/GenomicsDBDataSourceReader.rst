.. java:import:: org.apache.spark.sql.catalyst InternalRow

.. java:import:: org.apache.spark.sql.sources.v2 DataSourceOptions

.. java:import:: org.apache.spark.sql.sources.v2.reader DataSourceReader

.. java:import:: org.apache.spark.sql.sources.v2.reader InputPartition

.. java:import:: org.json.simple JSONArray

.. java:import:: org.json.simple JSONObject

.. java:import:: org.json.simple.parser JSONParser

.. java:import:: org.json.simple.parser ParseException

.. java:import:: scala.collection JavaConverters

.. java:import:: java.io FileNotFoundException

.. java:import:: java.io FileReader

.. java:import:: java.io IOException

.. java:import:: java.util HashMap

.. java:import:: java.util List

.. java:import:: java.util Map

GenomicsDBDataSourceReader
==========================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBDataSourceReader implements DataSourceReader

Fields
------
input
^^^^^

.. java:field::  GenomicsDBInput<GenomicsDBInputPartition> input
   :outertype: GenomicsDBDataSourceReader

Constructors
------------
GenomicsDBDataSourceReader
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBDataSourceReader()
   :outertype: GenomicsDBDataSourceReader

GenomicsDBDataSourceReader
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBDataSourceReader(StructType schema, DataSourceOptions options)
   :outertype: GenomicsDBDataSourceReader

GenomicsDBDataSourceReader
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBDataSourceReader(DataSourceOptions options)
   :outertype: GenomicsDBDataSourceReader

Methods
-------
planInputPartitions
^^^^^^^^^^^^^^^^^^^

.. java:method:: @Override @SuppressWarnings public List<InputPartition<InternalRow>> planInputPartitions()
   :outertype: GenomicsDBDataSourceReader

readSchema
^^^^^^^^^^

.. java:method:: @Override public StructType readSchema()
   :outertype: GenomicsDBDataSourceReader

