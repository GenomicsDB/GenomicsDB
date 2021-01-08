.. java:import:: org.apache.spark.sql.catalyst InternalRow

.. java:import:: org.apache.spark.sql.sources.v2.reader InputPartition

.. java:import:: org.apache.spark.sql.sources.v2.reader InputPartitionReader

.. java:import:: org.apache.spark.sql.types StructType

.. java:import:: java.util Map

.. java:import:: java.util ArrayList

GenomicsDBInputPartition
========================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBInputPartition implements InputPartition<InternalRow>, GenomicsDBInputInterface

Constructors
------------
GenomicsDBInputPartition
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputPartition()
   :outertype: GenomicsDBInputPartition

GenomicsDBInputPartition
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputPartition(String host, GenomicsDBConfiguration gConf, StructType schema)
   :outertype: GenomicsDBInputPartition

GenomicsDBInputPartition
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputPartition(GenomicsDBPartitionInfo partition, ArrayList<GenomicsDBQueryInfo> queryList, GenomicsDBConfiguration gConf, StructType schema)
   :outertype: GenomicsDBInputPartition

Methods
-------
createPartitionReader
^^^^^^^^^^^^^^^^^^^^^

.. java:method:: @Override public InputPartitionReader<InternalRow> createPartitionReader()
   :outertype: GenomicsDBInputPartition

getGenomicsDBVidSchema
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public Map<String, GenomicsDBVidSchema> getGenomicsDBVidSchema()
   :outertype: GenomicsDBInputPartition

getLoader
^^^^^^^^^

.. java:method:: public String getLoader()
   :outertype: GenomicsDBInputPartition

getLoaderIsPB
^^^^^^^^^^^^^

.. java:method:: public boolean getLoaderIsPB()
   :outertype: GenomicsDBInputPartition

getPartitionInfo
^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBPartitionInfo getPartitionInfo()
   :outertype: GenomicsDBInputPartition

getQuery
^^^^^^^^

.. java:method:: public String getQuery()
   :outertype: GenomicsDBInputPartition

getQueryInfoList
^^^^^^^^^^^^^^^^

.. java:method:: public ArrayList<GenomicsDBQueryInfo> getQueryInfoList()
   :outertype: GenomicsDBInputPartition

getQueryIsPB
^^^^^^^^^^^^

.. java:method:: public boolean getQueryIsPB()
   :outertype: GenomicsDBInputPartition

getSchema
^^^^^^^^^

.. java:method:: public StructType getSchema()
   :outertype: GenomicsDBInputPartition

preferredLocations
^^^^^^^^^^^^^^^^^^

.. java:method:: @Override public String[] preferredLocations()
   :outertype: GenomicsDBInputPartition

setGenomicsDBConf
^^^^^^^^^^^^^^^^^

.. java:method:: public void setGenomicsDBConf(GenomicsDBConfiguration g)
   :outertype: GenomicsDBInputPartition

setGenomicsDBSchema
^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setGenomicsDBSchema(StructType s)
   :outertype: GenomicsDBInputPartition

setGenomicsDBVidSchema
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setGenomicsDBVidSchema(Map<String, GenomicsDBVidSchema> v)
   :outertype: GenomicsDBInputPartition

setPartitionInfo
^^^^^^^^^^^^^^^^

.. java:method:: public void setPartitionInfo(GenomicsDBPartitionInfo p)
   :outertype: GenomicsDBInputPartition

setQueryInfoList
^^^^^^^^^^^^^^^^

.. java:method:: public void setQueryInfoList(ArrayList<GenomicsDBQueryInfo> q)
   :outertype: GenomicsDBInputPartition

