.. java:import:: org.apache.hadoop.io Writable

.. java:import:: org.apache.hadoop.io Text

.. java:import:: org.apache.hadoop.mapreduce InputSplit

.. java:import:: org.apache.log4j Logger

.. java:import:: org.apache.spark Partition

.. java:import:: java.util ArrayList

.. java:import:: java.io DataInput

.. java:import:: java.io DataOutput

.. java:import:: java.io IOException

GenomicsDBInputSplit
====================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBInputSplit extends InputSplit implements Writable, GenomicsDBInputInterface

Fields
------
hosts
^^^^^

.. java:field::  String[] hosts
   :outertype: GenomicsDBInputSplit

length
^^^^^^

.. java:field::  long length
   :outertype: GenomicsDBInputSplit

logger
^^^^^^

.. java:field::  Logger logger
   :outertype: GenomicsDBInputSplit

partition
^^^^^^^^^

.. java:field::  GenomicsDBPartitionInfo partition
   :outertype: GenomicsDBInputSplit

queryRangeList
^^^^^^^^^^^^^^

.. java:field::  ArrayList<GenomicsDBQueryInfo> queryRangeList
   :outertype: GenomicsDBInputSplit

Constructors
------------
GenomicsDBInputSplit
^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputSplit()
   :outertype: GenomicsDBInputSplit

   Default constructor

GenomicsDBInputSplit
^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputSplit(String loc)
   :outertype: GenomicsDBInputSplit

GenomicsDBInputSplit
^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputSplit(GenomicsDBPartitionInfo partition, ArrayList<GenomicsDBQueryInfo> queryRangeList)
   :outertype: GenomicsDBInputSplit

GenomicsDBInputSplit
^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputSplit(long length)
   :outertype: GenomicsDBInputSplit

Methods
-------
getLength
^^^^^^^^^

.. java:method:: public long getLength() throws IOException, InterruptedException
   :outertype: GenomicsDBInputSplit

getLocations
^^^^^^^^^^^^

.. java:method:: public String[] getLocations() throws IOException, InterruptedException
   :outertype: GenomicsDBInputSplit

   Called from \ :java:ref:`org.apache.spark.rdd.NewHadoopRDD.getPreferredLocations(Partition)`\

   :return: locations The values of locations or nodes are passed from host file in GenomicsDBConfiguration or hadoopConfiguration in SparkContext

getPartitionInfo
^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBPartitionInfo getPartitionInfo()
   :outertype: GenomicsDBInputSplit

getQueryInfoList
^^^^^^^^^^^^^^^^

.. java:method:: public ArrayList<GenomicsDBQueryInfo> getQueryInfoList()
   :outertype: GenomicsDBInputSplit

readFields
^^^^^^^^^^

.. java:method:: public void readFields(DataInput dataInput) throws IOException
   :outertype: GenomicsDBInputSplit

setHost
^^^^^^^

.. java:method:: public void setHost(String loc)
   :outertype: GenomicsDBInputSplit

setPartitionInfo
^^^^^^^^^^^^^^^^

.. java:method:: public void setPartitionInfo(GenomicsDBPartitionInfo p)
   :outertype: GenomicsDBInputSplit

setQueryInfoList
^^^^^^^^^^^^^^^^

.. java:method:: public void setQueryInfoList(ArrayList<GenomicsDBQueryInfo> q)
   :outertype: GenomicsDBInputSplit

write
^^^^^

.. java:method:: public void write(DataOutput dataOutput) throws IOException
   :outertype: GenomicsDBInputSplit

