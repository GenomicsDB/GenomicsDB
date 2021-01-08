.. java:import:: org.apache.spark.sql.types StructType

.. java:import:: java.util Map

.. java:import:: java.util ArrayList

GenomicsDBInputInterface
========================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type::  interface GenomicsDBInputInterface

   This interface declares some methods to be used by classes that define quantums of data that are divided from a large amount of data to be distributed to many workers. Currently intended to be implemeented by GenomicsDBInputSplit and GenomicsDBInputPartition

Methods
-------
setGenomicsDBConf
^^^^^^^^^^^^^^^^^

.. java:method::  void setGenomicsDBConf(GenomicsDBConfiguration g)
   :outertype: GenomicsDBInputInterface

setGenomicsDBSchema
^^^^^^^^^^^^^^^^^^^

.. java:method::  void setGenomicsDBSchema(StructType s)
   :outertype: GenomicsDBInputInterface

setGenomicsDBVidSchema
^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  void setGenomicsDBVidSchema(Map<String, GenomicsDBVidSchema> vMap)
   :outertype: GenomicsDBInputInterface

setHost
^^^^^^^

.. java:method::  void setHost(String h)
   :outertype: GenomicsDBInputInterface

setPartitionInfo
^^^^^^^^^^^^^^^^

.. java:method::  void setPartitionInfo(GenomicsDBPartitionInfo p)
   :outertype: GenomicsDBInputInterface

setQueryInfoList
^^^^^^^^^^^^^^^^

.. java:method::  void setQueryInfoList(ArrayList<GenomicsDBQueryInfo> q)
   :outertype: GenomicsDBInputInterface

