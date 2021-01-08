.. java:import:: com.fasterxml.jackson.databind ObjectMapper

.. java:import:: java.io IOException

.. java:import:: java.io Serializable

GenomicsDBPartitionInfo
=======================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type::  class GenomicsDBPartitionInfo implements Serializable

   Maintain global information on how variant data is partitioned across multiple nodes (if loaded with mpiexec) The assumption is that every object described one partition as: {"begin": 0, "workspace":"/home/forall/tiledb-ws", "array": "test0", "vcf_output_filename":"/home/forall/lihao/mytest/test0.vcf.gz" }

Constructors
------------
GenomicsDBPartitionInfo
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBPartitionInfo(long pos, String workspace, String arrayName, String vcfOutput)
   :outertype: GenomicsDBPartitionInfo

GenomicsDBPartitionInfo
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBPartitionInfo(GenomicsDBPartitionInfo copy)
   :outertype: GenomicsDBPartitionInfo

Methods
-------
equals
^^^^^^

.. java:method:: public boolean equals(Object obj)
   :outertype: GenomicsDBPartitionInfo

getArrayName
^^^^^^^^^^^^

.. java:method:: public String getArrayName()
   :outertype: GenomicsDBPartitionInfo

getBeginPosition
^^^^^^^^^^^^^^^^

.. java:method:: public long getBeginPosition()
   :outertype: GenomicsDBPartitionInfo

getVcfOutputFileName
^^^^^^^^^^^^^^^^^^^^

.. java:method:: public String getVcfOutputFileName()
   :outertype: GenomicsDBPartitionInfo

getWorkspace
^^^^^^^^^^^^

.. java:method:: public String getWorkspace()
   :outertype: GenomicsDBPartitionInfo

hashCode
^^^^^^^^

.. java:method:: @Override public int hashCode()
   :outertype: GenomicsDBPartitionInfo

toJSON
^^^^^^

.. java:method:: public String toJSON() throws IOException
   :outertype: GenomicsDBPartitionInfo

   Create a JSON string from the partition info object

   :throws IOException: write throws IOException
   :return: A JSON string of the partition info

toString
^^^^^^^^

.. java:method:: @Override public String toString()
   :outertype: GenomicsDBPartitionInfo

