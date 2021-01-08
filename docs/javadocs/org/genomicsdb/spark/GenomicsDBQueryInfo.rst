.. java:import:: java.io Serializable

GenomicsDBQueryInfo
===================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type::  class GenomicsDBQueryInfo implements Serializable

   Maintain global information on query ranges which is used to create GenomicsDBInputSplits

Constructors
------------
GenomicsDBQueryInfo
^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryInfo(long start, long end)
   :outertype: GenomicsDBQueryInfo

GenomicsDBQueryInfo
^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryInfo(GenomicsDBQueryInfo copy)
   :outertype: GenomicsDBQueryInfo

Methods
-------
equals
^^^^^^

.. java:method:: public boolean equals(Object obj)
   :outertype: GenomicsDBQueryInfo

getBeginPosition
^^^^^^^^^^^^^^^^

.. java:method:: public long getBeginPosition()
   :outertype: GenomicsDBQueryInfo

getEndPosition
^^^^^^^^^^^^^^

.. java:method:: public long getEndPosition()
   :outertype: GenomicsDBQueryInfo

hashCode
^^^^^^^^

.. java:method:: public int hashCode()
   :outertype: GenomicsDBQueryInfo

