.. java:import:: org.genomicsdb GenomicsDBLibLoader

.. java:import:: org.genomicsdb.exception GenomicsDBException

GenomicsDBQuery.VariantCall
===========================

.. java:package:: org.genomicsdb.reader
   :noindex:

.. java:type:: public static class VariantCall
   :outertype: GenomicsDBQuery

Fields
------
colIndex
^^^^^^^^

.. java:field::  long colIndex
   :outertype: GenomicsDBQuery.VariantCall

contigName
^^^^^^^^^^

.. java:field::  String contigName
   :outertype: GenomicsDBQuery.VariantCall

genomicFields
^^^^^^^^^^^^^

.. java:field::  Map<String, Object> genomicFields
   :outertype: GenomicsDBQuery.VariantCall

genomic_interval
^^^^^^^^^^^^^^^^

.. java:field::  Pair genomic_interval
   :outertype: GenomicsDBQuery.VariantCall

rowIndex
^^^^^^^^

.. java:field::  long rowIndex
   :outertype: GenomicsDBQuery.VariantCall

sampleName
^^^^^^^^^^

.. java:field::  String sampleName
   :outertype: GenomicsDBQuery.VariantCall

Constructors
------------
VariantCall
^^^^^^^^^^^

.. java:constructor:: public VariantCall(long rowIndex, long colIndex, String sampleName, String contigName, long start, long end, Map<String, Object> genomicFields)
   :outertype: GenomicsDBQuery.VariantCall

Methods
-------
getColIndex
^^^^^^^^^^^

.. java:method:: public long getColIndex()
   :outertype: GenomicsDBQuery.VariantCall

getContigName
^^^^^^^^^^^^^

.. java:method:: public String getContigName()
   :outertype: GenomicsDBQuery.VariantCall

getGenomicFields
^^^^^^^^^^^^^^^^

.. java:method:: public Map<String, Object> getGenomicFields()
   :outertype: GenomicsDBQuery.VariantCall

getGenomic_interval
^^^^^^^^^^^^^^^^^^^

.. java:method:: public Pair getGenomic_interval()
   :outertype: GenomicsDBQuery.VariantCall

getRowIndex
^^^^^^^^^^^

.. java:method:: public long getRowIndex()
   :outertype: GenomicsDBQuery.VariantCall

getSampleName
^^^^^^^^^^^^^

.. java:method:: public String getSampleName()
   :outertype: GenomicsDBQuery.VariantCall

