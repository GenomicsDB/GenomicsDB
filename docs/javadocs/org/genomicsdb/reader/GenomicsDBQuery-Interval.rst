.. java:import:: org.genomicsdb GenomicsDBLibLoader

.. java:import:: org.genomicsdb.exception GenomicsDBException

GenomicsDBQuery.Interval
========================

.. java:package:: org.genomicsdb.reader
   :noindex:

.. java:type:: public static class Interval
   :outertype: GenomicsDBQuery

Fields
------
calls
^^^^^

.. java:field::  List<VariantCall> calls
   :outertype: GenomicsDBQuery.Interval

interval
^^^^^^^^

.. java:field::  Pair interval
   :outertype: GenomicsDBQuery.Interval

Constructors
------------
Interval
^^^^^^^^

.. java:constructor::  Interval()
   :outertype: GenomicsDBQuery.Interval

Interval
^^^^^^^^

.. java:constructor::  Interval(long start, long end)
   :outertype: GenomicsDBQuery.Interval

Methods
-------
addCall
^^^^^^^

.. java:method:: public void addCall(VariantCall call)
   :outertype: GenomicsDBQuery.Interval

getCalls
^^^^^^^^

.. java:method:: public List<VariantCall> getCalls()
   :outertype: GenomicsDBQuery.Interval

getInterval
^^^^^^^^^^^

.. java:method:: public Pair getInterval()
   :outertype: GenomicsDBQuery.Interval

