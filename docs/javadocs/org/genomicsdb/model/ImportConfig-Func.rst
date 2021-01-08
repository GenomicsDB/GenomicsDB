.. java:import:: htsjdk.tribble FeatureReader

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.vcf VCFHeaderLine

.. java:import:: java.util.function Function

.. java:import:: java.lang Integer

.. java:import:: java.util.stream IntStream

.. java:import:: java.util.stream Collectors

.. java:import:: java.net URI

ImportConfig.Func
=================

.. java:package:: org.genomicsdb.model
   :noindex:

.. java:type:: @FunctionalInterface public interface Func<T1, T2, T3, R>
   :outertype: ImportConfig

Methods
-------
apply
^^^^^

.. java:method::  R apply(T1 t1, T2 t2, T3 t3)
   :outertype: ImportConfig.Func

