.. java:import:: org.genomicsdb.importer.model ChromosomeInterval

.. java:import:: htsjdk.samtools SAMSequenceDictionary

.. java:import:: htsjdk.tribble AbstractFeatureReader

.. java:import:: htsjdk.tribble CloseableTribbleIterator

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: java.io IOException

.. java:import:: java.util ArrayList

.. java:import:: java.util Iterator

.. java:import:: java.util List

.. java:import:: java.util NoSuchElementException

MultiChromosomeIterator
=======================

.. java:package:: org.genomicsdb.importer
   :noindex:

.. java:type::  class MultiChromosomeIterator<SOURCE> implements Iterator<VariantContext>

   Given an AbstractFeatureReader over VariantContext objects and a list of ChromosomeInterval objects, this class is an Iterator over VariantContext for all the chromosome intervals in the list

   :param <SOURCE>: LineIterator for VCFs, PositionalBufferedStream for BCFs

Constructors
------------
MultiChromosomeIterator
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor::  MultiChromosomeIterator(AbstractFeatureReader<VariantContext, SOURCE> reader, List<ChromosomeInterval> chromosomeIntervals) throws IOException
   :outertype: MultiChromosomeIterator

   Constructor

   :param reader: AbstractFeatureReader over VariantContext objects - SOURCE can vary - BCF v/s VCF for example
   :param chromosomeIntervals: chromosome intervals over which to iterate
   :throws IOException: when the reader's query method throws an IOException

Methods
-------
hasNext
^^^^^^^

.. java:method:: @Override public boolean hasNext()
   :outertype: MultiChromosomeIterator

next
^^^^

.. java:method:: @Override public VariantContext next() throws NoSuchElementException
   :outertype: MultiChromosomeIterator

