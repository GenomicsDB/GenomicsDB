.. java:import:: org.genomicsdb.exception GenomicsDBException

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.variantcontext.writer Options

.. java:import:: htsjdk.variant.variantcontext.writer VariantContextWriter

.. java:import:: htsjdk.variant.variantcontext.writer VariantContextWriterBuilder

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: java.util Iterator

GenomicsDBImporterStreamWrapper
===============================

.. java:package:: org.genomicsdb.importer
   :noindex:

.. java:type::  class GenomicsDBImporterStreamWrapper

   Utility class wrapping a stream and a VariantContextWriter for a given stream Each GenomicsDB import stream consists of a buffer stream and a writer object If the caller provides an iterator, then currentVC points to the VariantContext object to be written, if any.

Constructors
------------
GenomicsDBImporterStreamWrapper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBImporterStreamWrapper(VCFHeader vcfHeader, long bufferCapacity, VariantContextWriterBuilder.OutputType streamType, Iterator<VariantContext> vcIterator) throws GenomicsDBException
   :outertype: GenomicsDBImporterStreamWrapper

   Constructor

   :param vcfHeader: VCF header for the stream
   :param bufferCapacity: Capacity of the stream buffer in bytes
   :param streamType: BCF_STREAM or VCF_STREAM
   :param vcIterator: iterator over VariantContext objects, can be null if the caller is managing the buffer explicitly

Methods
-------
getCurrentVC
^^^^^^^^^^^^

.. java:method::  VariantContext getCurrentVC()
   :outertype: GenomicsDBImporterStreamWrapper

   Returns currentVC - could be null if iterator is null or !iterator.hasNext()

   :return: VariantContext object to be written

getStream
^^^^^^^^^

.. java:method:: public SilentByteBufferStream getStream()
   :outertype: GenomicsDBImporterStreamWrapper

getVcWriter
^^^^^^^^^^^

.. java:method:: public VariantContextWriter getVcWriter()
   :outertype: GenomicsDBImporterStreamWrapper

hasIterator
^^^^^^^^^^^

.. java:method::  boolean hasIterator()
   :outertype: GenomicsDBImporterStreamWrapper

   Returns true if a non-null Iterator over VariantContext objects was provided for this stream

   :return: True if iterator was provided, False otherwise

next
^^^^

.. java:method:: public VariantContext next()
   :outertype: GenomicsDBImporterStreamWrapper

   Returns the next VariantContext object iff the Iterator over VariantContext objects is non-null and has a next() object, else returns null. Stores the result in currentVC

   :return: the next VariantContext object or null

