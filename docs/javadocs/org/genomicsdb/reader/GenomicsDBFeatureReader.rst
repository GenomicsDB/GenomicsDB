.. java:import:: htsjdk.variant.bcf2 BCF2Codec

.. java:import:: htsjdk.variant.vcf VCFContigHeaderLine

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: org.genomicsdb.model Coordinates

.. java:import:: org.genomicsdb.model GenomicsDBExportConfiguration

.. java:import:: java.io IOException

GenomicsDBFeatureReader
=======================

.. java:package:: org.genomicsdb.reader
   :noindex:

.. java:type:: public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T>

   A reader for GenomicsDB that implements \ :java:ref:`htsjdk.tribble.FeatureReader`\  Currently, the reader only return \ :java:ref:`htsjdk.variant.variantcontext.VariantContext`\

Constructors
------------
GenomicsDBFeatureReader
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBFeatureReader(GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration, FeatureCodec<T, SOURCE> codec, Optional<String> loaderJSONFile) throws IOException
   :outertype: GenomicsDBFeatureReader

   Constructor

   :param exportConfiguration: query parameters
   :param codec: FeatureCodec, currently only \ :java:ref:`htsjdk.variant.bcf2.BCF2Codec`\  and \ :java:ref:`htsjdk.variant.vcf.VCFCodec`\  are tested
   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :throws IOException: when data cannot be read from the stream

Methods
-------
close
^^^^^

.. java:method:: public void close() throws IOException
   :outertype: GenomicsDBFeatureReader

getHeader
^^^^^^^^^

.. java:method:: public Object getHeader()
   :outertype: GenomicsDBFeatureReader

   Return the VCF header of the combined gVCF stream

   :return: the VCF header of the combined gVCF stream

getSequenceNames
^^^^^^^^^^^^^^^^

.. java:method:: public List<String> getSequenceNames()
   :outertype: GenomicsDBFeatureReader

   Return the list of contigs in the combined VCF header

   :return: list of strings of the contig names

iterator
^^^^^^^^

.. java:method:: public CloseableTribbleIterator<T> iterator() throws IOException
   :outertype: GenomicsDBFeatureReader

   Return an iterator over \ :java:ref:`htsjdk.variant.variantcontext.VariantContext`\  objects for the specified TileDB array and query configuration

   :return: iterator over \ :java:ref:`htsjdk.variant.variantcontext.VariantContext`\  objects

query
^^^^^

.. java:method:: public CloseableTribbleIterator<T> query(String chr, int start, int end) throws IOException
   :outertype: GenomicsDBFeatureReader

   Return an iterator over \ :java:ref:`htsjdk.variant.variantcontext.VariantContext`\  objects for the specified TileDB array and queried position

   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :return: iterator over \ :java:ref:`htsjdk.variant.variantcontext.VariantContext`\  objects

