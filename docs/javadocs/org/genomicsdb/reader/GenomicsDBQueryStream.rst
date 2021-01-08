.. java:import:: org.genomicsdb GenomicsDBLibLoader

.. java:import:: org.genomicsdb.exception GenomicsDBException

.. java:import:: org.genomicsdb.model GenomicsDBExportConfiguration

.. java:import:: java.io IOException

.. java:import:: java.io InputStream

GenomicsDBQueryStream
=====================

.. java:package:: org.genomicsdb.reader
   :noindex:

.. java:type:: public class GenomicsDBQueryStream extends InputStream

   Provides a java.io.InputStream interface for the GenomicsDB combine gVCF operation. It can be used as to construct a \ ` PositionalBufferedStream <https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/tribble/readers/PositionalBufferedStream.html>`_\  object.The PositionalBufferedStream object can then be used by FeatureCodecs such as BCF2Codec to construct VariantContext objects

Constructors
------------
GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, boolean readAsBCF) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param readAsBCF: serialize-deserialize VCF records as BCF2 records
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, boolean readAsBCF, boolean produceHeaderOnly) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param readAsBCF: serialize-deserialize VCF records as BCF2 records
   :param produceHeaderOnly: produce only the header
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, String chr, int start, int end) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, String chr, int start, int end, boolean readAsBCF) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :param readAsBCF: serialize-deserialize VCF records as BCF2 records
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, String chr, int start, int end, int rank) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :param rank: rank of this object if launched from within an MPI context (not used)
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, String chr, int start, int end, int rank, long bufferCapacity, long segmentSize) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :param rank: rank of this object if launched from within an MPI context (not used)
   :param bufferCapacity: size of buffer in bytes to be used by the native layer to store combined BCF2 records
   :param segmentSize: buffer to be used for querying TileDB
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, String chr, int start, int end, int rank, long bufferCapacity, long segmentSize, boolean readAsBCF, boolean produceHeaderOnly) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :param rank: rank of this object if launched from within an MPI context (not used)
   :param bufferCapacity: size of buffer in bytes to be used by the native layer to store combined BCF2 records
   :param segmentSize: buffer to be used for querying TileDB
   :param readAsBCF: serialize-deserialize VCF records as BCF2 records
   :param produceHeaderOnly: produce VCF/BCF header only - no records (minor optimization)
   :throws IOException: when data cannot be read from the stream

GenomicsDBQueryStream
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBQueryStream(String loaderJSONFile, GenomicsDBExportConfiguration.ExportConfiguration queryPB, String chr, int start, int end, int rank, long bufferCapacity, long segmentSize, boolean readAsBCF, boolean produceHeaderOnly, boolean useMissingValuesOnlyNotVectorEnd, boolean keepIDXFieldsInHeader) throws IOException
   :outertype: GenomicsDBQueryStream

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param queryPB: GenomicsDB query protobuf
   :param chr: contig name
   :param start: start position (1-based)
   :param end: end position, inclusive (1-based)
   :param rank: rank of this object if launched from within an MPI context (not used)
   :param bufferCapacity: size of buffer in bytes to be used by the native layer to store combined BCF2 records
   :param segmentSize: buffer to be used for querying TileDB
   :param readAsBCF: serialize-deserialize VCF records as BCF2 records
   :param produceHeaderOnly: produce VCF/BCF header only - no records (minor optimization)
   :param useMissingValuesOnlyNotVectorEnd: don't add BCF2.2 vector end values
   :param keepIDXFieldsInHeader: keep BCF IDX fields in header
   :throws IOException: when data cannot be read from the stream

Methods
-------
available
^^^^^^^^^

.. java:method:: @Override public int available() throws IOException
   :outertype: GenomicsDBQueryStream

close
^^^^^

.. java:method:: @Override public void close() throws IOException
   :outertype: GenomicsDBQueryStream

markSupported
^^^^^^^^^^^^^

.. java:method:: @Override public boolean markSupported()
   :outertype: GenomicsDBQueryStream

read
^^^^

.. java:method:: @Override public int read() throws IOException
   :outertype: GenomicsDBQueryStream

read
^^^^

.. java:method:: @Override public int read(byte[] buffer, int off, int len) throws IOException
   :outertype: GenomicsDBQueryStream

skip
^^^^

.. java:method:: @Override public long skip(long n) throws IOException
   :outertype: GenomicsDBQueryStream

