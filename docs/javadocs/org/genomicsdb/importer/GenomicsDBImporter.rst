.. java:import:: htsjdk.samtools.util CloseableIterator

.. java:import:: htsjdk.samtools.util RuntimeIOException

.. java:import:: htsjdk.tribble AbstractFeatureReader

.. java:import:: htsjdk.tribble FeatureReader

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.variantcontext.writer VariantContextWriterBuilder

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: org.genomicsdb Constants

.. java:import:: org.genomicsdb GenomicsDBLibLoader

.. java:import:: org.genomicsdb.exception GenomicsDBException

.. java:import:: org.genomicsdb.importer.extensions CallSetMapExtensions

.. java:import:: org.genomicsdb.importer.extensions JsonFileExtensions

.. java:import:: org.genomicsdb.importer.extensions VidMapExtensions

.. java:import:: org.genomicsdb.importer.model ChromosomeInterval

.. java:import:: org.genomicsdb.importer.model SampleInfo

.. java:import:: org.json.simple JSONArray

.. java:import:: org.json.simple JSONObject

.. java:import:: org.json.simple.parser JSONParser

.. java:import:: org.json.simple.parser ParseException

.. java:import:: java.io File

.. java:import:: java.io FileNotFoundException

.. java:import:: java.io IOException

.. java:import:: java.io StringWriter

.. java:import:: java.util ArrayList

.. java:import:: java.util Iterator

.. java:import:: java.util List

.. java:import:: java.util Map

.. java:import:: java.util.concurrent CompletableFuture

.. java:import:: java.util.concurrent ExecutorService

.. java:import:: java.util.concurrent Executors

.. java:import:: java.util.concurrent ForkJoinPool

.. java:import:: java.util.stream Collectors

.. java:import:: java.util.stream IntStream

GenomicsDBImporter
==================

.. java:package:: org.genomicsdb.importer
   :noindex:

.. java:type:: public class GenomicsDBImporter extends GenomicsDBImporterJni implements JsonFileExtensions, CallSetMapExtensions, VidMapExtensions

   Java wrapper for vcf2genomicsdb - imports VCFs into GenomicsDB. All vid information is assumed to be set correctly by the user (JSON files)

Constructors
------------
GenomicsDBImporter
^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBImporter(String loaderJSONFile)
   :outertype: GenomicsDBImporter

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file

GenomicsDBImporter
^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBImporter(String loaderJSONFile, int rank)
   :outertype: GenomicsDBImporter

   Constructor

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param rank: Rank of this process (TileDB/GenomicsDB partition idx)

GenomicsDBImporter
^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBImporter(ImportConfig config) throws FileNotFoundException, com.googlecode.protobuf.format.JsonFormat.ParseException
   :outertype: GenomicsDBImporter

   Constructor to create required data structures from a list of GVCF files and a chromosome interval. This constructor is developed specifically for running Chromosome intervals imports in parallel.

   :param config: Parallel import configuration
   :throws FileNotFoundException: when files could not be read/written
   :throws com.googlecode.protobuf.format.JsonFormat.ParseException: when existing callset jsons are invalid. incremental case

Methods
-------
add
^^^

.. java:method:: public boolean add(VariantContext vc, int streamIdx) throws GenomicsDBException, RuntimeIOException
   :outertype: GenomicsDBImporter

   Write VariantContext object to stream - may fail if the buffer is full It's the caller's responsibility keep track of the VC object that's not written

   :param vc: VariantContext object
   :param streamIdx: index of the stream returned by the addBufferStream() call
   :return: true if the vc object was written successfully, false otherwise

addBufferStream
^^^^^^^^^^^^^^^

.. java:method:: @SuppressWarnings public int addBufferStream(String streamName, VCFHeader vcfHeader, long bufferCapacity, VariantContextWriterBuilder.OutputType streamType, Iterator<VariantContext> vcIterator, Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException
   :outertype: GenomicsDBImporter

   Add a buffer stream or VC iterator - internal function

   :param streamName: Name of the stream being added - must be unique with respect to this GenomicsDBImporter object
   :param vcfHeader: VCF header for the stream
   :param bufferCapacity: Capacity of the stream buffer in bytes
   :param streamType: BCF_STREAM or VCF_STREAM
   :param vcIterator: Iterator over VariantContext objects - can be null
   :param sampleIndexToInfo: map from sample index in the vcfHeader to SampleInfo object which contains row index and globally unique name can be set to null, which implies that the mapping is stored in a callsets JSON file
   :return: returns the stream index

addSortedVariantContextIterator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public int addSortedVariantContextIterator(String streamName, VCFHeader vcfHeader, Iterator<VariantContext> vcIterator, long bufferCapacity, VariantContextWriterBuilder.OutputType streamType, Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException
   :outertype: GenomicsDBImporter

   Add a sorted VC iterator as the data source - caller must: 1. Call setupGenomicsDBImporter() once all iterators are added 2. Call doSingleImport() 3. Done!

   :param streamName: Name of the stream being added - must be unique with respect to this GenomicsDBImporter object
   :param vcfHeader: VCF header for the stream
   :param vcIterator: Iterator over VariantContext objects
   :param bufferCapacity: Capacity of the stream buffer in bytes
   :param streamType: BCF_STREAM or VCF_STREAM
   :param sampleIndexToInfo: map from sample index in the vcfHeader to SampleInfo object which contains row index and globally unique name can be set to null, which implies that the mapping is stored in a callsets JSON file
   :return: returns the stream index

columnPartitionIterator
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static <SOURCE> MultiChromosomeIterator<SOURCE> columnPartitionIterator(AbstractFeatureReader<VariantContext, SOURCE> reader, String loaderJSONFile, int partitionIdx) throws ParseException, IOException
   :outertype: GenomicsDBImporter

   Utility function that returns a MultiChromosomeIterator given an AbstractFeatureReader that will iterate over the VariantContext objects provided by the reader belonging to the column partition specified by the loader JSON file and rank/partition index

   :param <SOURCE>: LineIterator for VCFs, PositionalBufferedStream for BCFs
   :param reader: AbstractFeatureReader over VariantContext objects - SOURCE can vary - BCF v/s VCF for example
   :param loaderJSONFile: path to loader JSON file
   :param partitionIdx: rank/partition index
   :throws IOException: when the reader's query method throws an IOException
   :throws ParseException: when there is a bug in the JNI interface and a faulty JSON is returned
   :return: MultiChromosomeIterator that iterates over VariantContext objects in the reader belonging to the specified column partition

columnPartitionIterator
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public <SOURCE> MultiChromosomeIterator<SOURCE> columnPartitionIterator(AbstractFeatureReader<VariantContext, SOURCE> reader) throws ParseException, IOException
   :outertype: GenomicsDBImporter

   Utility function that returns a MultiChromosomeIterator given an AbstractFeatureReader that will iterate over the VariantContext objects provided by the reader belonging to the column partition specified by this object's loader JSON file and rank/partition index

   :param <SOURCE>: LineIterator for VCFs, PositionalBufferedStream for BCFs
   :param reader: AbstractFeatureReader over VariantContext objects - SOURCE can vary - BCF v/s VCF for example
   :throws IOException: when the reader's query method throws an IOException
   :throws ParseException: when there is a bug in the JNI interface and a faulty JSON is returned
   :return: MultiChromosomeIterator that iterates over VariantContext objects in the reader belonging to the specified column partition

doSingleImport
^^^^^^^^^^^^^^

.. java:method:: public boolean doSingleImport() throws IOException
   :outertype: GenomicsDBImporter

   Only to be used in cases where iterator of VariantContext are not used. The data is written to buffers directly after which this function is called. See TestBufferStreamGenomicsDBImporter.java for an example

   :throws IOException: if the import fails
   :return: true if the import process is done

executeImport
^^^^^^^^^^^^^

.. java:method:: public void executeImport() throws InterruptedException
   :outertype: GenomicsDBImporter

   Import multiple chromosome interval

   :throws InterruptedException: when there is an exception in any of the threads in the stream

executeImport
^^^^^^^^^^^^^

.. java:method:: public void executeImport(int numThreads)
   :outertype: GenomicsDBImporter

   Import multiple chromosome interval

   :param numThreads: number of threads used to import partitions

getExhaustedBufferStreamIndex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public int getExhaustedBufferStreamIndex(long i)
   :outertype: GenomicsDBImporter

   Get buffer stream index of i-th exhausted stream There are mNumExhaustedBufferStreams and the caller must provide data for streams with indexes getExhaustedBufferStreamIndex(0), getExhaustedBufferStreamIndex(1),..., getExhaustedBufferStreamIndex(mNumExhaustedBufferStreams-1)

   :param i: i-th exhausted buffer stream
   :return: the buffer stream index of the i-th exhausted stream

getNumExhaustedBufferStreams
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public long getNumExhaustedBufferStreams()
   :outertype: GenomicsDBImporter

   :return: get number of buffer streams for which new data must be supplied

getProtobufVidMapping
^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBVidMapProto.VidMappingPB getProtobufVidMapping()
   :outertype: GenomicsDBImporter

   Function to return the vid mapping protobuf object.

   :return: protobuf object for vid mapping

isDone
^^^^^^

.. java:method:: public boolean isDone()
   :outertype: GenomicsDBImporter

   Is the import process completed

   :return: true if complete, false otherwise

setupGenomicsDBImporter
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: @SuppressWarnings public void setupGenomicsDBImporter() throws IOException
   :outertype: GenomicsDBImporter

   Setup the importer after all the buffer streams are added, but before any data is inserted into any stream No more buffer streams can be added once setupGenomicsDBImporter() is called

   :throws IOException: throws IOException if modified callsets JSON cannot be written

updateProtobufVidMapping
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void updateProtobufVidMapping(GenomicsDBVidMapProto.VidMappingPB vidMapPB)
   :outertype: GenomicsDBImporter

   Function to update vid mapping protobuf object in the top level config object. Used in cases where the VCF header doesn't contain accurate information about how to parse fields. For instance, allele specific annotations

   :param vidMapPB: vid mapping protobuf object to use as new

write
^^^^^

.. java:method:: public void write() throws GenomicsDBException
   :outertype: GenomicsDBImporter

   Write to TileDB/GenomicsDB using the configuration specified in the loader file passed to constructor

write
^^^^^

.. java:method:: public void write(int rank) throws GenomicsDBException
   :outertype: GenomicsDBImporter

   Write to TileDB/GenomicsDB using the configuration specified in the loader file passed to constructor

   :param rank: Rank of this process (TileDB/GenomicsDB partition idx)

write
^^^^^

.. java:method:: public void write(String loaderJSONFile, int rank) throws GenomicsDBException
   :outertype: GenomicsDBImporter

   Write to TileDB/GenomicsDB

   :param loaderJSONFile: GenomicsDB loader JSON configuration file
   :param rank: Rank of this process (TileDB/GenomicsDB partition idx)

