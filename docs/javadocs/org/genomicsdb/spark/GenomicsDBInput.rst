.. java:import:: org.apache.spark.sql.types StructType

.. java:import:: org.json.simple JSONObject

.. java:import:: org.json.simple JSONArray

.. java:import:: org.json.simple.parser JSONParser

.. java:import:: org.json.simple.parser ParseException

.. java:import:: java.io File

.. java:import:: java.io FileWriter

.. java:import:: java.io FileReader

.. java:import:: java.io FileNotFoundException

.. java:import:: java.io IOException

.. java:import:: java.lang RuntimeException

.. java:import:: java.lang InstantiationException

.. java:import:: java.lang IllegalAccessException

.. java:import:: java.lang ClassNotFoundException

GenomicsDBInput
===============

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBInput<T extends GenomicsDBInputInterface>

   The input class represents all the data being queried from GenomicsDB. This can be used both for hadoop map-reduce style InputSplits and for Spark (datasourcev2) InputPartitions

Fields
------
genomicsDBConfiguration
^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public GenomicsDBConfiguration genomicsDBConfiguration
   :outertype: GenomicsDBInput

Constructors
------------
GenomicsDBInput
^^^^^^^^^^^^^^^

.. java:constructor::  GenomicsDBInput(GenomicsDBConfiguration gdbconf, StructType schema, Map<String, GenomicsDBVidSchema> vMap, long minQBS, long maxQBS, Class<T> clazz)
   :outertype: GenomicsDBInput

   constructor for GenomicsDBInput

   :param gdbconf: GenomicsDBConfiguration object
   :param schema: schema used for the datasource API
   :param vMap: map of attribute to vid mapping information
   :param minQBS: minimum query block size used for partitioning query
   :param maxQBS: maximum query block size used for partitioning query
   :param clazz: Class object used to decide how to instantiate partitions

Methods
-------
createTargetExportConfigurationPB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static GenomicsDBExportConfiguration.ExportConfiguration createTargetExportConfigurationPB(String queryFileOrPB, GenomicsDBPartitionInfo partition, ArrayList<GenomicsDBQueryInfo> queryList, boolean isPB) throws IOException, ParseException
   :outertype: GenomicsDBInput

   Creates export configuration protobuf object based on partition, query and existing query file or protobuf

   :param queryFileOrPB: Existing query json file or base64 encoded protobuf byte data
   :param partition: used to populate array
   :param queryList: used to bound query column ranges
   :param isPB: boolean parameter that denotes if queryFileOrPB is protobuf
   :throws IOException: Thrown if other IO exception while handling file operations
   :throws ParseException: Thrown if JSON parsing fails
   :return: Returns export configuration protobuf object

divideInput
^^^^^^^^^^^

.. java:method:: public List<T> divideInput()
   :outertype: GenomicsDBInput

   Divide input data/datasource into chunks so that we can distribute the work amongst many workers. Called by InputFormat::getSplits and DataSourceReader::planInputPartitions

   :return: Returns list of "work chunks" (either InputSplits or InputPartitions)

getExportConfigurationFromJsonFile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: @Deprecated public static GenomicsDBExportConfiguration.ExportConfiguration.Builder getExportConfigurationFromJsonFile(String queryFile) throws IOException, ParseException
   :outertype: GenomicsDBInput

getGenomicsDBConfiguration
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBConfiguration getGenomicsDBConfiguration()
   :outertype: GenomicsDBInput

   get the GenomicsDBConfiguration object

   :return: GenomicsDBConfiguration object

getSchema
^^^^^^^^^

.. java:method:: public StructType getSchema()
   :outertype: GenomicsDBInput

   Get the schema for DataSource. Only relevant for Spark's Datasourcev2 implementation

   :return: Returns schema associated with DataSourceReaderV2

setGenomicsDBConfiguration
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public void setGenomicsDBConfiguration(GenomicsDBConfiguration gdbconf)
   :outertype: GenomicsDBInput

   set the GenomicsDBConfiguration object

   :param gdbconf: GenomcisDBConfiguration object

