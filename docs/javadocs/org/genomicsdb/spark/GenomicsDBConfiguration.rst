.. java:import:: org.apache.hadoop.conf Configuration

.. java:import:: org.json.simple JSONArray

.. java:import:: org.json.simple JSONObject

.. java:import:: org.json.simple.parser JSONParser

.. java:import:: org.json.simple.parser ParseException

.. java:import:: java.io FileInputStream

.. java:import:: java.io FileNotFoundException

.. java:import:: java.io FileReader

.. java:import:: java.io IOException

.. java:import:: java.io Serializable

.. java:import:: java.util ArrayList

.. java:import:: java.util List

.. java:import:: java.util Scanner

.. java:import:: java.util Iterator

.. java:import:: java.util Base64

.. java:import:: java.lang RuntimeException

GenomicsDBConfiguration
=======================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBConfiguration extends Configuration implements Serializable

   The configuration class enables users to use Java/Scala to populate the input parameters of GenomicsDB. Input parameters can be passed in as json files or base64 encoded protobuf byte data

Fields
------
LOADERJSON
^^^^^^^^^^

.. java:field:: public static final String LOADERJSON
   :outertype: GenomicsDBConfiguration

LOADERPB
^^^^^^^^

.. java:field:: public static final String LOADERPB
   :outertype: GenomicsDBConfiguration

MPIHOSTFILE
^^^^^^^^^^^

.. java:field:: public static final String MPIHOSTFILE
   :outertype: GenomicsDBConfiguration

PARTITION_STRATEGY
^^^^^^^^^^^^^^^^^^

.. java:field:: public static final String PARTITION_STRATEGY
   :outertype: GenomicsDBConfiguration

QUERYJSON
^^^^^^^^^

.. java:field:: @Deprecated public static final String QUERYJSON
   :outertype: GenomicsDBConfiguration

QUERYPB
^^^^^^^

.. java:field:: public static final String QUERYPB
   :outertype: GenomicsDBConfiguration

Constructors
------------
GenomicsDBConfiguration
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBConfiguration()
   :outertype: GenomicsDBConfiguration

GenomicsDBConfiguration
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBConfiguration(Configuration configuration) throws FileNotFoundException
   :outertype: GenomicsDBConfiguration

GenomicsDBConfiguration
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBConfiguration(Configuration configuration, List<GenomicsDBPartitionInfo> list) throws FileNotFoundException
   :outertype: GenomicsDBConfiguration

   Constructor with partition information

   :param configuration: Existing configuration object (can contain Hadoop config values)
   :param list: Contains partition information
   :throws FileNotFoundException: thrown when loader file not found

Methods
-------
getHosts
^^^^^^^^

.. java:method::  List<String> getHosts()
   :outertype: GenomicsDBConfiguration

getPartitions
^^^^^^^^^^^^^

.. java:method::  ArrayList<GenomicsDBPartitionInfo> getPartitions()
   :outertype: GenomicsDBConfiguration

   Return partition list; used when creating input splits.

   :return: Returns ArrayList of PartitionInfo objects

getQueryBlockSize
^^^^^^^^^^^^^^^^^

.. java:method::  long getQueryBlockSize()
   :outertype: GenomicsDBConfiguration

   Return value used to determine optimal query size when creating InputSplits

   :return: Returns QueryBlockSize

getQueryBlockSizeMargin
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method::  long getQueryBlockSizeMargin()
   :outertype: GenomicsDBConfiguration

   Return value used to determine "slop" for optimal query size when creating InputSplits

   :return: Returns QueryBlockSizeMargin

getQueryRanges
^^^^^^^^^^^^^^

.. java:method::  ArrayList<GenomicsDBQueryInfo> getQueryRanges()
   :outertype: GenomicsDBConfiguration

   Return query range list; used when creating input splits.

   :return: Returns ArrayList of QueryRange objects

populateListFromJson
^^^^^^^^^^^^^^^^^^^^

.. java:method::  void populateListFromJson(String jsonType) throws FileNotFoundException, IOException, ParseException
   :outertype: GenomicsDBConfiguration

   parse json file to populate partition and query lists Assuming here that using the Spark interface implies column partitions

   :param jsonType: json file to use while loading - either LOADERJSON or QUERYJSON
   :throws FileNotFoundException: Thrown if queryJson file isn't found
   :throws IOException: Thrown if other IO exception while handling file operations
   :throws ParseException: Thrown if JSON parsing fails

populateListFromPB
^^^^^^^^^^^^^^^^^^

.. java:method::  void populateListFromPB(String pbType)
   :outertype: GenomicsDBConfiguration

   parse protobuf to populate partition and query lists Assuming here that using the Spark interface implies column partitions

   :param pbType: protobuf type to use - either LOADERPB or QUERYPB

setHostFile
^^^^^^^^^^^

.. java:method:: public GenomicsDBConfiguration setHostFile(String path) throws FileNotFoundException
   :outertype: GenomicsDBConfiguration

   Host file contains the hosts where GenomicsDB instances reside. This file can be a replica of the slaves file. For now, we have kept it separate, can be merged later

   :param path: Full path of the host file
   :throws FileNotFoundException: If file not found, throw exception
   :return: GenomicsDBConfiguration object

setLoaderJsonFile
^^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBConfiguration setLoaderJsonFile(String path)
   :outertype: GenomicsDBConfiguration

setQueryJsonFile
^^^^^^^^^^^^^^^^

.. java:method:: @Deprecated public GenomicsDBConfiguration setQueryJsonFile(String path)
   :outertype: GenomicsDBConfiguration

setQueryPB
^^^^^^^^^^

.. java:method:: public GenomicsDBConfiguration setQueryPB(String pb)
   :outertype: GenomicsDBConfiguration

