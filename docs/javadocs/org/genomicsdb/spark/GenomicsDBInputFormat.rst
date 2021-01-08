.. java:import:: org.genomicsdb.reader GenomicsDBFeatureReader

.. java:import:: org.genomicsdb.model Coordinates

.. java:import:: org.genomicsdb.model GenomicsDBExportConfiguration

.. java:import:: htsjdk.tribble Feature

.. java:import:: htsjdk.tribble FeatureCodec

.. java:import:: htsjdk.variant.bcf2 BCF2Codec

.. java:import:: org.apache.hadoop.conf Configurable

.. java:import:: org.apache.hadoop.conf Configuration

.. java:import:: org.apache.log4j Logger

.. java:import:: org.json.simple JSONObject

.. java:import:: org.json.simple.parser JSONParser

.. java:import:: org.json.simple.parser ParseException

.. java:import:: java.io File

.. java:import:: java.io FileWriter

.. java:import:: java.io FileReader

.. java:import:: java.io FileNotFoundException

.. java:import:: java.io IOException

.. java:import:: java.util ArrayList

.. java:import:: java.util List

.. java:import:: java.util Optional

GenomicsDBInputFormat
=====================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBInputFormat<VCONTEXT extends Feature, SOURCE> extends InputFormat<String, VCONTEXT> implements Configurable

Fields
------
logger
^^^^^^

.. java:field::  Logger logger
   :outertype: GenomicsDBInputFormat

Constructors
------------
GenomicsDBInputFormat
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputFormat()
   :outertype: GenomicsDBInputFormat

   default constructor

GenomicsDBInputFormat
^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBInputFormat(GenomicsDBConfiguration conf)
   :outertype: GenomicsDBInputFormat

Methods
-------
createRecordReader
^^^^^^^^^^^^^^^^^^

.. java:method:: public RecordReader<String, VCONTEXT> createRecordReader(InputSplit inputSplit, TaskAttemptContext taskAttemptContext) throws IOException, InterruptedException
   :outertype: GenomicsDBInputFormat

getConf
^^^^^^^

.. java:method:: @Override public Configuration getConf()
   :outertype: GenomicsDBInputFormat

getSplits
^^^^^^^^^

.. java:method:: @SuppressWarnings public List<InputSplit> getSplits(JobContext jobContext) throws FileNotFoundException
   :outertype: GenomicsDBInputFormat

   When this function is called, it is already assumed that configuration object is set

   :param jobContext: Hadoop Job context passed from newAPIHadoopRDD defined in SparkContext
   :throws FileNotFoundException: Thrown if creaing configuration object fails
   :return: Returns a list of input splits

setConf
^^^^^^^

.. java:method:: @Override public void setConf(Configuration configuration)
   :outertype: GenomicsDBInputFormat

setHostFile
^^^^^^^^^^^

.. java:method:: public GenomicsDBInputFormat<VCONTEXT, SOURCE> setHostFile(String hostFile) throws FileNotFoundException
   :outertype: GenomicsDBInputFormat

   Set the host file path

   :param hostFile: Full qualified path of the hosts file
   :throws FileNotFoundException: thrown if the hosts file is not found
   :return: Returns the same object for forward function calls

setLoaderJsonFile
^^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBInputFormat<VCONTEXT, SOURCE> setLoaderJsonFile(String jsonFile)
   :outertype: GenomicsDBInputFormat

   Set the loader JSON file path

   :param jsonFile: Full qualified path of the loader JSON file
   :return: Returns the same object for forward function calls

setQueryJsonFile
^^^^^^^^^^^^^^^^

.. java:method:: public GenomicsDBInputFormat<VCONTEXT, SOURCE> setQueryJsonFile(String jsonFile)
   :outertype: GenomicsDBInputFormat

   Set the query JSON file path

   :param jsonFile: Full qualified path of the query JSON file
   :return: Returns the same object for forward function calls

