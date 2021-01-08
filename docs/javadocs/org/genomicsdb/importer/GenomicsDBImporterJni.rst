GenomicsDBImporterJni
=====================

.. java:package:: org.genomicsdb.importer
   :noindex:

.. java:type::  class GenomicsDBImporterJni

Methods
-------
jniAddBufferStream
^^^^^^^^^^^^^^^^^^

.. java:method:: native void jniAddBufferStream(long genomicsDBImporterHandle, String streamName, boolean isBCF, long bufferCapacity, byte[] buffer, long numValidBytesInBuffer)
   :outertype: GenomicsDBImporterJni

   Notify importer object that a new stream is to be added

   :param genomicsDBImporterHandle: "pointer" returned by jniInitializeGenomicsDBImporterObject
   :param streamName: name of the stream
   :param isBCF: use BCF format to pass data to C++ layer
   :param bufferCapacity: in bytes
   :param buffer: initialization buffer containing the VCF/BCF header
   :param numValidBytesInBuffer: num valid bytes in the buffer (length of the header)

jniConsolidateTileDBArray
^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: static native void jniConsolidateTileDBArray(String workspace, String arrayName)
   :outertype: GenomicsDBImporterJni

   Consolidate TileDB array

   :param workspace: path to workspace directory
   :param arrayName: array name

jniCopyCallsetMap
^^^^^^^^^^^^^^^^^

.. java:method:: native long jniCopyCallsetMap(long genomicsDBImporterHandle, byte[] callsetMapAsByteArray)
   :outertype: GenomicsDBImporterJni

   Copy the callset map protocol buffer to C++ through JNI

   :param genomicsDBImporterHandle: Reference to a C++ GenomicsDBImporter object
   :param callsetMapAsByteArray: Callset name and row index map
   :return: Reference to a C++ GenomicsDBImporter object as long

jniCopyVidMap
^^^^^^^^^^^^^

.. java:method:: native long jniCopyVidMap(long genomicsDBImporterHandle, byte[] vidMapAsByteArray)
   :outertype: GenomicsDBImporterJni

   Copy the vid map protocol buffer to C++ through JNI

   :param genomicsDBImporterHandle: Reference to a C++ GenomicsDBImporter object
   :param vidMapAsByteArray: INFO, FORMAT, FILTER header lines and contig positions
   :return: Reference to a C++ GenomicsDBImporter object as long

jniGenomicsDBImporter
^^^^^^^^^^^^^^^^^^^^^

.. java:method:: native int jniGenomicsDBImporter(String loaderJSONFile, int rank)
   :outertype: GenomicsDBImporterJni

   Creates GenomicsDBImporter object when importing VCF files (no streams)

   :param loaderJSONFile: Path to loader JSON file
   :param rank: Rank of object - corresponds to the partition index in the loader for which this object will import data
   :return: status - 0 if everything was ok, -1 otherwise

jniGetChromosomeIntervalsForColumnPartition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: static native String jniGetChromosomeIntervalsForColumnPartition(String loaderJSONFile, int rank)
   :outertype: GenomicsDBImporterJni

   Obtain the chromosome intervals for the column partition specified in the loader JSON file identified by the rank. The information is returned as a string in JSON format { "contigs": [ { "chr1": [ 100, 200] }, { "chr2": [ 500, 600] } ] }

   :param loaderJSONFile: path to loader JSON file
   :param rank: rank/partition index
   :return: chromosome intervals for the queried column partition in JSON format

jniImportBatch
^^^^^^^^^^^^^^

.. java:method:: native boolean jniImportBatch(long genomicsDBImporterHandle, long[] exhaustedBufferIdentifiers)
   :outertype: GenomicsDBImporterJni

   Import the next batch of data into TileDB/GenomicsDB

   :param genomicsDBImporterHandle: "pointer" returned by jniInitializeGenomicsDBImporterObject
   :param exhaustedBufferIdentifiers: contains the list of exhausted buffer stream identifiers - the number of exhausted streams is stored in the last element of the array
   :return: true if the whole import process is completed, false otherwise

jniInitializeGenomicsDBImporterObject
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: native long jniInitializeGenomicsDBImporterObject(String loaderJSONFile, int rank)
   :outertype: GenomicsDBImporterJni

   Creates GenomicsDBImporter object when importing VCF files (no streams)

   :param loaderJSONFile: Path to loader JSON file
   :param rank: Rank of object - corresponds to the partition index in the loader for which this object will import data
   :return: "pointer"/address to GenomicsDBImporter object in memory, if 0, then something went wrong

jniSetupGenomicsDBLoader
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: native long jniSetupGenomicsDBLoader(long genomicsDBImporterHandle, String callsetMappingJSON)
   :outertype: GenomicsDBImporterJni

   Setup loader after all the buffer streams are added

   :param genomicsDBImporterHandle: "pointer" returned by jniInitializeGenomicsDBImporterObject
   :param callsetMappingJSON: JSON formatted string containing globally consistent callset name to row index mapping
   :return: maximum number of buffer stream identifiers that can be returned in mExhaustedBufferStreamIdentifiers later (this depends on the number of partitions and the number of buffer streams)

jniWriteDataToBufferStream
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: native void jniWriteDataToBufferStream(long handle, int streamIdx, int partitionIdx, byte[] buffer, long numValidBytesInBuffer)
   :outertype: GenomicsDBImporterJni

   :param handle: "pointer" returned by jniInitializeGenomicsDBImporterObject
   :param streamIdx: stream index
   :param partitionIdx: partition index (unused now)
   :param buffer: buffer containing data
   :param numValidBytesInBuffer: num valid bytes in the buffer

