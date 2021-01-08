.. java:import:: org.genomicsdb.exception GenomicsDBException

GenomicsDBUtilsJni
==================

.. java:package:: org.genomicsdb
   :noindex:

.. java:type:: public class GenomicsDBUtilsJni

Methods
-------
jniConsolidateTileDBArray
^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static native void jniConsolidateTileDBArray(String workspace, String arrayName)
   :outertype: GenomicsDBUtilsJni

   Consolidate TileDB array

   :param workspace: path to workspace directory
   :param arrayName: array name

jniCreateTileDBWorkspace
^^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static native int jniCreateTileDBWorkspace(String workspace, boolean replace)
   :outertype: GenomicsDBUtilsJni

   Create TileDB workspace

   :param workspace: path to workspace directory
   :param replace: if set existing directory is first deleted
   :return: status 0 = workspace created, -1 = path was not a directory, -2 = failed to create workspace, 1 = existing directory, nothing changed

jniDeleteFile
^^^^^^^^^^^^^

.. java:method:: public static native int jniDeleteFile(String filename)
   :outertype: GenomicsDBUtilsJni

   Delete file

   :param filename: path to file, can be cloud URL
   :return: status 0 = OK

jniGetMaxValidRowIndex
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static native int jniGetMaxValidRowIndex(String workspace, String array)
   :outertype: GenomicsDBUtilsJni

   Get max valid row index

   :param workspace: path to workspace
   :param array: name of workspace
   :return: max valid row idx, -1 for error

jniIsTileDBArray
^^^^^^^^^^^^^^^^

.. java:method:: public static native boolean jniIsTileDBArray(String workspace, String arrayName)
   :outertype: GenomicsDBUtilsJni

   Checks if GenomicsDB array exists.

   :param workspace: workspace
   :param arrayName: array
   :return: \ ``true``\  if workspace with arrayName exists else return \ ``false``\

jniListTileDBArrays
^^^^^^^^^^^^^^^^^^^

.. java:method:: public static native String[] jniListTileDBArrays(String workspace)
   :outertype: GenomicsDBUtilsJni

   List Arrays in given workspace.

   :param workspace: workspace
   :return: list of arrays in given workspace.

jniListTileDBFragments
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static native String[] jniListTileDBFragments(String workspace)
   :outertype: GenomicsDBUtilsJni

   List Fragments in given workspace.

   :param workspace: workspace
   :return: list of fragments in given workspace.

jniMoveFile
^^^^^^^^^^^

.. java:method:: public static native int jniMoveFile(String source, String destination)
   :outertype: GenomicsDBUtilsJni

   Copy source path contents to destination

   :param source: path to source file, can be cloud URL
   :param destination: path to destination file, can be cloud URL
   :return: status 0 = OK

jniReadEntireFile
^^^^^^^^^^^^^^^^^

.. java:method:: public static native String jniReadEntireFile(String filename)
   :outertype: GenomicsDBUtilsJni

   Read entire file as string

   :param filename: path to file, can be cloud URL
   :return: contents of file as string

jniWriteToFile
^^^^^^^^^^^^^^

.. java:method:: public static native int jniWriteToFile(String filename, String contents, long length)
   :outertype: GenomicsDBUtilsJni

   Write contents into file

   :param filename: path to file, can be cloud URL
   :param contents: buffer to be written out
   :param length: of buffer to be written out
   :return: status 0 = OK

