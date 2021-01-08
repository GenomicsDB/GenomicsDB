GenomicsDBUtils
===============

.. java:package:: org.genomicsdb
   :noindex:

.. java:type:: public class GenomicsDBUtils

Methods
-------
consolidateTileDBArray
^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static void consolidateTileDBArray(String workspace, String arrayName)
   :outertype: GenomicsDBUtils

   Consolidate TileDB array

   :param workspace: path to workspace directory
   :param arrayName: array name

createTileDBWorkspace
^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static int createTileDBWorkspace(String workspace, boolean replace)
   :outertype: GenomicsDBUtils

   Create TileDB workspace

   :param workspace: path to workspace directory
   :param replace: when set, the directory is deleted first if it exists
   :return: status 0 = workspace created, -1 = path was not a directory, -2 = failed to create workspace, 1 = existing directory, nothing changed

deleteFile
^^^^^^^^^^

.. java:method:: public static int deleteFile(String filename)
   :outertype: GenomicsDBUtils

   Delete file

   :param filename: path to file
   :return: status 0 = OK

getMaxValidRowIndex
^^^^^^^^^^^^^^^^^^^

.. java:method:: public static int getMaxValidRowIndex(String workspace, String array)
   :outertype: GenomicsDBUtils

   Get max valid row index

   :param workspace: path to workspace
   :param array: name of workspace
   :return: max valid row idx, -1 for error

isGenomicsDBArray
^^^^^^^^^^^^^^^^^

.. java:method:: public static boolean isGenomicsDBArray(String workspace, String arrayName)
   :outertype: GenomicsDBUtils

   Checks if GenomicsDB array exists.

   :param workspace: workspace
   :param arrayName: arrayName
   :return: \ ``true``\  if workspace with arrayName exists else return \ ``false``\

listGenomicsDBArrays
^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static String[] listGenomicsDBArrays(String workspace)
   :outertype: GenomicsDBUtils

   List the GenomicsDB arrays in the given workspace

   :param workspace: workspace
   :return: names of GenomicsDB arrays if they exist

listGenomicsDBFragments
^^^^^^^^^^^^^^^^^^^^^^^

.. java:method:: public static String[] listGenomicsDBFragments(String workspace)
   :outertype: GenomicsDBUtils

   List the GenomicsDB fragments in the given workspace

   :param workspace: workspace
   :return: names of GenomicsDB fragments if they exist

moveFile
^^^^^^^^

.. java:method:: public static int moveFile(String source, String destination)
   :outertype: GenomicsDBUtils

   Copy source path contents to destination

   :param source: local filesystem path
   :param destination: local or cloud filesystem URI
   :return: status 0 = OK

readEntireFile
^^^^^^^^^^^^^^

.. java:method:: public static String readEntireFile(String filename)
   :outertype: GenomicsDBUtils

   Read entire file as string

   :param filename: path to file
   :return: contents of file as string

writeToFile
^^^^^^^^^^^

.. java:method:: public static int writeToFile(String filename, String contents)
   :outertype: GenomicsDBUtils

   Write contents into file

   :param filename: path to file
   :param contents: buffer to be written out
   :return: status 0 = OK

