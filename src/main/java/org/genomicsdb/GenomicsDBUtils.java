/*
 * The MIT License (MIT)
 * Copyright (c) 2018 Omics Data Automation Inc. and Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.genomicsdb;

public class GenomicsDBUtils {
    /**
     * Create TileDB workspace
     *
     * @param workspace path to workspace directory
     * @param replace when set, the directory is deleted first if it exists
     * @return status 0 = workspace created,
     * -1 = path was not a directory,
     * -2 = failed to create workspace,
     * 1 = existing directory, nothing changed
     */
    public static int createTileDBWorkspace(final String workspace, final boolean replace) {
        return GenomicsDBUtilsJni.jniCreateTileDBWorkspace(workspace, replace);
    }

    /**
     * Write contents into file
     * @param filename path to file
     * @param contents buffer to be written out
     * @return status 0 = OK
     */
    public static int writeToFile(final String filename, final String contents) {
	return GenomicsDBUtilsJni.jniWriteToFile(filename, contents, (long)contents.length());
    }

    /**
     * Read entire file as string
     * @param filename path to file
     * @return contents of file as string
     */
    public static String readEntireFile(final String filename) {
	return GenomicsDBUtilsJni.jniReadEntireFile(filename);
    }

    /**
     * Copy source path contents to destination
     * @param source local filesystem path
     * @param destination local or cloud filesystem URI
     * @return status 0 = OK
     */
    public static int moveFile(final String source, final String destination) {
	return GenomicsDBUtilsJni.jniMoveFile(source, destination);
    } 

    /**
     * Delete file
     * @param filename path to file
     * @return status 0 = OK
     */
    public static int deleteFile(final String filename) {
	return GenomicsDBUtilsJni.jniDeleteFile(filename);
    } 

    /**
     * Consolidate TileDB array
     *
     * @param workspace path to workspace directory
     * @param arrayName array name
     */
    public static void consolidateTileDBArray(final String workspace, final String arrayName) {
        GenomicsDBUtilsJni.jniConsolidateTileDBArray(workspace, arrayName);
    }

    /**
     * Checks if GenomicsDB array exists.
     * @param workspace workspace
     * @param arrayName arrayName
     * @return <code>true</code> if workspace with arrayName exists else return <code>false</code>
     */
    public static boolean isGenomicsDBArray(final String workspace, final String arrayName) {
	return GenomicsDBUtilsJni.jniIsTileDBArray(workspace, arrayName);
    }

    /**
     * List the GenomicsDB arrays in the given workspace
     * @param workspace workspace
     * @return names of GenomicsDB arrays if they exist
     */
    public static String[] listGenomicsDBArrays(final String workspace) {
        return GenomicsDBUtilsJni.jniListTileDBArrays(workspace);
    }

    /**
     * List the GenomicsDB fragments in the given workspace
     * @param workspace workspace
     * @return names of GenomicsDB fragments if they exist
     */
    public static String[] listGenomicsDBFragments(final String workspace) {
        return GenomicsDBUtilsJni.jniListTileDBFragments(workspace);
    }

    /**
     * Get max valid row index
     * @param workspace path to workspace
     * @param array name of workspace
     * @return max valid row idx, -1 for error
     */
    public static int getMaxValidRowIndex(final String workspace, final String array) {
	return GenomicsDBUtilsJni.jniGetMaxValidRowIndex(workspace, array);
    } 

}

