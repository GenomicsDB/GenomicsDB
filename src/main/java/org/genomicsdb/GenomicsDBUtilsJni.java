/*
 * The MIT License (MIT)
 * Copyright (c) 2018 Omics Data Automation Inc. and Intel Corporation
 * Copyright (c) 2021 Omics Data Automation Inc.
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

import org.genomicsdb.exception.GenomicsDBException;

public class GenomicsDBUtilsJni {
    static {
        try {
            boolean loaded = GenomicsDBLibLoader.loadLibrary();
            if (!loaded) throw new GenomicsDBException("Could not load genomicsdb native library");
        } catch (UnsatisfiedLinkError ule) {
            throw new GenomicsDBException("Could not load genomicsdb native library", ule);
        }
    }

    public static native String jniLibraryVersion();

    /**
     * Create TileDB workspace
     *
     * @param workspace path to workspace directory
     * @param replace if set existing directory is first deleted
     * @return status 0 = workspace created,
     * -1 = path was not a directory,
     * -2 = failed to create workspace,
     * 1 = existing directory, nothing changed
     */
    public static native int jniCreateTileDBWorkspace(final String workspace, final boolean replace);

    /**
     * Checks if GenomicsDB array exists.
     * @param workspace workspace
     * @param arrayName array
     * @return <code>true</code> if workspace with arrayName exists else return <code>false</code>
     */
    public static native boolean jniIsTileDBArray(final String workspace, final String arrayName);

    /**
     * List Arrays in given workspace.
     * @param workspace workspace
     * @return list of arrays in given workspace.
     */
    public static native String[] jniListTileDBArrays(final String workspace);

    /**
     * List Fragments in given workspace.
     * @param workspace workspace
     * @return list of fragments in given workspace.
     */
    public static native String[] jniListTileDBFragments(final String workspace);

    /**
     * Write contents into file
     * @param filename path to file, can be cloud URL
     * @param contents buffer to be written out
     * @param length of buffer to be written out
     * @return status 0 = OK
     */
    public static native int jniWriteToFile(final String filename, final String contents, final long length);

    /**
     * Copy source path contents to destination
     * @param source path to source file, can be cloud URL
     * @param destination path to destination file, can be cloud URL
     * @return status 0 = OK
     */
    public static native int jniMoveFile(final String source, final String destination);

    /**
     * Read entire file as string
     * @param filename path to file, can be cloud URL
     * @return contents of file as string
     */
    public static native String jniReadEntireFile(final String filename);

    /**
     * Delete file
     * @param filename path to file, can be cloud URL
     * @return status 0 = OK
     */
    public static native int jniDeleteFile(final String filename);

    /**
     * Delete directory
     * @param dirname path to directory, can be cloud URL
     * @return status 0 = OK
     */
    public static native int jniDeleteDir(final String dirname);

    /**
     * Get max valid row index
     * @param workspace path to workspace
     * @param array name of workspace
     * @return max valid row idx, -1 for error
     */
    public static native int jniGetMaxValidRowIndex(final String workspace, final String array);

    /**
     * Get array column bounds
     * @param workspace path to workspace
     * @param array name of workspace
     * @return array of length 2 with with [0] = min column, [1] = max column
     */
    public static native long[] jniGetArrayColumnBounds(final String workspace, final String array);

    /**
     * The default behavior uses the GCS native SDK client for processing gs:// URI. This can be overridden to use the
     * GCS HDFS Connector instead when the given option is true
     * @param option = true allows use, otherwise disallow use of the GCS HDFS Connector.
     */
    public static native void jniUseGcsHdfsConnector(final boolean option);
}

