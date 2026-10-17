/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.util.Locale;

public class GenomicsDBLibLoader {
    public final static String GENOMICSDB_LIBRARY_PATH = "genomicsdb.library.path";
    private final static String GENOMICSDB_LIBRARY_NAME = "tiledbgenomicsdb";
    // Separates the parts of a path to a resource in the jar, on every platform
    private final static String RESOURCE_SEPARATOR = "/";
    private static boolean mIsLibraryLoaded = false;

    private static Logger logger = LogManager.getLogger(GenomicsDBLibLoader.class);

    private static native int jniGenomicsDBOneTimeInitialize();

    public static synchronized boolean loadLibrary() {
        if (mIsLibraryLoaded) return true;

        try {
            // If GENOMICSDB_LIBRARY_PATH is specified via -Dgenomicsdb.library.path=/path/to/lib then
            // load library from GENOMICSDB_LIBRARY_PATH
            String genomicsDBLibraryPath =  System.getProperty(GENOMICSDB_LIBRARY_PATH);
            if (genomicsDBLibraryPath != null && !genomicsDBLibraryPath.isEmpty()) {
                File genomicsdbLibraryFile = new File(genomicsDBLibraryPath, System.mapLibraryName(GENOMICSDB_LIBRARY_NAME));
                if (!genomicsdbLibraryFile.exists()) {
                    throw new RuntimeException("GenomicsDB library not found at " + genomicsDBLibraryPath);
                }
                System.load(genomicsdbLibraryFile.getAbsolutePath());
                logger.info("GenomicsDB native library has been loaded from " + genomicsdbLibraryFile.getAbsolutePath());
            } else {
                String libraryName = System.mapLibraryName(GENOMICSDB_LIBRARY_NAME);
                String platform = nativeLibraryPlatform(System.getProperty("os.name"), System.getProperty("os.arch"));
                String platformPath = RESOURCE_SEPARATOR + platform + RESOURCE_SEPARATOR + libraryName;
                // Jars built for a single platform may hold the library at the root instead
                boolean hasPlatformLibrary = GenomicsDBLibLoader.class.getResource(platformPath) != null;
                loadLibraryFromJar(hasPlatformLibrary ? platformPath : RESOURCE_SEPARATOR + libraryName);
            }
        } catch (IOException e) {
            logger.fatal("", e);
            throw new RuntimeException("GenomicsDB native library could not be loaded", e);
        }

        jniGenomicsDBOneTimeInitialize();
        mIsLibraryLoaded = true;
        logger.info("GenomicsDB native library version : " + GenomicsDBUtils.nativeLibraryVersion());
        return true;
    }

    /**
     * Returns the directory inside the jar that holds the native library for a platform, as
     * {@code <os>-<arch>}, e.g. {@code linux-x86_64}, {@code linux-aarch64} or {@code macos-aarch64}.
     * The build writes the library to the same directory.
     *
     * @param osName the {@code os.name} system property
     * @param osArch the {@code os.arch} system property
     */
    static String nativeLibraryPlatform(String osName, String osArch) {
        String os = osName.toLowerCase(Locale.ROOT);
        os = os.startsWith("mac") ? "macos" : os.split(" ")[0];
        String arch = osArch.toLowerCase(Locale.ROOT);
        if (arch.equals("amd64")) {
            arch = "x86_64";
        } else if (arch.equals("arm64")) {
            arch = "aarch64";
        }
        return os + "-" + arch;
    }

    //Copied from http://frommyplayground.com/how-to-load-native-jni-library-from-jar 

    /**
     * Loads library from current JAR archive
     * The file from JAR is copied into system temporary directory and then loaded. The temporary file is deleted after exiting.
     * Method uses String as filename because the pathname is "abstract", not system-dependent.
     *
     * @param path The filename inside JAR as absolute path (beginning with '/'), e.g. /package/File.ext
     * @throws IOException              If temporary file creation or read/write operation fails
     * @throws IllegalArgumentException If source file (param path) does not exist
     * @throws IllegalArgumentException If the path is not absolute or if the filename is shorter than three characters (restriction of @see File#createTempFile(java.lang.String, java.lang.String)).
     */
    private static void loadLibraryFromJar(String path) throws IOException {
        if (!path.startsWith("/")) {
            throw new IllegalArgumentException("The path should be absolute (start with '/').");
        }

        // Obtain filename from path
        String[] parts = path.split("/");
        String filename = (parts.length >= 1) ? parts[parts.length - 1] : null;

        // Split filename to prefix and suffix (extension)
        String prefix = "";
        String suffix = null;
        if (filename != null) {
            parts = filename.split("\\.", 2);
            prefix = parts[0];
            suffix = (parts.length > 1) ? "." + parts[parts.length - 1] : null;
        }

        // Check if the filename is okay
        if (filename == null || prefix.length() < 3) {
            throw new IllegalArgumentException("The filename has to be at least 3 characters long.");
        }

        // Prepare temporary file
        File temp = File.createTempFile(prefix, suffix);
        temp.deleteOnExit();

        if (!temp.exists()) {
            throw new FileNotFoundException("File " + temp.getAbsolutePath() + " does not exist.");
        }

        // Prepare buffer for data copying
        byte[] buffer = new byte[1024];
        int readBytes;

        // Open and check input stream
        InputStream is = GenomicsDBLibLoader.class.getResourceAsStream(path);
        if (is == null) {
          throw new FileNotFoundException("File " + path + " was not found inside JAR.");
        }

        OutputStream os = null;
        try {
          os = new FileOutputStream(temp);
        }
        catch(FileNotFoundException e) {
          is.close();
          throw e;
        }
        try {
          // Open output stream and copy data between source file in JAR and the temporary file
          while ((readBytes = is.read(buffer)) != -1) {
            os.write(buffer, 0, readBytes);
          }
        } finally {
            IOException thrownException = null;
            try {
              os.flush();
            }
            catch(IOException e) {
              thrownException = e;
            }
            // If read/write fails, close streams safely before throwing an exception
            try {
              os.close();
            }
            catch(IOException e) {
              if(thrownException != null)
                thrownException = e;
            }
            is.close();
            if(thrownException != null)
              throw thrownException;
        }

        // Finally, load the library
        System.load(temp.getAbsolutePath());
    }
}
