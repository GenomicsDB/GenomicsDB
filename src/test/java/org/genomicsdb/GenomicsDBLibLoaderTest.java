/*
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation, Inc.
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
import org.junit.internal.runners.statements.Fail;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.stream.Collectors;

public class GenomicsDBLibLoaderTest {
  @Test
  void testGenomicsDBLoader() {
    Assert.assertTrue(GenomicsDBLibLoader.loadLibrary());
  }

  public int runTestInSeparateProcess(String libraryPath) throws Exception {
    String classpath = Arrays.stream(((URLClassLoader) Thread.currentThread().getContextClassLoader()).getURLs())
            .map(URL::getFile)
            .collect(Collectors.joining(File.pathSeparator));
    ProcessBuilder processBuilder = null;
    if (libraryPath == null) {
      processBuilder = new ProcessBuilder(System.getProperty("java.home") + "/bin/java",
                      "-classpath", System.getProperty("java.class.path"),
                      GenomicsDBLibLoaderTest.class.getName());
    } else {
      processBuilder = new ProcessBuilder(System.getProperty("java.home") + "/bin/java",
              "-D"+GenomicsDBLibLoader.GENOMICSDB_LIBRARY_PATH+"="+libraryPath,
              "-classpath", System.getProperty("java.class.path"),
              GenomicsDBLibLoaderTest.class.getName());
    }
    Process process = processBuilder.inheritIO().start();
    return process.waitFor();
  }

  @Test
  public void testGenomicsDBLibLoaderInSeparateProcess() throws Exception {
    Assert.assertEquals(runTestInSeparateProcess(null), 0);
  }

  @Test
  public void testGenomicsDBLibLoaderFromPath() throws Exception {
    String buildDir = Paths.get("target", "classes").toAbsolutePath().toString();
    if (!new File(buildDir).exists()) {
      buildDir = Paths.get("build", "target", "classes").toAbsolutePath().toString();
    }
    if (new File(buildDir).exists()) {
      Assert.assertEquals(runTestInSeparateProcess(buildDir), 0);
    }
  }

  @Test
  public void testGenomicsDBLibLoaderFromNonExistentPath() {
    String buildDir = "non-existent-path";
    Assert.assertFalse(new File(buildDir).exists());
    try {
      Assert.assertEquals(runTestInSeparateProcess(buildDir), 1);
    } catch (Exception e) {
      // Pass
    }
  }

  public static void main(String[] args) {
    if (!GenomicsDBLibLoader.loadLibrary()) {
      throw new RuntimeException("GenomicsDB Native Library could not be loaded");
    }
  }
}
