/*
 * The MIT License (MIT)
 * Copyright (c) 2021 Omics Data Automation, Inc.
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

import com.google.common.io.Files;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.UUID;

public class GenomicsDBUtilsTest {
  @Test void testGenomicsDBCreateTileDBWorkspace() throws IOException {
    File tmpDir = Files.createTempDir();
    tmpDir.deleteOnExit();
    String workspacePath = String.format("%s/%s", tmpDir.getAbsolutePath(), "workspace");
    Assert.assertEquals(GenomicsDBUtils.listGenomicsDBArrays(workspacePath).length, 0);
    Assert.assertEquals(GenomicsDBUtils.createTileDBWorkspace(workspacePath, false), 0);
    Assert.assertEquals(GenomicsDBUtils.createTileDBWorkspace(workspacePath, false), 1); // Workspace exists nothing changed
    Assert.assertEquals(GenomicsDBUtils.createTileDBWorkspace(workspacePath, true), 0);
    File tmpFile = File.createTempFile("genomicsdb", "file");
    Assert.assertEquals(GenomicsDBUtils.deleteDir(workspacePath), 0);
    Assert.assertEquals(GenomicsDBUtils.deleteDir("non-existent-workspace"), -1);
  }

  boolean isEnvSet(String name) {
    String value = System.getenv(name);
    return value != null && !value.isEmpty();
  }

  String getEnv(String name) {
    return System.getenv(name);
  }

  @Test
  public void testGenomicsDBUseHdfsConnectorEnvCheck() {
    GenomicsDBUtils.useGcsHdfsConnector(true);
    Assert.assertTrue(GenomicsDBUtils.isUseGcsHdfsConnectorSet());
    GenomicsDBUtils.useGcsHdfsConnector(false);
    Assert.assertFalse(GenomicsDBUtils.isUseGcsHdfsConnectorSet());
  }

  @Test
  public void testGenomicDBUseGcsHdfsConnectorOption() throws ClassNotFoundException {
    if (isEnvSet("GOOGLE_APPLICATION_CREDENTIALS") && isEnvSet("GCS_BUCKET_NAME")) {
      String gcs_bucket = getEnv("GCS_BUCKET_NAME");
      String workspace = String.format("gs://%s/genomicsdb_%s_ws", gcs_bucket, UUID.randomUUID());
      Assert.assertEquals(GenomicsDBUtils.createTileDBWorkspace(workspace, false), 0);
      Assert.assertEquals(GenomicsDBUtils.deleteDir(workspace), 0);

      workspace += "_with_connector";
      GenomicsDBUtils.useGcsHdfsConnector(true);
      try {
        Class.forName("com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem");
        Assert.assertEquals(GenomicsDBUtils.createTileDBWorkspace(workspace, false), 0);
        Assert.assertEquals(GenomicsDBUtils.deleteDir(workspace), 0);
      } catch (Exception e) {
        Assert.assertEquals(GenomicsDBUtils.createTileDBWorkspace(workspace, false), -2);
      }
    }
  }
}
