/*
 * The MIT License (MIT)
 * Copyright (c) 2018 University of California, Los Angeles and Intel Corporation
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
package org.genomicsdb.spark;

import org.genomicsdb.GenomicsDBTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.util.ReflectionUtils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.InterruptedException;
import java.util.ArrayList;
import java.util.List;

public class GenomicsDBInputFormatTest {

  @Test(testName = "Testcase0 for creating InputSplits",
      dataProvider = "loaderQueryHostFilesTest0",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testGetSplits0(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException{
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBPartitionInfo> pList = new ArrayList<>(3);
    ArrayList<GenomicsDBQueryInfo> qList = new ArrayList<>(3);
    // test with queries that will not get split up and each
    // query maps to a separate partition
    for(int i=0; i<3; i++) {
      GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(i*500000, "hdfs://tmp/ws", "part"+i, "/tmp/test0.vcf.gz");
      GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(i*500000+1000, i*500000+1100);
      pList.add(p);
      qList.add(q);
    }

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<GenomicsDBInputSplit> splits = format.getSplits(job);
    Assert.assertEquals(splits.size(), 3);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), pList.get(i));
      Assert.assertEquals(gSplit.getQueryInfo(), qList.get(i));
    }
  }

  @Test(testName = "Test each query maps to same partition, no split",
      dataProvider = "loaderQueryHostFilesTest1",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testAllQueriesMapsToSamePartitionNoSplit(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException{
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBQueryInfo> qList = new ArrayList<>(3);
    // test with queries that will not get split up and all
    // query maps to a single partition
    GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(0, "hdfs://tmp/ws", "part", "/tmp/test0.vcf.gz");
    for(int i=0; i<3; i++) {
      GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(i*500000+2000, i*500000+2100);
      qList.add(q);
    }

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    Assert.assertEquals(splits.size(), 3);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), p);
      Assert.assertEquals(gSplit.getQueryInfo(), qList.get(i));
    }
  }

  @Test(testName = "Test single query splits up by query block over single partition",
      dataProvider = "loaderQueryHostFilesTest2",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testQuerySplitUpSinglePartition(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException {
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    // test with single query that will get split up and maps to single partition
    GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(0, "hdfs://tmp/ws", "part", "/tmp/test0.vcf.gz");
    int qstart = 500;
    int qend = 25000;
    GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(qstart, qend);

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    Assert.assertEquals(splits.size(), 3);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), p);
      if (i==0) {
        Assert.assertEquals(gSplit.getQueryInfo().getBeginPosition(), qstart);
      }
      else {
        GenomicsDBInputSplit gSplitPrev = (GenomicsDBInputSplit)splits.get(i-1);
        Assert.assertEquals(gSplit.getQueryInfo().getBeginPosition(), 
			gSplitPrev.getQueryInfo().getEndPosition()+1);
      }
      if (i==splits.size()-1) {
        Assert.assertEquals(gSplit.getQueryInfo().getEndPosition(), qend);
      }
    }
  }

  @Test(testName = "Test local filesystem splits same as hdfs-compliant",
      dataProvider = "loaderQueryHostFilesTest3",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testGetSplits3(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException {
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBPartitionInfo> pList = new ArrayList<>(3);
    // test with single query and non-hdfs compliant data store
    GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(0, "/tmp/ws", "part", "/tmp/test0.vcf.gz");
    int qstart = 500;
    int qend = 25000;
    GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(qstart, qend);

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    // local disk will also split in the same way as hdfs compliant stores now
    Assert.assertEquals(splits.size(), 1);
    GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(0);
    Assert.assertEquals(gSplit.getPartitionInfo().getBeginPosition(), 0);
    Assert.assertEquals(gSplit.getQueryInfo().getBeginPosition(), 500);
    Assert.assertEquals(gSplit.getQueryInfo().getEndPosition(), 25000);
  }

  @Test(testName = "Test query larger than partitions",
      dataProvider = "loaderQueryHostFilesTest4",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testQueryLargerThanPartitions(String queryPath, 
              String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException {
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBPartitionInfo> pList = new ArrayList<>(3);
    // test with single query that will get split up and maps to single partition
    for(int i=0; i<3; i++) {
      GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(i*10000, "hdfs://tmp/ws", "part"+i, "/tmp/test0.vcf.gz");
      pList.add(p);
    }
    int qstart = 500;
    int qend = 25000;
    GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(qstart, qend);

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    Assert.assertEquals(splits.size(), 3);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), pList.get(i));
    }
  }

}
