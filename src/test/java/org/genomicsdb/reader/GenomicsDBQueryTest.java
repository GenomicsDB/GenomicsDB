/**
 * The MIT License (MIT) Copyright (c) 2019 Omics Data Automation, Inc.
 *
 * <p>Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 * and associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * <p>The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * <p>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package org.genomicsdb.reader;

import org.genomicsdb.exception.GenomicsDBException;
import org.genomicsdb.reader.GenomicsDBQuery.*;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

public class GenomicsDBQueryTest {

  static String inputsDir = Paths.get("target", "test", "inputs").toAbsolutePath().toString();

  static String workspace = Paths.get(inputsDir, "ws").toAbsolutePath().toString();
  static String callsetMapping = Paths.get(inputsDir, "callset_t0_1_2.json").toAbsolutePath().toString();
  static String vidMapping = Paths.get(inputsDir, "vid.json").toAbsolutePath().toString();
  static String referenceGenome = Paths.get(inputsDir, "chr1_10MB.fasta.gz").toAbsolutePath().toString();

  static String queryJSON = Paths.get(inputsDir, "query.json").toAbsolutePath().toString();
  static String loaderJSON = Paths.get(inputsDir, "loader.json").toAbsolutePath().toString();

  static String arrayName = "t0_1_2";

  @BeforeClass
  public void setUp() {
    File inputsDirFile = new File(inputsDir);
    if (!inputsDirFile.exists() || inputsDirFile.list().length == 0) {
      Assert.fail("Aborting GenomicsDBQueryTest. Could not find test inputs folder " + inputsDir);
    }
  }

  @Test
  public void testGenomicsDBQueryVersion() {
    Assert.assertNotNull(new GenomicsDBQuery().version());
  }

  @Test
  void testGenomicsDBConnectWithEmptyRequiredParams() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    try {
      query.connect("", "", "", "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
    try {
      query.connect(workspace, "", "", "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
    try {
      query.connect(workspace, vidMapping, "", "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
    try {
      query.connect(workspace, vidMapping, callsetMapping, "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
  }

  private long connect() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long handle = query.connect(workspace, vidMapping, callsetMapping, referenceGenome, new ArrayList<>());
    Assert.assertTrue(handle > 0);
    return handle;
  }

  private long connectWithDPAttribute() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    List<String> attributes = new ArrayList<>();
    attributes.add("DP");
    long handle = query.connect(workspace, vidMapping, callsetMapping, referenceGenome, attributes);
    Assert.assertTrue(handle > 0);
    return handle;
  }

  private long connectWithSegmentSize() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long handle = query.connect(workspace, vidMapping, callsetMapping, referenceGenome, new ArrayList<>(), 40);
    Assert.assertTrue(handle > 0);
    return handle;
  }

  @Test
  void testGenomicsDBBasicConnectDisconnect() throws GenomicsDBException {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connect();
    query.disconnect(genomicsDBHandle);

    genomicsDBHandle = connectWithDPAttribute();
    query.disconnect(genomicsDBHandle);

    genomicsDBHandle = connectWithSegmentSize();
    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBVariantCallQuery() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithDPAttribute();

    List<VariantCalls> calls = query.queryVariantCalls(genomicsDBHandle, arrayName);
    assert(calls.size() == 0);

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    calls = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges);
    assert(calls.size() == 0);

    List<Pair>rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));
    calls = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    //assert(calls.size() == columnRanges.size());

    query.disconnect(genomicsDBHandle);
  }
/*
  @Test
  void testGenomicsDBVariantCallQueryWithSegmentSize() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    List<Pair>rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    List<VariantCalls> calls = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    assert(calls.size() == 1);

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBVariantCallQueryWithMultipleIntervals() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 50000L));
    columnRanges.add(new Pair(50000L, 1000000000L));
    List<Pair>rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    List<VariantCalls> calls = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    assert(calls.size() == 1);

    query.disconnect(genomicsDBHandle);
  }*/
}
