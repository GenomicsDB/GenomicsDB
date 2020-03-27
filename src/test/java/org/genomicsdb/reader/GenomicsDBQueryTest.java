/**
 * The MIT License (MIT) Copyright (c) 2019-2020 Omics Data Automation, Inc.
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
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.reader.GenomicsDBQuery.Interval;
import org.genomicsdb.reader.GenomicsDBQuery.Pair;
import org.genomicsdb.reader.GenomicsDBQuery.VariantCall;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class GenomicsDBQueryTest {

  static String inputsDir = Paths.get("target", "test", "inputs").toAbsolutePath().toString();

  static String workspace = Paths.get(inputsDir, "ws").toAbsolutePath().toString();
  static String callsetMapping = Paths.get(inputsDir, "callset_t0_1_2.json").toAbsolutePath().toString();
  static String vidMapping = Paths.get(inputsDir, "vid.json").toAbsolutePath().toString();
  static String referenceGenome = Paths.get(inputsDir, "chr1_10MB.fasta.gz").toAbsolutePath().toString();

  static String queryJSONFile = Paths.get(inputsDir, "query.json").toAbsolutePath().toString();
  static String loaderJSONFile = Paths.get(inputsDir, "loader.json").toAbsolutePath().toString();

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

  private long connectWithDPandGTAttribute() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    List<String> attributes = new ArrayList<>();
    attributes.add("DP");
    attributes.add("GT");
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

  private long connectWithJson() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    Assert.assertTrue(new File(queryJSONFile).exists());
    Assert.assertTrue(new File(loaderJSONFile).exists());
    long handle = query.connectJSON(queryJSONFile, loaderJSONFile);
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

  void checkVariantCall(VariantCall variantCall) {
    assert(variantCall != null);
    assert(variantCall.getRowIndex() >= 0 && variantCall.getRowIndex() <=2);
    assert(!variantCall.getContigName().isEmpty());
    assert(variantCall.getGenomic_interval().getStart() > 0);
    assert(variantCall.getGenomic_interval().getEnd() > 0);
    assert(variantCall.getGenomicFields().size() > 1);
  }

  @Test
  void testGenomicsDBVariantCallQuery() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithDPAttribute();

    List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle, arrayName);
    assert(intervals.size() == 0);

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges);
    assert(intervals.size() == 0);

    List<Pair>rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));
    intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    assert(intervals.size() == 1);

    Interval interval = intervals.get(0);
    assert(interval.getInterval().getStart() == 0L);
    assert(interval.getInterval().getEnd() == 1000000000L);
    assert(interval.getCalls().size() == 5);

    checkVariantCall(interval.getCalls().get(0));

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBVariantCallQueryWithMultipleAttributes() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithDPandGTAttribute();

    List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle, arrayName);
    Assert.assertEquals(intervals.size(), 0);

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    List<Pair>rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    Assert.assertEquals(intervals.size(), 1);

    Pair interval = intervals.get(0).getInterval();
    List<VariantCall> calls = intervals.get(0).getCalls();

    Assert.assertEquals(interval.getStart(), 0L);
    Assert.assertEquals(interval.getEnd(), 1000000000L);

    Assert.assertEquals(calls.size(), 5);

    boolean foundVariantCall = false;
    for (VariantCall call : calls) {
      if (call.sampleName.equals("HG00141") && call.contigName.equals("1") && call.genomic_interval.getStart() == 12141 && call.genomic_interval.getEnd() == 12295) {
        foundVariantCall = true;
        Assert.assertEquals(call.genomicFields.size(), 3);
        Assert.assertEquals(call.genomicFields.get("REF"), "C");
        Assert.assertEquals(call.genomicFields.get("ALT"), "[<NON_REF>]");
        Assert.assertEquals(call.genomicFields.get("GT"), "0/0");
      }
    }
    Assert.assertTrue(foundVariantCall, "One Variant Call should have been found");
  }
  /* GenomicsDBColumnarField::print_data_in_buffer_at_index  genomicsdb_columnar_field.cc:101*/

  @Test
  void testGenomicsDBVariantCallQueryWithSegmentSize() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    List<Pair>rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    assert(intervals.size() == 1);
    checkVariantCall(intervals.get(0).getCalls().get(0));

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
    
    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBGenerateVCF() throws IOException {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair>columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 50000L));
    columnRanges.add(new Pair(50000L, 1000000000L));
    List<Pair>rowRanges = new ArrayList<>();

    File vcfFile = File.createTempFile("GenomicsDBQueryTest", "");
    File vcfIndexFile = new File(vcfFile+".tbi");

    query.generateVCF(genomicsDBHandle, arrayName, columnRanges, new ArrayList<>(), vcfFile.toString(), "z", true);

    Assert.assertTrue(vcfFile.exists());
    Assert.assertTrue(vcfFile.length() > 0);
    Assert.assertTrue(vcfIndexFile.exists());
    Assert.assertTrue(vcfIndexFile.length() > 0);

    try {
      query.generateVCF(genomicsDBHandle, arrayName, columnRanges, new ArrayList<>(), vcfFile.toString(), "z", false);
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected exception
    }

    vcfFile.delete();
  }

  @Test
  void testGenomicsDBGenerateVCFWithJson() throws IOException {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithJson();

    File vcfFile = File.createTempFile("GenomicsDBQueryTest", "");
    File vcfIndexFile = new File(vcfFile+".tbi");

    query.generateVCF(genomicsDBHandle, vcfFile.toString(), "z", true);

    Assert.assertTrue(vcfFile.exists());
    Assert.assertTrue(vcfFile.length() > 0);
    Assert.assertTrue(vcfIndexFile.exists());
    Assert.assertTrue(vcfIndexFile.length() > 0);
  }

  @Test
  void testGenomicsDBGenerateVCFWithPBExportConfig() throws IOException {
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
            .setWorkspace(workspace)
            .setVidMappingFile(vidMapping)
            .setCallsetMappingFile(callsetMapping)
            .setReferenceGenome(referenceGenome)
            .setArrayName(arrayName)
            .setSegmentSize(40)
            .setScanFull(true)
            .build();

    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);
    Assert.assertTrue(genomicsDBHandle > 0);

    File vcfFile = File.createTempFile("GenomicsDBQueryTest", "");
    File vcfIndexFile = new File(vcfFile+".tbi");

    query.generateVCF(genomicsDBHandle, vcfFile.toString(), "z", true);

    Assert.assertTrue(vcfFile.exists());
    Assert.assertTrue(vcfFile.length() > 0);
    Assert.assertTrue(vcfIndexFile.exists());
    Assert.assertTrue(vcfIndexFile.length() > 0);
  }
}
