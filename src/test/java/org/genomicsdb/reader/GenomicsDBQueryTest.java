/**
 * The MIT License (MIT)
 * Copyright (c) 2019-2020 Omics Data Automation, Inc.
 * Copyright (c) 2023 dātma, inc™
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

import com.google.protobuf.util.JsonFormat;
import org.apache.commons.lang3.math.NumberUtils;
import org.genomicsdb.exception.GenomicsDBException;
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.Coordinates.GenomicsDBColumnInterval;
import org.genomicsdb.model.Coordinates.GenomicsDBColumnOrInterval;
import org.genomicsdb.model.Coordinates.TileDBColumnInterval;
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
import java.util.Map;
import java.util.logging.Level;

import static org.genomicsdb.model.GenomicsDBExportConfiguration.*;

public class GenomicsDBQueryTest {

  static String inputsDir = Paths.get("target", "test", "inputs").toAbsolutePath().toString();

  static String workspace;
  static String callsetMapping;
  static String vidMapping;
  static String referenceGenome;

  static String queryJSONFile;
  static String loaderJSONFile;

  static String arrayName = "t0_1_2";

  @BeforeClass
  public void setUp() {
    File inputsDirFile = new File(inputsDir);
    if (!inputsDirFile.exists() || inputsDirFile.list().length == 0) {
      inputsDir = Paths.get("build", "target", "test", "inputs").toAbsolutePath().toString();
      inputsDirFile = new File(inputsDir);
      if (!inputsDirFile.exists() || inputsDirFile.list().length == 0) {
        Assert.fail("Aborting GenomicsDBQueryTest. Could not find test inputs folder " + inputsDir);
      }
    }
    workspace = Paths.get(inputsDir, "ws").toAbsolutePath().toString();
    callsetMapping = Paths.get(inputsDir, "callset_t0_1_2.json").toAbsolutePath().toString();
    vidMapping = Paths.get(inputsDir, "vid.json").toAbsolutePath().toString();
    referenceGenome = Paths.get(inputsDir, "chr1_10MB.fasta.gz").toAbsolutePath().toString();

    queryJSONFile = Paths.get(inputsDir, "query.json").toAbsolutePath().toString();
    loaderJSONFile = Paths.get(inputsDir, "loader.json").toAbsolutePath().toString();
  }

  @Test
  public void testGenomicsDBQueryVersion() {
    Assert.assertNotNull(new GenomicsDBQuery().version());
  }

  @Test
  void testGenomicsDBConnectWithEmptyRequiredParams() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    try {
              query.connect("", "", "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
    try {
      query.connect(workspace, "", "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
    try {
      query.connect(workspace, vidMapping, "", new ArrayList<String>());
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
  }

  private long connect() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long handle = query.connect(workspace, vidMapping, callsetMapping, new ArrayList<>());
    Assert.assertTrue(handle > 0);
    return handle;
  }

  private long connectWithDPAttribute() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    List<String> attributes = new ArrayList<>();
    attributes.add("DP");
    long handle = query.connect(workspace, vidMapping, callsetMapping, attributes);
    Assert.assertTrue(handle > 0);
    return handle;
  }

  private long connectWithDPandGTAttribute() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    List<String> attributes = new ArrayList<>();
    attributes.add("DP");
    attributes.add("GT");
    long handle = query.connect(workspace, vidMapping, callsetMapping, attributes);
    Assert.assertTrue(handle > 0);
    return handle;
  }

  private long connectWithSegmentSize() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long handle = query.connect(workspace, vidMapping, callsetMapping, new ArrayList<>(), 40);
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
    Assert.assertFalse(variantCall.toString().isEmpty());
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

    List<Pair> columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    List<Pair> rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    Assert.assertEquals(intervals.size(), 1);
    Assert.assertFalse(intervals.toString().isEmpty());

    Pair interval = intervals.get(0).getInterval();
    List<VariantCall> calls = intervals.get(0).getCalls();

    Assert.assertEquals(interval.getStart(), 0L);
    Assert.assertEquals(interval.getEnd(), 1000000000L);

    Assert.assertEquals(calls.size(), 5);

    boolean foundVariantCall = false;
    for (VariantCall call : calls) {
      if (call.sampleName.equals("HG00141") && call.contigName.equals("1") && call.genomic_interval.getStart() == 12141
          && call.genomic_interval.getEnd() == 12295) {
        foundVariantCall = true;
        Assert.assertEquals(call.genomicFields.size(), 3);
        Assert.assertEquals(call.genomicFields.get("REF"), "C");
        Assert.assertEquals(call.genomicFields.get("ALT"), "[<NON_REF>]");
        Assert.assertEquals(call.genomicFields.get("GT"), "0/0");
      }
    }
    Assert.assertTrue(foundVariantCall, "One Variant Call should have been found");
  }

  @Test
  void testGenomicsDBVariantCallQueryWithSegmentSize() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair> columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));
    List<Pair> rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    assert (intervals.size() == 1);
    checkVariantCall(intervals.get(0).getCalls().get(0));

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBVariantCallQueryWithMultipleIntervals() {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair> columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 50000L));
    columnRanges.add(new Pair(50000L, 1000000000L));
    List<Pair> rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBVariantCallQueryWithPBExportConfig() throws IOException {
    Coordinates.ContigInterval interval = Coordinates.ContigInterval.newBuilder()
        .setContig("1").setBegin(1).setEnd(100000).build();
    RowRangeList rowRangeList = RowRangeList.newBuilder().addRangeList(RowRange.newBuilder()
            .setLow(0L).setHigh(3L).build()).build();
    ExportConfiguration exportConfiguration = ExportConfiguration.newBuilder()
        .setWorkspace(workspace)
        .setVidMappingFile(vidMapping)
        .setCallsetMappingFile(callsetMapping)
        .setArrayName(arrayName)
        .setSegmentSize(40)
        .addQueryContigIntervals(interval)
        .addQueryRowRanges(rowRangeList)
        .build();

    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);
    Assert.assertTrue(genomicsDBHandle > 0);

    List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle);
    Assert.assertEquals(intervals.size(), 1);
    Assert.assertEquals(intervals.get(0).calls.size(), 5);
    query.disconnect(genomicsDBHandle);

    // With filter
    exportConfiguration = ExportConfiguration.newBuilder()
        .setWorkspace(workspace)
        .setVidMappingFile(vidMapping)
        .setCallsetMappingFile(callsetMapping)
        .setArrayName(arrayName)
        .setQueryFilter("END==17384 && REF==\"G\" && ALT|=\"A\" && GT&=\"0/1\"")
        .setSegmentSize(40)
        .addQueryContigIntervals(interval)
        .addQueryRowRanges(rowRangeList)
        .build();
    query = new GenomicsDBQuery();
    genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);
    Assert.assertTrue(genomicsDBHandle > 0);

    intervals = query.queryVariantCalls(genomicsDBHandle);
    Assert.assertEquals(intervals.size(), 1);
    Assert.assertEquals(intervals.get(0).calls.size(), 2);
    query.disconnect(genomicsDBHandle);

    // With bad filter
    exportConfiguration = ExportConfiguration.newBuilder()
        .setWorkspace(workspace)
        .setVidMappingFile(vidMapping)
        .setCallsetMappingFile(callsetMapping)
        .setArrayName(arrayName)
        .setQueryFilter("POS=17384 && REF==\"G\" && ALT|=\"A\" && GT&=\"0/1\"")
        .setSegmentSize(40)
        .addQueryContigIntervals(interval)
        .addQueryRowRanges(rowRangeList)
        .build();
    query = new GenomicsDBQuery();
    genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);
    Assert.assertTrue(genomicsDBHandle > 0);

    try {
      intervals = query.queryVariantCalls(genomicsDBHandle);
    } catch (GenomicsDBException e) {
      // Expected Exception for bad filter
    }
  }

  @Test
  void testGenomicsDBGenerateVCF() throws IOException {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithSegmentSize();

    List<Pair> columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 50000L));
    columnRanges.add(new Pair(50000L, 1000000000L));
    List<Pair> rowRanges = new ArrayList<>();

    File vcfFile = File.createTempFile("GenomicsDBQueryTest", "");
    File vcfIndexFile = new File(vcfFile + ".tbi");

    query.generateVCF(genomicsDBHandle, arrayName, columnRanges, new ArrayList<>(), referenceGenome, "", vcfFile.toString(), "z", true);

    Assert.assertTrue(vcfFile.exists());
    Assert.assertTrue(vcfFile.length() > 0);
    Assert.assertTrue(vcfIndexFile.exists());
    Assert.assertTrue(vcfIndexFile.length() > 0);

    try {
      query.generateVCF(genomicsDBHandle, arrayName, columnRanges, new ArrayList<>(), referenceGenome, "", vcfFile.toString(), "z", false);
      Assert.fail();
    } catch (GenomicsDBException e) {
      // Expected exception
    }

    vcfFile.delete();

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBGenerateVCFWithJson() throws IOException {
    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = connectWithJson();

    File vcfFile = File.createTempFile("GenomicsDBQueryTest", "");
    File vcfIndexFile = new File(vcfFile + ".tbi");

    query.generateVCF(genomicsDBHandle, vcfFile.toString(), "z", true);

    Assert.assertTrue(vcfFile.exists());
    Assert.assertTrue(vcfFile.length() > 0);
    Assert.assertTrue(vcfIndexFile.exists());
    Assert.assertTrue(vcfIndexFile.length() > 0);

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBGenerateVCFWithPBExportConfig() throws IOException {
    ExportConfiguration exportConfiguration = ExportConfiguration.newBuilder()
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
    File vcfIndexFile = new File(vcfFile + ".tbi");

    query.generateVCF(genomicsDBHandle, vcfFile.toString(), "z", true);

    Assert.assertTrue(vcfFile.exists());
    Assert.assertTrue(vcfFile.length() > 0);
    Assert.assertTrue(vcfIndexFile.exists());
    Assert.assertTrue(vcfIndexFile.length() > 0);

    query.disconnect(genomicsDBHandle);
  }

  @Test
  void testGenomicsDBDemoWorkspace() throws IOException {
    String ws = System.getenv("GENOMICSDB_DEMO_WS");
    if (ws == null || ws.isEmpty()) {
      return;
    }
    String arrayName = "allcontigs$1$3095677412";
    Coordinates.ContigInterval interval = Coordinates.ContigInterval.newBuilder()
            .setContig("17").setBegin(7571719).setEnd(7590868).build();
    RowRangeList rowRangeList = RowRangeList.newBuilder().addRangeList(RowRange.newBuilder()
            .setLow(0l).setHigh(200000l).build()).build();

    String filters[] = {"", "REF==\"A\"", "REF==\"A\" && ALT|=\"T\"", "REF==\"A\" && ALT|=\"T\" && GT&=\"1/1\""};
    Long expected_calls[] = new Long[]{ 2962039L, 400032L, 82245L, 82245L };
    Assert.assertEquals(filters.length, expected_calls.length);

    for (int i=0; i< filters.length; i++) {
      Long startTime = System.currentTimeMillis();
      ExportConfiguration exportConfiguration = ExportConfiguration.newBuilder()
              .setWorkspace(ws)
              .setVidMappingFile(ws + "/vidmap.json")
              .setCallsetMappingFile(ws + "/callset.json")
              .setArrayName(arrayName)
              .setEnableSharedPosixfsOptimizations(true)
              .setBypassIntersectingIntervalsPhase(true)
              .setQueryFilter(filters[i])
              .addQueryContigIntervals(interval)
              .addQueryRowRanges(rowRangeList)
              .build();
      GenomicsDBQuery query = new GenomicsDBQuery();
      long genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);
      Assert.assertTrue(genomicsDBHandle > 0);
      List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle);
      Assert.assertEquals(intervals.size(), 1);
      Assert.assertEquals(intervals.get(0).calls.size(), expected_calls[i]);
      System.out.println("Elapsed Time for "+filters[i]+(System.currentTimeMillis()-startTime)/1000+"s");
      query.disconnect(genomicsDBHandle);
    }
  }

  /**
   * Test the VCF annotation service. This function is based on the java
   * testGenomicsDBGenerateVCFWithPBExportConfig test and the
   * annotate_variant_calls_with_tds0 in test_genomicsdb_api.cc.
   */
  @Test
  void testGenomicsDbVcfAnnotationService() {
    String dataSourceName = "dataSourceZero";

    ExportConfiguration exportConfiguration = ExportConfiguration
        .newBuilder()
        .setWorkspace(workspace)
        .setVidMappingFile(vidMapping)
        .setCallsetMappingFile(callsetMapping)
        .setReferenceGenome(referenceGenome)
        .setArrayName(arrayName)
        .setSegmentSize(40)
        .addQueryColumnRanges(GenomicsDBColumnOrIntervalList.newBuilder()
            .addColumnOrIntervalList(GenomicsDBColumnOrInterval.newBuilder()
                .setColumnInterval(GenomicsDBColumnInterval.newBuilder()
                    .setTiledbColumnInterval(TileDBColumnInterval.newBuilder()
                        .setBegin(0)
                        .setEnd(1000000000)))))
        .addAnnotationSource(AnnotationSource.newBuilder()
            .setDataSource(dataSourceName)
            .setFilename(inputsDir + "/test_datasource0.vcf.bgz")
            .addAttributes("ID")
            .addAttributes("field0")
            .addAttributes("field1")
            .addAttributes("field2")
            .addAttributes("field3")
            .addAttributes("field4")
            .addAttributes("field5")
            .addAttributes("field6")
            .addAttributes("field7")
            .setIsVcf(true))
        .build();

    GenomicsDBQuery query = new GenomicsDBQuery();
    long genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);
    Assert.assertTrue(genomicsDBHandle > 0);

    List<Interval> intervals = query.queryVariantCalls(genomicsDBHandle, arrayName);
    assert (intervals.size() == 0);

    List<Pair> columnRanges = new ArrayList<>();
    columnRanges.add(new Pair(0L, 1000000000L));

    List<Pair> rowRanges = new ArrayList<>();
    rowRanges.add(new Pair(0L, 3L));

    intervals = query.queryVariantCalls(genomicsDBHandle, arrayName, columnRanges, rowRanges);
    assert (intervals.size() == 2);

    int variantsAnnotated = 0;

    for (Interval interval : intervals) {
      for (VariantCall call : interval.getCalls()) {
        if (call.getSampleName().equals("HG00141")
            && call.getContigName().equals("1")
            && call.getGenomic_interval().getStart() == 17385
            && String.valueOf(call.getGenomicFields().get("REF")).equals("G")
            && getVariantAlts((String) call.getGenomicFields().get("ALT"))[0].equals("A")) {
          verifyAnnotations(call.getGenomicFields(), dataSourceName, "id001",
              2345, 5678, 3.141500, "noot", null, null, null, true);
          ++variantsAnnotated;
        } else if (call.getSampleName().equals("HG01958")
            && call.getContigName().equals("1")
            && call.getGenomic_interval().getStart() == 17385
            && String.valueOf(call.getGenomicFields().get("REF")).equals("G")
            && getVariantAlts((String) call.getGenomicFields().get("ALT"))[0].equals("T")) {

          verifyAnnotations(call.getGenomicFields(), dataSourceName, "id002", null, null, null, "waldo", 5679,
              6.660000, null, false);
          ++variantsAnnotated;
        } else if (call.getSampleName().equals("HG01530")
            && call.getContigName().equals("1")
            && call.getGenomic_interval().getStart() == 17385
            && String.valueOf(call.getGenomicFields().get("REF")).equals("G")
            && getVariantAlts((String) call.getGenomicFields().get("ALT"))[0].equals("A")) {

          verifyAnnotations(call.getGenomicFields(), dataSourceName, "id001",
              2345, 5678, 3.141500, "noot", null, null, null, true);
          ++variantsAnnotated;
        }

      }
    }

    Assert.assertEquals(variantsAnnotated, 6, "Expected six variants to have annotations");

    query.disconnect(genomicsDBHandle);
  }

  /**
   * Compare the annotations in the genomicFields map to the expected values
   *
   * @param genomicFields
   * @param field0
   * @param field1
   * @param field2
   * @param field3
   * @param field4
   * @param field5
   * @param field6
   * @param field7
   */
  private void verifyAnnotations(Map<String, Object> genomicFields, String dataSourceName, String fieldId,
      Integer field0, Integer field1, Double field2,
      String field3, Integer field4, Double field5, String field6, boolean field7) {

    if (fieldId != null) {
      Assert.assertEquals(String.valueOf(genomicFields.get(dataSourceName + "_ID")), fieldId);
    }

    if (field0 != null) {
      Assert.assertEquals(NumberUtils.toInt((String) genomicFields.get(dataSourceName + "_field0")), field0.intValue());
    }

    if (field1 != null) {
      Assert.assertEquals(NumberUtils.toInt((String) genomicFields.get(dataSourceName + "_field1")), field1.intValue());
    }

    if (field2 != null) {
      Assert.assertEquals(NumberUtils.toDouble((String) genomicFields.get(dataSourceName + "_field2")),
          field2.doubleValue());
    }

    if (field3 != null) {
      Assert.assertEquals(String.valueOf(genomicFields.get(dataSourceName + "_field3")), field3);
    }

    if (field4 != null) {
      Assert.assertEquals(NumberUtils.toInt((String) genomicFields.get(dataSourceName + "_field4")), field4.intValue());
    }

    if (field5 != null) {
      Assert.assertEquals(NumberUtils.toDouble((String) genomicFields.get(dataSourceName + "_field5")),
          field5.doubleValue());
    }

    if (field6 != null) {
      Assert.assertEquals(String.valueOf(genomicFields.get(dataSourceName + "_field6")), field6);
    }

    Assert.assertEquals(genomicFields.containsKey(dataSourceName + "_field7"), field7);

  }

  /**
   * Convert the comma separated GenomicFields["ALT"] value to an array of strings
   *
   * @param alt
   * @return
   */
  private String[] getVariantAlts(String alt) {
    if (alt == null) {
      return null;
    } else {
      return alt.substring(1, alt.length() - 1).split(",");
    }
  }

  /**
   * Expect GenomicsDBException to be thrown when
   */
  @Test
  void testGenomicsDbVcfAnnotationServiceMissingVcfException() {
    String dataSourceName = "dataSource123";

    ExportConfiguration exportConfiguration = ExportConfiguration
        .newBuilder()
        .setWorkspace(workspace)
        .setVidMappingFile(vidMapping)
        .setCallsetMappingFile(callsetMapping)
        .setReferenceGenome(referenceGenome)
        .setArrayName(arrayName)
        .setSegmentSize(40)
        .addQueryColumnRanges(GenomicsDBColumnOrIntervalList.newBuilder()
            .addColumnOrIntervalList(GenomicsDBColumnOrInterval.newBuilder()
                .setColumnInterval(GenomicsDBColumnInterval.newBuilder()
                    .setTiledbColumnInterval(TileDBColumnInterval.newBuilder()
                        .setBegin(0)
                        .setEnd(1000000000)))))
        .addAnnotationSource(AnnotationSource.newBuilder()
            .setDataSource(dataSourceName)
            .setFilename(inputsDir + "/file_does_not_exist.vcf.bgz")
            .setIsVcf(true))
        .build();

    GenomicsDBQuery query = new GenomicsDBQuery();

    try {
      query.connectExportConfiguration(exportConfiguration);
      Assert.fail("Expected GenomicsDBException due to missing VCF");
    } catch (GenomicsDBException e) {
      // Expected Exception
    }
  }
}
