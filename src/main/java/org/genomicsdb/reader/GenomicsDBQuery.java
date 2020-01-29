/*
 * The MIT License (MIT)
 * Copyright (c) 2019 Omics Data Automation, Inc.
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

package org.genomicsdb.reader;

import org.genomicsdb.GenomicsDBLibLoader;
import org.genomicsdb.exception.GenomicsDBException;

import java.util.*;

public class GenomicsDBQuery {

  static {
    try {
      if (!GenomicsDBLibLoader.loadLibrary()) {
        throw new GenomicsDBException("Could not load genomicsdb native library");
      }
    } catch(UnsatisfiedLinkError ule) {
      throw new GenomicsDBException("Could not load genomicsdb native library", ule);
    }
    try {
      jniInitialize();
    } catch(Exception e) {
      throw new GenomicsDBException("Could not initialize GenomicsDBQuery native bindings");
    }
  }

  public static class Pair {
    private long start;
    private long end;
    public Pair(long start, long end){
      this.start = start;
      this.end = end;
    }
    public long getStart(){
      return start;
    }
    public long getEnd(){
      return end;
    }
  }

  public static class VariantCall {
    long rowIndex;
    long colIndex;
    String sampleName;
    String contigName;
    Pair genomic_interval;
    Map<String, Object> genomicFields;
    public VariantCall(long rowIndex, long colIndex, String sampleName, String contigName, long start, long end, Map<String, Object> genomicFields) {
      this.rowIndex = rowIndex;
      this.colIndex = colIndex;
      this.sampleName = sampleName;
      this.contigName = contigName;
      this.genomic_interval = new Pair(start, end);
      this.genomicFields = genomicFields;
    }
    public long getRowIndex() {
      return rowIndex;
    }
    public long getColIndex() {
      return colIndex;
    }
    public String getSampleName() {
      return sampleName;
    }
    public String getContigName() {
      return contigName;
    }
    public Pair getGenomic_interval() {
      return genomic_interval;
    }
    public Map<String, Object> getGenomicFields() {
      return genomicFields;
    }
  }

  public static class Interval {
    Pair interval = null;
    List<VariantCall> calls = new ArrayList<>();
    Interval() {}
    Interval(long start, long end) {
      interval = new Pair(start, end);
    }
    public List<VariantCall> getCalls() {
      return calls;
    }
    public void addCall(VariantCall call) {
      calls.add(call);
    }
    public Pair getInterval() {
        return interval;
    }
  }

  public String version() {
    return jniVersion();
  }

  public static final long defaultSegmentSize=10*1024*1024;

  public long connect(final String workspace,
                      final String vidMappingFile,
                      final String callsetMappingFile,
                      final String referenceGenome,
                      final List<String> attributes) throws GenomicsDBException {
    return connect(workspace, vidMappingFile, callsetMappingFile, referenceGenome, attributes, defaultSegmentSize);
  }

  public long connect(final String workspace,
                      final String vidMappingFile,
                      final String callsetMappingFile,
                      final String referenceGenome,
                      final List<String> attributes,
                      final long segmentSize) throws GenomicsDBException {
    return jniConnect(workspace, vidMappingFile, callsetMappingFile, referenceGenome, attributes, segmentSize);
  }

  public void disconnect(long handle) {
    jniDisconnect(handle);
  }

  public List<Interval> queryVariantCalls(long handle,
                                          final String arrayName) {
    return jniQueryVariantCalls(handle, arrayName, new ArrayList<>(), new ArrayList<>());
  }

  public List<Interval> queryVariantCalls(long handle,
                                          final String arrayName,
                                          final List<Pair> columnRanges) {
    return jniQueryVariantCalls(handle, arrayName, columnRanges, new ArrayList<>());
  }

  public List<Interval> queryVariantCalls(long handle,
                                          final String arrayName,
                                          final List<Pair> columnRanges,
                                          final List<Pair> rowRanges) {
    return jniQueryVariantCalls(handle, arrayName, columnRanges, rowRanges);
  }

  // Native Bindings
  private static native void jniInitialize();

  private static native String jniVersion();

  private static native long jniConnect(final String workspace,
                                 final String vidMappingFile,
                                 final String callsetMappingFile,
                                 final String referenceGenome,
                                 final List<String> attributes,
                                 final long segmentSize) throws GenomicsDBException;

  private static native void jniDisconnect(long handle);

  private static native List<Interval> jniQueryVariantCalls(long handle,
                                                            final String arrayName,
                                                            final List<Pair> columnRanges,
                                                            final List<Pair> rowRanges);

}
