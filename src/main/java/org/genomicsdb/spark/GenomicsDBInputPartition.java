/*
 * The MIT License (MIT)
 * Copyright (c) 2018 Omics Data Automation
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

import org.apache.spark.sql.catalyst.InternalRow;
import org.apache.spark.sql.sources.v2.reader.InputPartition;
import org.apache.spark.sql.sources.v2.reader.InputPartitionReader;
import org.apache.spark.sql.types.StructType;

import java.util.Map;

public class GenomicsDBInputPartition
    implements InputPartition<InternalRow>, GenomicsDBInputInterface {

  private String loaderFile;
  private GenomicsDBPartitionInfo partition;
  private GenomicsDBQueryInfo queryRange;
  private String queryFile;
  private String[] hosts;
  private StructType schema;
  private Map<String, GenomicsDBVidSchema> vMap;

  public GenomicsDBInputPartition() {}

  public GenomicsDBInputPartition(String host, GenomicsDBConfiguration gConf, StructType schema) {
    hosts = new String[1];
    hosts[0] = host;
    loaderFile = gConf.get(GenomicsDBConfiguration.LOADERJSON);
    queryFile = gConf.get(GenomicsDBConfiguration.QUERYJSON);
    partition = null;
    queryRange = null;
    this.schema = schema;
  }

  public GenomicsDBInputPartition(GenomicsDBPartitionInfo partition,
      GenomicsDBQueryInfo query,
      GenomicsDBConfiguration gConf,
      StructType schema) {
    hosts = null;
    loaderFile = gConf.get(GenomicsDBConfiguration.LOADERJSON);
    queryFile = gConf.get(GenomicsDBConfiguration.QUERYJSON);
    this.partition = new GenomicsDBPartitionInfo(partition);
    this.queryRange = new GenomicsDBQueryInfo(query);
    this.schema = schema;
  }

  public InputPartitionReader createPartitionReader() {
    return new GenomicsDBInputPartitionReader(this);
  }

  public String[] preferredLocations() {
    if (hosts == null) {
      return new String[] {};
    }
    return hosts;
  }

  public GenomicsDBPartitionInfo getPartitionInfo() {
    return partition;
  }

  public GenomicsDBQueryInfo getQueryInfo() {
    return queryRange;
  }

  public StructType getSchema() {
    return schema;
  }

  public Map<String, GenomicsDBVidSchema> getGenomicsDBVidSchema() {
    return vMap;
  }

  public String getLoaderFile() {
    return loaderFile;
  }

  public String getQueryFile() {
    return queryFile;
  }

  public void setGenomicsDBConf(GenomicsDBConfiguration g) {
    loaderFile = g.get(GenomicsDBConfiguration.LOADERJSON);
    queryFile = g.get(GenomicsDBConfiguration.QUERYJSON);
  }

  public void setGenomicsDBSchema(StructType s) {
    schema = s;
  }

  public void setGenomicsDBVidSchema(Map<String, GenomicsDBVidSchema> v) {
    vMap = v;
  }

  public void setPartitionInfo(GenomicsDBPartitionInfo p) {
    partition = p;
  }

  public void setQueryInfo(GenomicsDBQueryInfo q) {
    queryRange = q;
  }
}
