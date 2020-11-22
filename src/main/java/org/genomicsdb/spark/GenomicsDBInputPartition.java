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
import org.apache.spark.sql.connector.read.InputPartition;
import org.apache.spark.sql.types.StructType;

import java.util.Map;
import java.util.ArrayList;

public class GenomicsDBInputPartition implements InputPartition, GenomicsDBInputInterface {

  private String loader;
  private String query;
  private boolean loaderIsPB;
  private boolean queryIsPB;
  private GenomicsDBPartitionInfo partition;
  private ArrayList<GenomicsDBQueryInfo> queryRangeList;
  private String[] hosts;
  private StructType schema;
  private Map<String, GenomicsDBVidSchema> vMap;

  public GenomicsDBInputPartition(){}

  public GenomicsDBInputPartition(String host, GenomicsDBConfiguration gConf, StructType schema) {
    hosts = new String[1];
    hosts[0] = host;
    setLoaderAndQuery(gConf);
    partition = null;
    queryRangeList = null;
    this.schema = schema;
  }

  public GenomicsDBInputPartition(GenomicsDBPartitionInfo partition,
      ArrayList<GenomicsDBQueryInfo> queryList,
      GenomicsDBConfiguration gConf,
      StructType schema) {
    hosts = null;
    setLoaderAndQuery(gConf);
    this.partition = new GenomicsDBPartitionInfo(partition);
    this.queryRangeList = new ArrayList<GenomicsDBQueryInfo>(queryList);
    this.schema = schema;
  }


  @Override
  public String[] preferredLocations() {
    if (hosts == null) {
      return new String[] {};
    }
    return hosts;
  }

  public GenomicsDBPartitionInfo getPartitionInfo() {
    return partition;
  }

  public ArrayList<GenomicsDBQueryInfo> getQueryInfoList() {
    return queryRangeList;
  }

  public StructType getSchema() {
    return schema;
  }

  public Map<String, GenomicsDBVidSchema> getGenomicsDBVidSchema() {
    return vMap;
  }

  public String getLoader() {
    return loader;
  }

  public String getQuery() {
    return query;
  }

  public boolean getLoaderIsPB() {
    return loaderIsPB;
  }

  public boolean getQueryIsPB() {
    return queryIsPB;
  }

  public void setGenomicsDBConf(GenomicsDBConfiguration g) {
    setLoaderAndQuery(g);
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

  public void setQueryInfoList(ArrayList<GenomicsDBQueryInfo> q) {
    queryRangeList = new ArrayList<>(q);
  }

  private void setLoaderAndQuery(GenomicsDBConfiguration g) {
    if (g.get(GenomicsDBConfiguration.LOADERPB) != null) {
      loader = g.get(GenomicsDBConfiguration.LOADERPB);
      loaderIsPB = true;
    }
    else {
      loader = g.get(GenomicsDBConfiguration.LOADERJSON);
      loaderIsPB = false;
    }
    if (g.get(GenomicsDBConfiguration.QUERYPB) != null) {
      query = g.get(GenomicsDBConfiguration.QUERYPB);
      queryIsPB = true;
    }
    else {
      query = g.get(GenomicsDBConfiguration.QUERYJSON);
      queryIsPB = false;
    }
  }
}
