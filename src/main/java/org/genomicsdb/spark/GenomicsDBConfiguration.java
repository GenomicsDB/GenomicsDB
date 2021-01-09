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

package org.genomicsdb.spark;

import org.apache.hadoop.conf.Configuration;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import org.genomicsdb.model.*;

import scala.collection.JavaConverters;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Iterator;
import java.util.Base64;
import java.lang.RuntimeException;
import java.io.FileNotFoundException;

/**
 * The configuration class enables users to use Java/Scala
 * to populate the input parameters of GenomicsDB. 
 * Input parameters can be passed in as json files or 
 * base64 encoded protobuf byte data
 */
public class GenomicsDBConfiguration extends Configuration implements Serializable {

  public static final String LOADERJSON = "genomicsdb.input.loaderjsonfile";
  @Deprecated
  public static final String QUERYJSON = "genomicsdb.input.queryjsonfile";
  public static final String MPIHOSTFILE = "genomicsdb.input.mpi.hostfile";
  public static final String PARTITION_STRATEGY = "genomicsdb.partition.strategy";
  public static final String LOADERPB = "genomicsdb.input.loaderprotobuf";
  public static final String QUERYPB = "genomicsdb.input.queryprotobuf";

  private Boolean produceCombinedVCF = false;
  private Boolean produceTileDBArray = false;
  private Integer segmentSize = 1000;
  private Integer nCellsPerTile = 1000;

  private ArrayList<GenomicsDBPartitionInfo> partitionInfoList = null;
  private ArrayList<GenomicsDBQueryInfo> queryInfoList = null;
  private long QueryBlockSize = 10000000;
  private long QueryBlockSizeMargin = 500000;

  public GenomicsDBConfiguration() {
    super();
  }
  
  public GenomicsDBConfiguration(Configuration configuration) throws FileNotFoundException {
    super(configuration);
  }

  /**
   * Constructor with partition information
   * @param configuration  Existing configuration object (can contain Hadoop config values)
   * @param list  Contains partition information
   * @throws FileNotFoundException thrown when loader file not found
   */
  public GenomicsDBConfiguration(Configuration configuration, List<GenomicsDBPartitionInfo> list)
    throws FileNotFoundException {
    super(configuration);
    for (GenomicsDBPartitionInfo info : list) {
      addPartitions(info);
    }
  }

  public GenomicsDBConfiguration(Map<String,String> options) throws RuntimeException{
    setOptions(options);
  }

  public void setOptions(Map<String,String> options){
    if(options.containsKey(LOADERJSON)) {
      this.setLoaderJsonFile(options.get(LOADERJSON));
    }
    else {
      throw new RuntimeException("Must specify "+LOADERJSON);
    }
    
    if(options.containsKey(QUERYPB)) {
      this.setQueryPB(options.get(QUERYPB));
    }
    else if(options.containsKey(QUERYJSON)) {
      this.setQueryJsonFile(options.get(QUERYJSON));
    }
    else {
      throw new RuntimeException("Must specify either "+QUERYJSON+" or "+QUERYPB);
    }
    
    if(options.containsKey(MPIHOSTFILE)) {
      try {
        this.setHostFile(options.get(MPIHOSTFILE));
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
    }
  }

  // <String> left for backward compatibility to Java 7
  private ArrayList<String> hosts = new ArrayList<>();

  public GenomicsDBConfiguration setLoaderJsonFile(String path) {
    set(LOADERJSON, path);
    return this;
  }

  @Deprecated
  public GenomicsDBConfiguration setQueryJsonFile(String path) {
    set(QUERYJSON, path);
    return this;
  }

  public GenomicsDBConfiguration setQueryPB(String pb) {
    set(QUERYPB, pb);
    return this;
  }

  public String getLoaderJsonFile() {
    return get(LOADERJSON);
  }

  public String getLoaderPB() {
    return get(LOADERPB);
  }

  public String getQueryJsonFile() {
    return get(QUERYJSON);
  }

  public String getQueryPB() {
    return get(QUERYPB);
  }

  public String getHostFile() {
    return get(MPIHOSTFILE);
  }

  public Boolean hasProtoLoader(){
    return this.get(LOADERPB) != null;
  }

  public Boolean hasProtoQuery(){
    return this.get(QUERYPB) != null;
  }

  /**
   * Host file contains the hosts where GenomicsDB instances reside.
   * This file can be a replica of the slaves file. For now, we have
   * kept it separate, can be merged later
   *
   * @param path  Full path of the host file
   * @return  GenomicsDBConfiguration object
   * @throws FileNotFoundException  If file not found, throw exception
   */
  public GenomicsDBConfiguration setHostFile(String path) throws FileNotFoundException {
    set(MPIHOSTFILE, path);

    FileInputStream pathStream = new FileInputStream(path);
    Scanner scanner = new Scanner(pathStream);
    while (scanner.hasNextLine()) {
      String host = scanner.nextLine();
      hosts.add(host);
    }
    try {
      pathStream.close();
    }
    catch(IOException e) {
      throw new RuntimeException(e);
    }
    return this;
  }

  List<String> getHosts() {
    return hosts;
  }

  private void addPartitions(GenomicsDBPartitionInfo genomicsDBPartitionInfo) {
    if (partitionInfoList==null) {
      partitionInfoList = new ArrayList<>();
    }
    partitionInfoList.add(genomicsDBPartitionInfo);
  }

  /**
   * Return partition list; used when creating input splits.
   * @return Returns ArrayList of PartitionInfo objects
   */
  ArrayList<GenomicsDBPartitionInfo> getPartitions() {
    return partitionInfoList;
  }

  /**
   * Return query range list; used when creating input splits.
   * @return Returns ArrayList of QueryRange objects
   */
  ArrayList<GenomicsDBQueryInfo> getQueryRanges() {
    return queryInfoList;
  }

  /**
   * Return value used to determine optimal query size when creating InputSplits
   * @return Returns QueryBlockSize
   */
  long getQueryBlockSize() {
    return QueryBlockSize;
  }

  /**
   * Return value used to determine "slop" for optimal query size when creating InputSplits
   * @return Returns QueryBlockSizeMargin
   */
  long getQueryBlockSizeMargin() {
    return QueryBlockSizeMargin;
  }

  private void readColumnPartitions(JSONObject obj) throws ParseException {
    if (partitionInfoList==null) {
      partitionInfoList = new ArrayList<>();
    }
    JSONArray colPar = (JSONArray)obj.get("column_partitions");
    if(colPar == null)
      throw new RuntimeException("Could not find attribute \"column_partitions\" in the JSON configuration");
    Iterator it = colPar.iterator();
    while (it.hasNext()) {
      JSONObject obj0 = (JSONObject)it.next();
      String workspace = null, array = null, vcf_output_filename = null;

      long begin = (long)obj0.get("begin");
      workspace = (String)obj0.get("workspace");
      array = (String)obj0.get("array");
      if (array == null) {
        array = (String)obj0.get("array_name");
      }
      vcf_output_filename = (String)obj0.get("vcf_output_filename");
      partitionInfoList.add(new GenomicsDBPartitionInfo(begin, workspace, array, vcf_output_filename));
    }
  }

  @Deprecated
  private void readQueryRanges(JSONObject obj) throws ParseException {
    if (queryInfoList==null) {
      queryInfoList = new ArrayList<>();
    }
    // query_column_ranges is a list of lists; we'll only grab first one
    // Assuming here that query file to Spark interface doesn't have a notion
    // of trying to assign certain queries to certain processes or ranks
    assert obj.containsKey("query_column_ranges");
    JSONArray array = (JSONArray)obj.get("query_column_ranges");
    JSONArray firstList = (JSONArray)array.get(0);
    for (Object currElement : firstList) {
      long start = 0, end = 0;
      if (currElement instanceof JSONArray) {
        JSONArray val = (JSONArray)currElement;
	assert val.size() == 2;
        start = (long)val.get(0);
        end = (long)val.get(1);
      }
      else if (currElement instanceof JSONObject) {
        JSONObject val = (JSONObject)currElement;
	assert val.size() == 1;
        start = end = (long)val.get(0);
      }
      else {
        start = end = (long)currElement;
      }
      queryInfoList.add(new GenomicsDBQueryInfo(start, end));
    }

    if(obj.containsKey("query_block_size")) {
      QueryBlockSize = (long)obj.get("query_block_size");
    }
    if(obj.containsKey("query_block_size_margin")) {
      QueryBlockSizeMargin = (long)obj.get("query_block_size_margin");
    }
  }

  /**
   * parse json file to populate partition and query lists
   * Assuming here that using the Spark interface implies column partitions
   * @param jsonType json file to use while loading - either LOADERJSON or QUERYJSON
   * @throws FileNotFoundException  Thrown if queryJson file isn't found
   * @throws IOException  Thrown if other IO exception while handling file operations
   * @throws ParseException  Thrown if JSON parsing fails
   */
  void populateListFromJson(String jsonType) 
		  throws FileNotFoundException, IOException, ParseException {
    JSONParser parser = new JSONParser();
    String filename = this.get(jsonType, "");
    if (filename.isEmpty()) {
      throw new IOException(String.format("No filename specified with type=%s in GenomicdDBConfiguration", jsonType));
    }
    if (!new File(filename).exists()) {
      throw new IOException(String.format("Could not find file=%s associated with type=%s", filename, jsonType));
    }
    FileReader jsonReader = new FileReader(get(jsonType));
    try {
      JSONObject obj = (JSONObject)parser.parse(jsonReader);

      if (jsonType.equals(LOADERJSON)) {
        readColumnPartitions(obj);
      }
      else if (jsonType.equals(QUERYJSON)) {
        readQueryRanges(obj);
      }
    }
    finally {
      jsonReader.close();
    }
  }

  private void readColumnPartitionsPB(String pb) throws
        com.google.protobuf.InvalidProtocolBufferException {
    GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder = 
             GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();
    byte[] pbDecoded = Base64.getDecoder().decode(pb);
    importConfigurationBuilder.mergeFrom(pbDecoded);
    GenomicsDBImportConfiguration.ImportConfiguration loaderPB = importConfigurationBuilder.build();

    if (partitionInfoList==null) {
      partitionInfoList = new ArrayList<>();
    }

    for (GenomicsDBImportConfiguration.Partition partition : 
            loaderPB.getColumnPartitionsList()) {
      if (!partition.getBegin().hasTiledbColumn()) {
        throw new RuntimeException("Spark layer doesn't support specifying Partitions as contig position");
      }
      partitionInfoList.add(new GenomicsDBPartitionInfo(
              partition.getBegin().getTiledbColumn(),
              partition.getWorkspace(),
              partition.getArrayName(),
              partition.getVcfOutputFilename()));
    }
  }

  private void readQueryRangesPB(String pb) throws 
        com.google.protobuf.InvalidProtocolBufferException {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder exportConfigurationBuilder = 
             GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
    byte[] pbDecoded = Base64.getDecoder().decode(pb);
    exportConfigurationBuilder.mergeFrom(pbDecoded);
    GenomicsDBExportConfiguration.ExportConfiguration queryPB = exportConfigurationBuilder.build();

    if (queryInfoList==null) {
      queryInfoList = new ArrayList<>();
    }

    for (GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList 
            qcrList : queryPB.getQueryColumnRangesList()) {
      for (Coordinates.GenomicsDBColumnOrInterval colOrInterval : qcrList.getColumnOrIntervalListList()) {
        long start = -1;
        long end = -1;
        if (colOrInterval.hasColumn()) {
          if(colOrInterval.getColumn().hasTiledbColumn()) {
            start = end = colOrInterval.getColumn().getTiledbColumn();
          }
        }
        else if (colOrInterval.hasColumnInterval()) {
          if(colOrInterval.getColumnInterval().hasTiledbColumnInterval()) {
            start = colOrInterval.getColumnInterval().getTiledbColumnInterval().getBegin();
            end = colOrInterval.getColumnInterval().getTiledbColumnInterval().getEnd();
          }
        }
        if (start == -1 || end == -1) {
          throw new RuntimeException("Query range must be specified as tiledb columns or intervals");
        }
        queryInfoList.add(new GenomicsDBQueryInfo(start, end));
      }
    }
    if(queryPB.hasSparkConfig()) {
      if(queryPB.getSparkConfig().hasQueryBlockSize()) {
        QueryBlockSize = queryPB.getSparkConfig().getQueryBlockSize();
      }
      if(queryPB.getSparkConfig().hasQueryBlockSizeMargin()) {
        QueryBlockSizeMargin = queryPB.getSparkConfig().getQueryBlockSizeMargin();
      }
    }
  }

  /**
   * parse protobuf to populate partition and query lists
   * Assuming here that using the Spark interface implies column partitions
   * @param pbType protobuf type to use - either LOADERPB or QUERYPB
   */
  void populateListFromPB(String pbType) {
    try {
      if (pbType.equals(LOADERPB)) {
        readColumnPartitionsPB(get(pbType));
      }
      else if (pbType.equals(QUERYPB)) {
        readQueryRangesPB(get(pbType));
      }
    }
    catch (Exception e) {
        System.err.println("Could not parse protobuf while populating partition and query lists");
        e.printStackTrace();
    }
  }
}

