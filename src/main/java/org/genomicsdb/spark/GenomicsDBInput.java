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

import org.apache.spark.sql.types.StructType;

import org.json.simple.JSONObject;
import org.json.simple.JSONArray;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.util.Map;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.lang.RuntimeException;
import java.lang.InstantiationException;
import java.lang.IllegalAccessException;

/**
 * The input class represents all the data being queried from GenomicsDB.
 * This can be used both for hadoop map-reduce style InputSplits and
 * for Spark (datasourcev2) InputPartitions
 */
public class GenomicsDBInput<T extends GenomicsDBInputInterface> {

  public GenomicsDBConfiguration genomicsDBConfiguration;
  private StructType schema;
  private Map<String, GenomicsDBVidSchema> vMap;
  private long minQueryBlockSize;
  private long maxQueryBlockSize;
  private Class<T> clazz;

  /**
    * constructor
    */
  GenomicsDBInput(GenomicsDBConfiguration gdbconf, StructType schema, 
      Map<String, GenomicsDBVidSchema> vMap, long minQBS, long maxQBS, Class<T> clazz) {
    genomicsDBConfiguration = gdbconf;
    this.clazz = clazz;
    this.schema = schema;
    this.vMap = vMap;
    minQueryBlockSize = minQBS;
    maxQueryBlockSize = maxQBS;
  }

  /**
   * set the GenomicsDBConfiguration object
   */
  public void setGenomicsDBConfiguration(GenomicsDBConfiguration gdbconf) {
    genomicsDBConfiguration = gdbconf;
  }

  /**
   * get the GenomicsDBConfiguration object
   */
  public GenomicsDBConfiguration getGenomicsDBConfiguration() {
    return genomicsDBConfiguration;
  }

  /**
   * create an object of the correct type based on clazz:
   * InputPartition or InputSplit
   * @return Returns the appropriate partition/split object
   */
  private T createInputInstance() {
    try {
      return clazz.newInstance();
    }
    catch (IllegalAccessException | InstantiationException e) {
      e.printStackTrace();
    }
    return null;
  }

  /**
   * create the appropriate partition/split object and initialize it
   * @param host host location for InputSplit
   * @return Returns initialized InputSplit or InputPartition
   * @throws RuntimeException Thrown if clazz is not GenomicsDBInputPartition or
   * GenomicsDBInputSplit
   */
  private T getInputInstance(String host) throws RuntimeException {
    T instance = createInputInstance();
    // seems sorta hacky to instantiate the splits/partitions this way
    // but should work and we shouldn't have to
    // scale this up to many more different classes...
    if (GenomicsDBInputPartition.class.isAssignableFrom(clazz)) {
      instance.setGenomicsDBConf(genomicsDBConfiguration);
      instance.setGenomicsDBSchema(schema);
      instance.setGenomicsDBVidSchema(vMap);
      return instance;
    }
    else if (GenomicsDBInputSplit.class.isAssignableFrom(clazz)) {
      instance.setHost(host);
      return instance;
    }
    else {
      throw new RuntimeException("Unsupported class for GenomicsDBInput:"+clazz.getName());
    }
  }

  /**
   * create the appropriate partition/split object and initialize it
   * @param part GenomicsDBPartitionInfo for given partition/split
   * @param qrange GenomicsDBQueryInfo for given partition/split
   * @return Returns initialized InputSplit or InputPartition
   * @throws RuntimeException Thrown if clazz is not GenomicsDBInputPartition or
   * GenomicsDBInputSplit
   */
  private T getInputInstance(GenomicsDBPartitionInfo part, GenomicsDBQueryInfo qrange) throws RuntimeException {
    T instance = createInputInstance();
    // seems sorta hacky to instantiate the splits/partitions this way
    // but should work and we shouldn't have to
    // scale this up to many more different classes...
    if (GenomicsDBInputPartition.class.isAssignableFrom(clazz)) {
      instance.setPartitionInfo(part);
      instance.setQueryInfo(qrange);
      instance.setGenomicsDBConf(genomicsDBConfiguration);
      instance.setGenomicsDBSchema(schema);
      instance.setGenomicsDBVidSchema(vMap);
      return instance;
    }
    else if (GenomicsDBInputSplit.class.isAssignableFrom(clazz)) {
      instance.setPartitionInfo(part);
      instance.setQueryInfo(qrange);
      return instance;
    }
    else {
      throw new RuntimeException("Unsupported class for GenomicsDBInput:"+clazz.getName());
    }
  }

  /**
    * Get the schema for DataSource.
    * Only relevant for Spark's Datasourcev2 implementation
    * @return Returns schema associated with DataSourceReaderV2
    */
  public StructType getSchema() {
    return schema;
  }

  /**
    * Divide input data/datasource into chunks so that
    * we can distribute the work amongst many workers.
    * Called by InputFormat::getSplits and 
    * DataSourceReader::planInputPartitions
    * @return Returns list of "work chunks" (either InputSplits or
    * InputPartitions)
    */
  public List<T> divideInput() {

    try {
      genomicsDBConfiguration.populateListFromJson(GenomicsDBConfiguration.LOADERJSON);
      genomicsDBConfiguration.populateListFromJson(GenomicsDBConfiguration.QUERYJSON);
    }
    catch (IOException | ParseException e) {
      e.printStackTrace();
      return null;
    }

    ArrayList<GenomicsDBPartitionInfo> partitionsList = genomicsDBConfiguration.getPartitions();
    ArrayList<GenomicsDBQueryInfo> queryRangeList = genomicsDBConfiguration.getQueryRanges();

    long goalBlockSize = Math.max(minQueryBlockSize, 
		           Math.min(genomicsDBConfiguration.getQueryBlockSize(), maxQueryBlockSize));

    ArrayList<T> inputPartitions = new ArrayList<T>();
    // For now, assuming that each of partitions and queryRange are sorted
    // by start position, and that query ranges don't overlap.
    // TODO: not sure anything in GenomicsDB enforces the above assumption....
    int pIndex = 0;
    int qIndex = 0;
    GenomicsDBPartitionInfo partition = null;
    if (!partitionsList.isEmpty())
      partition = partitionsList.get(pIndex);

    if (partition != null) {
      while (qIndex < queryRangeList.size() && partition != null) {
        GenomicsDBQueryInfo queryRange = queryRangeList.get(qIndex);
  
        // advance partition index if needed
        // i.e., advance till we find the partition that contains the query begin position
        while ((pIndex + 1) < partitionsList.size() && partition.getBeginPosition() < queryRange.getBeginPosition()) {
          pIndex++;
  	  partition = partitionsList.get(pIndex);
        }
        if (partition.getBeginPosition() > queryRange.getBeginPosition()) {
          pIndex--;
  	  partition = partitionsList.get(pIndex);
        }
  
        long queryBlockSize = queryRange.getEndPosition() - queryRange.getBeginPosition() + 1;
        if (queryBlockSize < goalBlockSize) {
          inputPartitions.add(getInputInstance(partition, queryRange));
        }
        else {
          // bigger than goalBlockSize, so break up into "query chunks"
  
  	  long queryBlockStart = queryRange.getBeginPosition();
	  long queryBlockMargin = genomicsDBConfiguration.getQueryBlockSizeMargin();
  	  while (queryBlockStart < queryRange.getEndPosition()) {
            long blockSize = (queryBlockSize > (goalBlockSize+queryBlockMargin)) ? goalBlockSize : queryBlockSize;
  	    GenomicsDBQueryInfo queryBlock = new GenomicsDBQueryInfo(queryBlockStart, queryBlockStart + blockSize - 1);
  	    inputPartitions.add(getInputInstance(partition, queryBlock));
  
  	    // if this queryBlock spans multiple partitions, need to add those as splits as well
  	    while ((pIndex + 1) < partitionsList.size() &&
                    queryBlockStart + blockSize - 1 >= partitionsList.get(pIndex+1).getBeginPosition()) {
  	      pIndex++;
              partition = partitionsList.get(pIndex);
  	      inputPartitions.add(getInputInstance(partition, queryBlock));
  	    }
  	    queryBlockStart += blockSize;
  	    queryBlockSize -= blockSize;
  	  }
        }
        qIndex++;
      }
    }
    return inputPartitions;
  }

  /**
   * Creates tmp query file based on partition, query and existing query file
   *
   * @param queryJson Existing query json file
   * @param partition used to populate array
   * @param query used to bound query column ranges
   * @return  Returns path to temporary query file
   * @throws FileNotFoundException  Thrown if queryJson file isn't found
   * @throws IOException  Thrown if other IO exception while handling file operations
   * @throws ParseException  Thrown if JSON parsing fails
   */
  public static String createTmpQueryFile(String queryJson, GenomicsDBPartitionInfo partition, GenomicsDBQueryInfo query) 
		  throws FileNotFoundException, IOException, ParseException {
    JSONObject amended = new JSONObject();
    amended.put("array",partition.getArrayName());
    if (query.getBeginPosition() == query.getEndPosition()) {
      JSONArray l1 = new JSONArray();
      ArrayList<Long> a = new ArrayList<>(1);
      a.add(query.getBeginPosition());
      l1.add(a);
      amended.put("query_column_ranges",l1);
    }
    else {
      JSONArray l1 = new JSONArray();
      ArrayList<ArrayList<Long>> a = new ArrayList<>(1);
      a.add(new ArrayList<Long>(2));
      a.get(0).add(query.getBeginPosition());
      a.get(0).add(query.getEndPosition());
      l1.add(a);
      amended.put("query_column_ranges", l1);
    }

    try {

      JSONParser parser = new JSONParser();
      FileReader queryJsonReader = new FileReader(queryJson);
      JSONObject obj = null;
      try {
        obj = (JSONObject)parser.parse(queryJsonReader);
      }
      catch(ParseException | IOException e) {
        queryJsonReader.close();
        throw e;
      }
  
      obj.forEach((k, v) -> {
        if (!(k.equals("query_column_ranges") || k.equals("array"))) {
          amended.put(k, v);
        }
      });
      queryJsonReader.close();
    }
    catch (FileNotFoundException e) {
      e.printStackTrace();
      return null;
    }
    File tmpQueryFile = File.createTempFile("queryJson", ".json");
    tmpQueryFile.deleteOnExit();
    FileWriter fptr = new FileWriter(tmpQueryFile);
    try {
        fptr.write(amended.toJSONString());
    }
    catch(IOException e) {
        fptr.close();
        throw new IOException(e);
    }
    fptr.close();
    return tmpQueryFile.getAbsolutePath();
  }

}
