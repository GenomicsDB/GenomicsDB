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

import org.genomicsdb.model.*;

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
import java.lang.ClassNotFoundException;

import com.googlecode.protobuf.format.JsonFormat;

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
    * constructor for GenomicsDBInput
    * @param gdbconf GenomicsDBConfiguration object
    * @param schema schema used for the datasource API
    * @param vMap map of attribute to vid mapping information
    * @param minQBS minimum query block size used for partitioning query
    * @param maxQBS maximum query block size used for partitioning query
    * @param clazz Class object used to decide how to instantiate partitions
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
   * @param gdbconf GenomcisDBConfiguration object
   */
  public void setGenomicsDBConfiguration(GenomicsDBConfiguration gdbconf) {
    genomicsDBConfiguration = gdbconf;
  }

  /**
   * get the GenomicsDBConfiguration object
   * @return GenomicsDBConfiguration object
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
    if (GenomicsDBInputSplit.class.isAssignableFrom(clazz)) {
      instance.setHost(host);
      return instance;
    }
    else {
      try {
        Class c = Class.forName("org.genomicsdb.spark.GenomicsDBInputPartition");
        if (GenomicsDBInputPartition.class.isAssignableFrom(clazz)) {
          instance.setGenomicsDBConf(genomicsDBConfiguration);
          instance.setGenomicsDBSchema(schema);
          instance.setGenomicsDBVidSchema(vMap);
          return instance;
        }
        else {
          throw new RuntimeException("Unsupported class for GenomicsDBInput:"+clazz.getName());
        }
      }
      catch (ClassNotFoundException ex) {
        throw new RuntimeException("Warning: Could not find GenomicsDBInputPartition. " +
                                   "Datasourceapi only works with >= Spark 2.4.0");
      }
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
    if (GenomicsDBInputSplit.class.isAssignableFrom(clazz)) {
      instance.setPartitionInfo(part);
      instance.setQueryInfo(qrange);
      return instance;
    }
    else {
      try {
        Class c = Class.forName("org.genomicsdb.spark.GenomicsDBInputPartition");
        if (GenomicsDBInputPartition.class.isAssignableFrom(clazz)) {
          instance.setPartitionInfo(part);
          instance.setQueryInfo(qrange);
          instance.setGenomicsDBConf(genomicsDBConfiguration);
          instance.setGenomicsDBSchema(schema);
          instance.setGenomicsDBVidSchema(vMap);
          return instance;
        }
        else {
          throw new RuntimeException("Unsupported class for GenomicsDBInput:"+clazz.getName());
        }
      }
      catch (ClassNotFoundException ex) {
        throw new RuntimeException("Warning: Could not find GenomicsDBInputPartition. " +
                                   "Datasourceapi only works with >= Spark 2.4.0");
      }
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

    if (genomicsDBConfiguration.get(GenomicsDBConfiguration.LOADERPB) != null) {
      genomicsDBConfiguration.populateListFromPB(GenomicsDBConfiguration.LOADERPB);
    } else {
      try {
        genomicsDBConfiguration.populateListFromJson(GenomicsDBConfiguration.LOADERJSON);
      }
      catch (IOException | ParseException e) {
        e.printStackTrace();
        return null;
      }
    }

    if (genomicsDBConfiguration.get(GenomicsDBConfiguration.QUERYPB) != null) {
      genomicsDBConfiguration.populateListFromPB(GenomicsDBConfiguration.QUERYPB);
    } else {
      try {
        genomicsDBConfiguration.populateListFromJson(GenomicsDBConfiguration.QUERYJSON);
      }
      catch (IOException | ParseException e) {
        e.printStackTrace();
        return null;
      }
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
  	  // if this queryRange spans multiple partitions, need to add those as splits as well
  	  while ((pIndex + 1) < partitionsList.size() &&
                  queryRange.getEndPosition() >= partitionsList.get(pIndex+1).getBeginPosition()) {
  	    pIndex++;
            partition = partitionsList.get(pIndex);
  	    inputPartitions.add(getInputInstance(partition, queryRange));
  	  }
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
   * Creates export configuration protobuf object 
   * based on partition, query and existing query file or protobuf
   *
   * @param queryFileOrPB Existing query json file or protobuf string
   * @param partition used to populate array
   * @param query used to bound query column ranges
   * @param isPB boolean parameter that denotes if queryFileOrPB is protobug string
   * @return  Returns path to temporary query file
   * @throws FileNotFoundException  Thrown if queryJson file isn't found
   * @throws IOException  Thrown if other IO exception while handling file operations
   * @throws ParseException  Thrown if JSON parsing fails
   */
  public static GenomicsDBExportConfiguration.ExportConfiguration 
        createTargetExportConfigurationPB(String queryFileOrPB, 
        GenomicsDBPartitionInfo partition, GenomicsDBQueryInfo query, 
        boolean isPB) 
        throws FileNotFoundException, IOException, ParseException {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder 
        exportConfigurationBuilder;

    if (isPB) {
      exportConfigurationBuilder = 
          GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
      JsonFormat.merge(queryFileOrPB, exportConfigurationBuilder);
    }
    else {
      exportConfigurationBuilder = 
          getExportConfigurationFromJsonFile(queryFileOrPB);
    }
    // set the array name per the targeted partition
    exportConfigurationBuilder.setArrayName(partition.getArrayName());
    exportConfigurationBuilder.clearQueryColumnRanges();
    // set query range per the part of the query this will handle
    GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.Builder 
        queryRangeBuilder = 
        GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.newBuilder();
    if (query.getBeginPosition() == query.getEndPosition()) {
      Coordinates.GenomicsDBColumnOrInterval.Builder columnBuilder = 
              Coordinates.GenomicsDBColumnOrInterval.newBuilder();
      Coordinates.GenomicsDBColumn.Builder gdbColumnBuilder = 
              Coordinates.GenomicsDBColumn.newBuilder();

      gdbColumnBuilder.setTiledbColumn(query.getBeginPosition());
      queryRangeBuilder.addColumnOrIntervalList(
          columnBuilder.setColumn(gdbColumnBuilder));
    }
    else {
      Coordinates.GenomicsDBColumnOrInterval.Builder intervalBuilder = 
              Coordinates.GenomicsDBColumnOrInterval.newBuilder();
      Coordinates.GenomicsDBColumnInterval.Builder gdbColumnIntervalBuilder = 
              Coordinates.GenomicsDBColumnInterval.newBuilder();
      Coordinates.TileDBColumnInterval.Builder tdbColumnIntervalBuilder = 
              Coordinates.TileDBColumnInterval.newBuilder();

      tdbColumnIntervalBuilder.setBegin(query.getBeginPosition())
              .setEnd(query.getEndPosition());
      queryRangeBuilder.addColumnOrIntervalList(
              intervalBuilder.setColumnInterval(
              gdbColumnIntervalBuilder.setTiledbColumnInterval(
              tdbColumnIntervalBuilder)));
    }

    return exportConfigurationBuilder
            .addQueryColumnRanges(queryRangeBuilder).build();
  }

  public static GenomicsDBExportConfiguration.ExportConfiguration.Builder 
        getExportConfigurationFromJsonFile(String queryFile) 
        throws FileNotFoundException, IOException, ParseException {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder 
        exportConfigurationBuilder = 
        GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();

    try {
      JSONParser parser = new JSONParser();
      FileReader queryJsonReader = new FileReader(queryFile);
      JSONObject obj = null;
      try {
        obj = (JSONObject)parser.parse(queryJsonReader);
      }
      catch(ParseException | IOException e) {
        queryJsonReader.close();
        throw e;
      }
  
      obj.forEach((key, val) -> {
        switch (key.toString()) {
          case "workspace":
            exportConfigurationBuilder.setWorkspace(val.toString()); break;
          case "reference_genome":
            exportConfigurationBuilder.setReferenceGenome(val.toString()); break;
          case "produce_GT_field":
            exportConfigurationBuilder.setProduceGTField(
                val.toString().equals("true")); 
            break;
          case "produce_FILTER_field":
            exportConfigurationBuilder.setProduceFILTERField(
                val.toString().equals("true")); 
            break;
          case "query_filter":
            exportConfigurationBuilder.setQueryFilter(val.toString()); break;
          case "query_attributes":
            for (Object element : (JSONArray)val) {
              exportConfigurationBuilder.addAttributes(element.toString());
            }
            break;
          case "query_row_ranges":
            JSONArray list = (JSONArray)((JSONArray)val).get(0);
            for (Object element : list) {
              JSONArray value = (JSONArray)element;
              GenomicsDBExportConfiguration.RowRangeList.Builder queryRowRangesBuilder = 
                      GenomicsDBExportConfiguration.RowRangeList.newBuilder();
              GenomicsDBExportConfiguration.RowRange.Builder rowRangeBuilder = 
                      GenomicsDBExportConfiguration.RowRange.newBuilder();
              rowRangeBuilder.setLow((long)value.get(0))
                      .setHigh((long)value.get(1));
    
              exportConfigurationBuilder.addQueryRowRanges(queryRowRangesBuilder
                      .addRangeList(rowRangeBuilder));
            }
            break;
          case "array":
          case "array_name":
            exportConfigurationBuilder.setArrayName(val.toString()); break;
          case "query_column_ranges":
            JSONArray columnList = (JSONArray)((JSONArray)val).get(0);
            GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.Builder 
                queryRangeBuilder = 
                GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.newBuilder();
            for (Object element : columnList) {
              long start = 0, end = 0;
              if (element instanceof JSONArray) {
                JSONArray v = (JSONArray)element;
                assert v.size() == 2;
                start = (long)v.get(0);
                end = (long)v.get(1);
              }
              else if (element instanceof JSONObject) {
                JSONObject v = (JSONObject)element;
                assert v.size() == 1;
                start = end = (long)v.get(0);
              }
              else {
                start = end = (long)element;
              }
              if (start == end) {
                Coordinates.GenomicsDBColumnOrInterval.Builder columnBuilder = 
                        Coordinates.GenomicsDBColumnOrInterval.newBuilder();
                Coordinates.GenomicsDBColumn.Builder gdbColumnBuilder = 
                        Coordinates.GenomicsDBColumn.newBuilder();
          
                gdbColumnBuilder.setTiledbColumn(start);
                queryRangeBuilder.addColumnOrIntervalList(
                    columnBuilder.setColumn(gdbColumnBuilder));
              }
              else {
                Coordinates.GenomicsDBColumnOrInterval.Builder intervalBuilder = 
                        Coordinates.GenomicsDBColumnOrInterval.newBuilder();
                Coordinates.GenomicsDBColumnInterval.Builder gdbColumnIntervalBuilder = 
                        Coordinates.GenomicsDBColumnInterval.newBuilder();
                Coordinates.TileDBColumnInterval.Builder tdbColumnIntervalBuilder = 
                        Coordinates.TileDBColumnInterval.newBuilder();
          
                tdbColumnIntervalBuilder.setBegin(start)
                        .setEnd(end);
                queryRangeBuilder.addColumnOrIntervalList(
                        intervalBuilder.setColumnInterval(
                        gdbColumnIntervalBuilder.setTiledbColumnInterval(
                        tdbColumnIntervalBuilder)));
              }
            }
            exportConfigurationBuilder
                    .addQueryColumnRanges(queryRangeBuilder);
            break;

          case "query_block_size":
            JSONObject qbs = (JSONObject)val;
            exportConfigurationBuilder.setQueryBlockSize((long)qbs.get(0));
            break;
          case "query_block_size_margin":
            JSONObject qbsm = (JSONObject)val;
            exportConfigurationBuilder.setQueryBlockSizeMargin((long)qbsm.get(0));
            break;
          default:
            // we'll ignore all other fields
        }
      });
      queryJsonReader.close();
    }
    catch (FileNotFoundException e) {
      e.printStackTrace();
      return null;
    }
    return exportConfigurationBuilder;
  }
}
