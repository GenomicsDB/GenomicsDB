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
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.spark.sources.GenomicsDBInputPartition;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.*;

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
  public GenomicsDBInput(GenomicsDBConfiguration gdbconf, StructType schema,
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
      return clazz.getConstructor().newInstance();
    }
    catch (NoSuchMethodException | InvocationTargetException | InstantiationException | IllegalAccessException e) {
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
  }

  /**
   * create the appropriate partition/split object and initialize it
   * @param part GenomicsDBPartitionInfo for given partition/split
   * @param qrangeList GenomicsDBQueryInfo for given partition/split
   * @return Returns initialized InputSplit or InputPartition
   * @throws RuntimeException Thrown if clazz is not GenomicsDBInputPartition or
   * GenomicsDBInputSplit
   */
  private T getInputInstance(GenomicsDBPartitionInfo part, 
      ArrayList<GenomicsDBQueryInfo> qrangeList) throws RuntimeException {
    T instance = createInputInstance();
    // seems sorta hacky to instantiate the splits/partitions this way
    // but should work and we shouldn't have to
    // scale this up to many more different classes...
    if (GenomicsDBInputSplit.class.isAssignableFrom(clazz)) {
      instance.setPartitionInfo(part);
      instance.setQueryInfoList(qrangeList);
      return instance;
    }
    else {
      if (GenomicsDBInputPartition.class.isAssignableFrom(clazz)) {
        instance.setPartitionInfo(part);
        instance.setQueryInfoList(qrangeList);
        instance.setGenomicsDBConf(genomicsDBConfiguration);
        instance.setGenomicsDBSchema(schema);
        instance.setGenomicsDBVidSchema(vMap);
        return instance;
      }
      else {
        throw new RuntimeException("Unsupported class for GenomicsDBInput:"+clazz.getName());
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
  @SuppressWarnings("deprecation")
  public List<T> divideInput() {

    if (genomicsDBConfiguration.hasProtoLoader()){
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

    if (genomicsDBConfiguration.hasProtoQuery()) {
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
      // create a temporary arraylist that we'll use if we glom queries
      ArrayList<GenomicsDBQueryInfo> glomQuerys = new ArrayList<GenomicsDBQueryInfo>();
      // we'll use glommedSize to separate inputsplits if accumulated
      // size of glommed queries goes above goalBlockSize
      long glommedSize = 0;
      while (qIndex < queryRangeList.size() && partition != null) {
        GenomicsDBQueryInfo queryRange = queryRangeList.get(qIndex);
  
        // advance partition index if needed
        // i.e., advance till we find the partition that contains the query begin position
        while ((pIndex + 1) < partitionsList.size() && 
              partitionsList.get(pIndex+1).getBeginPosition() <= 
              queryRange.getBeginPosition()) {
          // add glommed queries to inputsplit using previous partition since
          // we're moving on to new partitions
          if (!glomQuerys.isEmpty()) {
            inputPartitions.add(getInputInstance(partition, glomQuerys));
            glomQuerys.clear();
          }
          glommedSize = 0;
          pIndex++;
  	  partition = partitionsList.get(pIndex);
        }
  
        long queryBlockSize = queryRange.getEndPosition() - queryRange.getBeginPosition() + 1;
        if (queryBlockSize < goalBlockSize) {
          // create glommed inputsplit
          if(doesQuerySpanPartitions(pIndex+1, partitionsList, queryRange.getEndPosition())) {
            // if current query spans multiple partitions then add to previous inputsplit as well
            glomQuerys.add(queryRange);
            inputPartitions.add(getInputInstance(partition, glomQuerys));
            glomQuerys.clear();
            glommedSize = 0;
  	    // if this queryBlock spans multiple partitions, need to add those as splits as well
  	    // can use the same ArrayList of queries since each inputsplit will only care
  	    // about the section that is relevant to its partition
            glomQuerys.add(queryRange);
            pIndex = addSplitsIfQuerySpansPartitions(inputPartitions, pIndex, 
                    queryRange.getEndPosition(), partitionsList, glomQuerys);
            partition = (pIndex < partitionsList.size()) ? partitionsList.get(pIndex) : null;
            glomQuerys.clear();
          }
          else {
            // add query to glom
            glomQuerys.add(queryRange);
            glommedSize += queryBlockSize;
            // add to inputsplit and reset glom if accumulated glom size hits goal
            if (glommedSize >= goalBlockSize) {
    	      inputPartitions.add(getInputInstance(partition, glomQuerys));
              glomQuerys.clear();
              glommedSize = 0;
            }
          }
        }
        else {
          if (!glomQuerys.isEmpty()) {
            inputPartitions.add(getInputInstance(partition, glomQuerys));
            glomQuerys.clear();
            glommedSize = 0;
          }
          // bigger than goalBlockSize, so break up into "query chunks"
  
  	  long queryBlockStart = queryRange.getBeginPosition();
	  long queryBlockMargin = genomicsDBConfiguration.getQueryBlockSizeMargin();
  	  while (queryBlockStart < queryRange.getEndPosition()) {
            long blockSize = (queryBlockSize > (goalBlockSize+queryBlockMargin)) ? goalBlockSize : queryBlockSize;
  	    GenomicsDBQueryInfo queryBlock = new GenomicsDBQueryInfo(queryBlockStart, queryBlockStart + blockSize - 1);
            glomQuerys.add(queryBlock);
  	    inputPartitions.add(getInputInstance(partition, glomQuerys));
  
  	    // if this queryBlock spans multiple partitions, need to add those as splits as well
  	    pIndex = addSplitsIfQuerySpansPartitions(inputPartitions, pIndex,
                    queryBlockStart+blockSize-1, partitionsList, glomQuerys);
            partition = (pIndex < partitionsList.size()) ? partitionsList.get(pIndex) : null;
            glomQuerys.clear();
  	    queryBlockStart += blockSize;
  	    queryBlockSize -= blockSize;
  	  }
        }
        qIndex++;
      }
      // if we still have glommed queries that haven't been assigned to
      // an inputsplit, they go with the final partition
      if(!glomQuerys.isEmpty()) {
        inputPartitions.add(getInputInstance(partition, glomQuerys));
      }
    }
    return inputPartitions;
  }

  private boolean doesQuerySpanPartitions(final int index, 
          final ArrayList<GenomicsDBPartitionInfo> list,
          final long queryEnd) {
    return index < list.size() &&
            queryEnd >= list.get(index).getBeginPosition();
  }

  private int addSplitsIfQuerySpansPartitions(
          ArrayList<T> list, int index,
          final long queryEnd,
          final ArrayList<GenomicsDBPartitionInfo> partitionList, 
          final ArrayList<GenomicsDBQueryInfo> queryList) {
    GenomicsDBPartitionInfo partition = partitionList.get(index);
    while(doesQuerySpanPartitions(index+1, partitionList, queryEnd)) {
      index++;
      partition = partitionList.get(index);
      list.add(getInputInstance(partition, queryList));
    }
    return index;
  }

  /**
   * Creates export configuration protobuf object 
   * based on partition, query and existing query file or protobuf
   *
   * @param queryFileOrPB Existing query json file or base64 encoded protobuf byte data
   * @param partition used to populate array
   * @param queryList used to bound query column ranges
   * @param isPB boolean parameter that denotes if queryFileOrPB is protobuf
   * @return  Returns export configuration protobuf object
   * @throws IOException  Thrown if other IO exception while handling file operations
   * @throws ParseException  Thrown if JSON parsing fails
   */
  public static GenomicsDBExportConfiguration.ExportConfiguration 
        createTargetExportConfigurationPB(String queryFileOrPB, 
        GenomicsDBPartitionInfo partition, 
        ArrayList<GenomicsDBQueryInfo> queryList, 
        boolean isPB) 
        throws IOException, ParseException {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder 
        exportConfigurationBuilder;

    if (isPB) {
      exportConfigurationBuilder = 
          GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
      byte[] queryPB = Base64.getDecoder().decode(queryFileOrPB);
      exportConfigurationBuilder.mergeFrom(queryPB);
    }
    else {
      exportConfigurationBuilder =
          getExportConfigurationFromJsonFile(queryFileOrPB);
    }
    // set the array name per the targeted partition
    exportConfigurationBuilder.setArrayName(partition.getArrayName());
    exportConfigurationBuilder.clearQueryColumnRanges();

    GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.Builder 
        queryRangeBuilder = 
        GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.newBuilder();
    for (GenomicsDBQueryInfo query : queryList) {
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
    }

    return exportConfigurationBuilder
            .addQueryColumnRanges(queryRangeBuilder).build();
  }

  @Deprecated
  public static GenomicsDBExportConfiguration.ExportConfiguration.Builder 
        getExportConfigurationFromJsonFile(String queryFile)
        throws IOException, ParseException {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder 
        exportConfigurationBuilder = 
        GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();

    try {
      JSONParser parser = new JSONParser();
      FileReader queryJsonReader = new FileReader(queryFile);
      HashMap<?,?> obj = null;
      try {
        obj = (HashMap<?,?>)parser.parse(queryJsonReader);
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
          case "scan_full":
            exportConfigurationBuilder.setScanFull(
                val.toString().equals("true")); 
            break;
          case "sites_only_query":
            exportConfigurationBuilder.setSitesOnlyQuery(
                val.toString().equals("true")); 
            break;
          case "enable_shared_posixfs_optimizations":
            exportConfigurationBuilder.setEnableSharedPosixfsOptimizations(
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
          case "vcf_header_filename":
            if (val instanceof JSONArray) {
              exportConfigurationBuilder.setVcfHeaderFilename(((JSONArray)val).get(0).toString());
            }
            else {
              exportConfigurationBuilder.setVcfHeaderFilename(val.toString());
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
          case "query_sample_names":
            for (Object element : (JSONArray)val) {
              exportConfigurationBuilder.addQuerySampleNames(element.toString());
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
          case "segment_size":
            exportConfigurationBuilder.setSegmentSize(((Long)val).intValue());
            break;
          case "query_block_size":
            {
              if(exportConfigurationBuilder.hasSparkConfig()) {
                exportConfigurationBuilder.getSparkConfigBuilder().setQueryBlockSize((long)((Object)val));
              }
              else {
                GenomicsDBExportConfiguration.SparkConfig.Builder configBuilder = 
                        GenomicsDBExportConfiguration.SparkConfig.newBuilder();
                exportConfigurationBuilder.setSparkConfig(configBuilder.setQueryBlockSize((long)((Object)val)).build());
              }
            }
            break;
          case "query_block_size_margin":
            {
              if(exportConfigurationBuilder.hasSparkConfig()) {
                exportConfigurationBuilder.getSparkConfigBuilder().setQueryBlockSizeMargin((long)((Object)val));
              }
              else {
                GenomicsDBExportConfiguration.SparkConfig.Builder configBuilder = 
                        GenomicsDBExportConfiguration.SparkConfig.newBuilder();
                exportConfigurationBuilder.setSparkConfig(configBuilder.setQueryBlockSizeMargin((long)((Object)val)).build());
              }
            }
            break;
          case "callset_mapping_file":
            exportConfigurationBuilder.setCallsetMappingFile(val.toString()); break;
          case "vid_mapping_file":
            exportConfigurationBuilder.setVidMappingFile(val.toString()); break;
          default:
            System.err.println("Ignoring attribute:"+key.toString()+
                    " with val:"+val.toString());
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
