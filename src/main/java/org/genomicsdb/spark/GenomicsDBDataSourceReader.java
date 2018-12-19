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

import org.apache.spark.sql.sources.v2.reader.DataSourceReader;
import org.apache.spark.sql.sources.v2.DataSourceOptions;
import org.apache.spark.sql.sources.v2.reader.InputPartition;
import org.apache.spark.sql.catalyst.InternalRow;
import org.apache.spark.sql.types.*;
import org.apache.spark.sql.types.DataTypes.*;
import org.apache.spark.sql.*;
import org.apache.spark.ml.linalg.SQLDataTypes;

import org.json.simple.parser.ParseException;

import scala.collection.JavaConverters;
import java.util.List;
import java.util.ArrayList;
import java.io.IOException;
import java.util.Optional;
import java.io.FileNotFoundException;

public class GenomicsDBDataSourceReader implements DataSourceReader {

  private GenomicsDBConfiguration genomicsDBConfiguration;

  private long minQueryBlockSize;
  private long maxQueryBlockSize;
  private StructType schema;

  public GenomicsDBDataSourceReader() {
  }

  public GenomicsDBDataSourceReader(StructType schema, DataSourceOptions options) {
    setSchemaOptions(options, schema);
  }

  public GenomicsDBDataSourceReader(DataSourceOptions options) {
    setSchemaOptions(options, null);
  }

  private void setSchemaOptions(DataSourceOptions options,  StructType schema) {
    StructField [] fields = new StructField[] {
                                 new StructField("contig", DataTypes.StringType, false, Metadata.empty()),
                                 new StructField("startPos", DataTypes.IntegerType, false, Metadata.empty()),
                                 new StructField("variantType", DataTypes.StringType, true, Metadata.empty()),
                                 new StructField("alleles", new ArrayType(DataTypes.StringType, true), false, Metadata.empty()),
                                 new StructField("alternateAlleles", new ArrayType(DataTypes.StringType, true), true, Metadata.empty()),
                                 new StructField("sampleNames", new ArrayType(DataTypes.StringType, false), false, Metadata.empty()),
                                 new StructField("genotypes", new ArrayType(DataTypes.StringType, true), true, Metadata.empty())
                            };
    StructType defaultSchema = new StructType(fields);
    // comment out below for now...revisit if needed
    //                             StructField("encodedGenotypes", VectorType, false))

    StructType finalSchema = defaultSchema;
    if (schema!=null) {
      for (StructField sf: JavaConverters.asJavaIterableConverter(schema.toIterable()).asJava()) {
        finalSchema = finalSchema.add(sf);
      }
    }
    genomicsDBConfiguration = new GenomicsDBConfiguration();
    genomicsDBConfiguration.setLoaderJsonFile(options.get(GenomicsDBConfiguration.LOADERJSON).get());
    genomicsDBConfiguration.setQueryJsonFile(options.get(GenomicsDBConfiguration.QUERYJSON).get());
    try {
      genomicsDBConfiguration.setHostFile(options.get(GenomicsDBConfiguration.MPIHOSTFILE).get());
    }
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    this.minQueryBlockSize = options.getLong("genomicsdb.minqueryblocksize", 1);
    this.maxQueryBlockSize = options.getLong("genomicsdb.maxqueryblocksize", Long.MAX_VALUE);
    this.schema = finalSchema;
  }

  public List<InputPartition<InternalRow>> planInputPartitions() {

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

    ArrayList<InputPartition<InternalRow>> inputPartitions = new ArrayList<InputPartition<InternalRow>>();
    // For now, assuming that each of partitions and queryRange are sorted
    // by start position, and that query ranges don't overlap.
    // TODO: not sure anything in GenomicsDB enforces the above assumption....
    int pIndex = 0;
    int qIndex = 0;
    GenomicsDBPartitionInfo partition = null;
    if (!partitionsList.isEmpty())
      partition = partitionsList.get(pIndex);

    // if workspace contains hdfs://, or s3:// or gc:// they're hdfs compliant and we support it
    if (partition != null && !(partition.getWorkspace().contains("s3://") ||
		partition.getWorkspace().contains("hdfs://") ||  
		partition.getWorkspace().contains("gs://"))) {
      List<String> hosts = genomicsDBConfiguration.getHosts();
      for (int i=0; i<hosts.size(); i++) {
        inputPartitions.add(new GenomicsDBInputPartition(hosts.get(i), genomicsDBConfiguration, schema));
      }
    }
    else if (partition != null) {
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
          inputPartitions.add(new GenomicsDBInputPartition(partition, queryRange, genomicsDBConfiguration, schema));
        }
        else {
          // bigger than goalBlockSize, so break up into "query chunks"
  
  	  long queryBlockStart = queryRange.getBeginPosition();
	  long queryBlockMargin = genomicsDBConfiguration.getQueryBlockSizeMargin();
  	  while (queryBlockStart < queryRange.getEndPosition()) {
            long blockSize = (queryBlockSize > (goalBlockSize+queryBlockMargin)) ? goalBlockSize : queryBlockSize;
  	    GenomicsDBQueryInfo queryBlock = new GenomicsDBQueryInfo(queryBlockStart, queryBlockStart + blockSize - 1);
  	    inputPartitions.add(new GenomicsDBInputPartition(partition, queryBlock, genomicsDBConfiguration, schema));
  
  	    // if this queryBlock spans multiple partitions, need to add those as splits as well
  	    while ((pIndex + 1) < partitionsList.size() &&
                    queryBlockStart + blockSize - 1 >= partitionsList.get(pIndex+1).getBeginPosition()) {
  	      pIndex++;
              partition = partitionsList.get(pIndex);
  	      inputPartitions.add(new GenomicsDBInputPartition(partition, queryBlock, genomicsDBConfiguration, schema));
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

  public StructType readSchema() {
    return schema;
  }
}
