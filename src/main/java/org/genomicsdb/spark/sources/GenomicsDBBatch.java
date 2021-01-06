/*
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation
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

package org.genomicsdb.spark.sources;

import org.apache.spark.sql.catalyst.InternalRow;
import org.apache.spark.sql.connector.read.Batch;
import org.apache.spark.sql.connector.read.InputPartition;
import org.apache.spark.sql.connector.read.PartitionReaderFactory;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

import org.apache.spark.sql.types.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import scala.collection.JavaConverters;

import org.genomicsdb.spark.GenomicsDBInput;
import org.genomicsdb.spark.GenomicsDBConfiguration;
import org.genomicsdb.spark.GenomicsDBSchemaFactory;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Represents a physical plan, logical table scan is converted to this at runtime.
 * Defines a factory that is sent to an executor, thus Datasource partitions are planned here.
 **/
public class GenomicsDBBatch implements Batch {

  GenomicsDBInput<GenomicsDBInputPartition> input;
  private final Map<String, String> properties;
  private final CaseInsensitiveStringMap options;

  public GenomicsDBBatch(StructType schema, Map<String, String> properties, 
      CaseInsensitiveStringMap options) {
    this.properties = properties;
    this.options = options;
    setSchemaOptions(options, schema);
  }

  private void setSchemaOptions(CaseInsensitiveStringMap options, StructType schema)
      throws RuntimeException {

    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration((Map<String,String>)options);

    GenomicsDBSchemaFactory schemaBuilder = 
      new GenomicsDBSchemaFactory(genomicsDBConfiguration.getLoaderJsonFile());
    StructType finalSchema = null;
    if (schema != null){ 
      finalSchema = schema;
    } else { 
      finalSchema = schemaBuilder.buildSchemaWithVid(schema.fields()); 
    }

    Long blocksize = new Long(1);
    Long maxblock = Long.MAX_VALUE;
    if (options.get("genomicsdb.minqueryblocksize") != null){
      blocksize = Long.valueOf(options.get("genomicsdb.minqueryblocksize"));
    }
    if (options.get("genomicsdb.maxqueryblocksize") != null){
      maxblock = Long.valueOf(options.get("genomicsdb.maxqueryblocksize"));
    }
    // TODO what is vMap used for here...
    input = new GenomicsDBInput<>(
            genomicsDBConfiguration,
            finalSchema,
            schemaBuilder.getVidMap(),
            blocksize,
            maxblock,
            GenomicsDBInputPartition.class);
  }

  @Override
  //@SuppressWarnings("unchecked")
  public InputPartition[] planInputPartitions(){
    List<GenomicsDBInputPartition> partitionsList = input.divideInput();
    GenomicsDBInputPartition[] ipartitions = new GenomicsDBInputPartition[partitionsList.size()];
    ipartitions = partitionsList.toArray(ipartitions);
    return ipartitions;
  }
  
  @Override
  public PartitionReaderFactory createReaderFactory(){
    return new GenomicsDBPartitionReaderFactory(); 
  }  

}
