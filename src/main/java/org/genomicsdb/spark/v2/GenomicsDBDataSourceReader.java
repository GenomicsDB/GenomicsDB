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

package org.genomicsdb.spark.v2;

import org.apache.spark.sql.catalyst.InternalRow;
import org.apache.spark.sql.sources.v2.DataSourceOptions;
import org.apache.spark.sql.sources.v2.reader.DataSourceReader;
import org.apache.spark.sql.sources.v2.reader.InputPartition;
import org.apache.spark.sql.types.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import scala.collection.JavaConverters;

import org.genomicsdb.spark.GenomicsDBConfiguration;
import org.genomicsdb.spark.GenomicsDBSchemaFactory;
import org.genomicsdb.spark.GenomicsDBInput;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GenomicsDBDataSourceReader implements DataSourceReader {

  GenomicsDBInput<GenomicsDBInputPartition> input;

  public GenomicsDBDataSourceReader() {}

  public GenomicsDBDataSourceReader(StructType schema, DataSourceOptions options) {
    setSchemaOptions(options, schema);
  }

  public GenomicsDBDataSourceReader(DataSourceOptions options) {
    setSchemaOptions(options, null);
  }

  private void setSchemaOptions(DataSourceOptions options, StructType schema)
      throws RuntimeException {

    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration(options.asMap());
    GenomicsDBSchemaFactory schemaBuilder = 
      new GenomicsDBSchemaFactory(genomicsDBConfiguration);
    StructType finalSchema = null;
    if (schema != null){ 
      finalSchema = schema;
    } else { 
      finalSchema = schemaBuilder.buildSchemaWithVid(schema.fields()); 
    }
    input = new GenomicsDBInput<>(
              genomicsDBConfiguration,
              finalSchema,
              schemaBuilder.getVidMap(),
              options.getLong("genomicsdb.minqueryblocksize", 1),
              options.getLong("genomicsdb.maxqueryblocksize", Long.MAX_VALUE),
              GenomicsDBInputPartition.class);
  }

  @Override
  @SuppressWarnings("unchecked")
  public List<InputPartition<InternalRow>> planInputPartitions() {
    return (List)input.divideInput();
  }

  @Override
  public StructType readSchema() {
    return input.getSchema();
  }
}
