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

  GenomicsDBInput input;

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
                                 new StructField("id", DataTypes.StringType, true, Metadata.empty()),
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
    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration();
    genomicsDBConfiguration.setLoaderJsonFile(options.get(GenomicsDBConfiguration.LOADERJSON).get());
    genomicsDBConfiguration.setQueryJsonFile(options.get(GenomicsDBConfiguration.QUERYJSON).get());
    try {
      genomicsDBConfiguration.setHostFile(options.get(GenomicsDBConfiguration.MPIHOSTFILE).get());
    }
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    input = new GenomicsDBInput<GenomicsDBInputPartition>(genomicsDBConfiguration, finalSchema,
                                options.getLong("genomicsdb.minqueryblocksize", 1),
                                options.getLong("genomicsdb.maxqueryblocksize", Long.MAX_VALUE),
                                GenomicsDBInputPartition.class);
  }

  public List<InputPartition<InternalRow>> planInputPartitions() {
    return input.divideInput();
  }

  public StructType readSchema() {
    return input.getSchema();
  }
}
