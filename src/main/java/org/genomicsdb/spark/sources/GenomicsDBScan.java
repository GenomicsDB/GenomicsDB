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
 * */

package org.genomicsdb.spark.sources;

import org.apache.spark.sql.connector.read.Batch;
import org.apache.spark.sql.connector.read.Scan;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

import java.util.Map;

/**
 * Scan is a logical plan of a GenomicsDB Scan.
 **/
public class GenomicsDBScan implements Scan {

  private final StructType schema;
  private final Map<String, String> properties;
  private final CaseInsensitiveStringMap options;

  public GenomicsDBScan(StructType schema, Map<String, String> properties,
    CaseInsensitiveStringMap options){

    this.schema = schema;
    this.properties = properties;
    this.options = options;

  }

  @Override 
  public StructType readSchema(){
    return schema;
  }

  @Override
  public String description(){
    return "GenomicsDBSource Scan";
  }

  @Override 
  public Batch toBatch(){
    return new GenomicsDBBatch(schema, properties, options);
  }

}
