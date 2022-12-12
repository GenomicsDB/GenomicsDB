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

import org.apache.spark.sql.connector.catalog.Table;
import org.apache.spark.sql.connector.catalog.TableProvider;
import org.apache.spark.sql.connector.expressions.Transform;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.types.*;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

import org.genomicsdb.spark.GenomicsDBSchemaFactory;

import java.util.Map;

/** 
 * The base interface for all custom datasources in Spark3. 
 * Replaces the primary data source interface in the DatasourceV2 API, 
 * which was removed as of Spark 3.0.x and replaced with TableProvider interface. 
 * This is the base interface for all custom datasource. 
 * 
 **/

public class GenomicsDBSource implements TableProvider {

  public GenomicsDBSource(){}

  /** 
   * This sets the default schema, which will contain all 
   * the default fields coming 
   * from the backend. 
   *
   **/
  public StructType inferSchema(CaseInsensitiveStringMap options){
    return GenomicsDBSchemaFactory.defaultSchema();
  }

  @Override
  public Table getTable(StructType schema, Transform[] partitioning, Map<String, String> properties){
    return new GenomicsDBTable(schema, properties);
  }

  @Override
  public boolean supportsExternalMetadata(){
    return true;
  }

}
