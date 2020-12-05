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

import org.apache.spark.sql.sources.BaseRelation;
import org.apache.spark.sql.sources.RelationProvider;
import org.apache.spark.sql.sources.SchemaRelationProvider;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.sources.DataSourceRegister;
import org.apache.spark.sql.SQLContext;

public class GenomicsDBDataSource implements RelationProvider, SchemaRelationProvider, DataSourceRegister {

  public GenomicsDBDataSource() {}

  @Override
  public BaseRelation createRelation(SQLContext sqlContext, scala.collection.immutable.Map<String,String> parameters){
    return new GenomicsDBRelation(sqlContext, parameters);
  }

  @Override
  public BaseRelation createRelation(SQLContext sqlContext, scala.collection.immutable.Map<String,String> parameters, 
      StructType schema){
    return new GenomicsDBRelation(sqlContext, schema, parameters);
  }

  @Override
  public String shortName() {
    // TODO need to add stuff to META-INF/services to make this work
    return "genomicsdb";
  }
}
