/*
 * The MIT License (MIT)
 * Copyright (c) 2020, 2023 Omics Data Automation, Inc.
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

package org.genomicsdb.spark.api;

import com.google.protobuf.util.JsonFormat;
import org.apache.commons.io.FileUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.reader.GenomicsDBQuery.Interval;
import org.genomicsdb.reader.GenomicsDBQuery.VariantCall;
import org.genomicsdb.spark.GenomicsDBConfiguration;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Base64;
import java.util.List;

/**
 * Example Invocation
 * spark-submit --class org.genomicsdb.spark.api.GenomicsDBSparkBindings genomicsdb-1.3.1-SNAPSHOT-allinone.jar loader.json querypb.json true
 *      querypb.json should be parseable by GenomicsDBExportConfiguration.ExportConfiguration
 *  OR
 *  spark-submit --class org.genomicsdb.spark.api.GenomicsDBSparkBindings genomicsdb-1.3.1-SNAPSHOT-allinone.jar loader.json query.json false
 *  OR
 *  spark-submit --class org.genomicsdb.spark.api.GenomicsDBSparkBindings genomicsdb-1.3.1-SNAPSHOT-allinone.jar loader.json query.json
 */
public class GenomicsDBSparkBindings {
  List<VariantCall> variantCalls;

  public static void main(String[] args) throws IOException, ClassNotFoundException {
    if (args.length < 2) {
      throw new RuntimeException("Usage: spark-submit --class org.genomicsdb.spark.api.GenomicsDBSparkBindings genomicsdb-<VERSION>-allinone.jar <loader.json> <query.json> [<is_serialized_pb>]"+
              "Optional Argument 2 - <is_serialized_pb=True|False, default is false, if is_serialized_pb then query.json is a protobuf serialized file.");
    }

    String loaderJsonFile = args[0];
    String queryJsonFile = args[1];
    boolean isPB = false;
    if (args.length == 3) {
      isPB = new Boolean(args[2]).booleanValue();
    }

    SparkConf conf = new SparkConf();
    conf.setAppName("GenomicsDB API Experimental Bindings");
    JavaSparkContext sc = new JavaSparkContext(conf);

    Configuration hadoopConf = sc.hadoopConfiguration();
    if (!loaderJsonFile.isEmpty()) {
      hadoopConf.set(GenomicsDBConfiguration.LOADERJSON, loaderJsonFile);
    }

    if (isPB) {
      String queryPBString = FileUtils.readFileToString(new File(queryJsonFile));
      final GenomicsDBExportConfiguration.ExportConfiguration.Builder builder = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
      JsonFormat.parser().merge(queryPBString, builder);
      queryPBString = Base64.getEncoder().encodeToString(builder.build().toByteArray());
      hadoopConf.set(GenomicsDBConfiguration.QUERYPB, queryPBString);
    } else {
      hadoopConf.set(GenomicsDBConfiguration.QUERYJSON, queryJsonFile);
    }

    Class variantCallListClass = Class.forName("java.util.List");
    JavaPairRDD<Interval, List<VariantCall>> variants = sc.newAPIHadoopRDD(hadoopConf,
            GenomicsDBQueryInputFormat.class, Interval.class, variantCallListClass);

    System.out.println("Number of variants "+variants.count());
    List variantList = variants.collect();
    for (Object variantObj : variantList) {
      System.out.println(variantObj);
    }
  }
}
