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
