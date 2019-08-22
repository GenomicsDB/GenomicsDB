package org.genomicsdb.spark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.genomicsdb.reader.GenomicsDBQuery;
import org.genomicsdb.reader.GenomicsDBQuery.*;

import java.lang.reflect.Field;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.List;

public class GenomicsDBSparkBindings {

  static Class getClass(Type type){
    if (type instanceof Class) {
      return type.getClass();
    } else if (type instanceof ParameterizedType) {
      ParameterizedType parameterizedType = (ParameterizedType)type;
      return parameterizedType.getRawType().getClass();
    } else {
      throw new RuntimeException("Cannot find class for type="+type.getTypeName());
    }
  }

  List<VariantCall> variantCalls;

  @SuppressWarnings("unchecked")
  public static void main(String[] args) {
    String loaderJsonFile = args[0];
    String queryJsonFile = args[1];
    String hostfile = args[2];

    SparkConf conf = new SparkConf();
    conf.setAppName("GenomicsDB API Experimental Bindings");
    JavaSparkContext sc = new JavaSparkContext(conf);

    Configuration hadoopConf = sc.hadoopConfiguration();
    hadoopConf.set(GenomicsDBConfiguration.LOADERJSON, loaderJsonFile);
    hadoopConf.set(GenomicsDBConfiguration.QUERYJSON, queryJsonFile);
    hadoopConf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostfile);

    Class variantCallListClass;
    try {
      List<VariantCall> variantCallClass = new ArrayList<>();
      variantCallListClass = Class.forName(variantCallClass.getClass().getGenericSuperclass().getTypeName());
    } catch (ClassNotFoundException e) {
      e.printStackTrace();
      throw new RuntimeException(e.getLocalizedMessage());
    }
    JavaPairRDD<Interval, List<VariantCall>> variants = sc.newAPIHadoopRDD(hadoopConf,
            GenomicsDBAPIInputFormat.class, Interval.class, variantCallListClass);


    System.out.println("Number of variants "+variants.count());
    List variantList = variants.collect();
    for (Object variantObj : variantList) {
      System.out.println(variantObj);
    }

  }
}
