/**
  * The MIT License (MIT)
  * Copyright (c) 2016-2017 Intel Corporation
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

package org.genomicsdb

import org.genomicsdb.spark.{GenomicsDBConfiguration, GenomicsDBInputFormat}
import org.genomicsdb.spark.GenomicsDBInputFormat
import htsjdk.tribble.readers.PositionalBufferedStream
import htsjdk.variant.variantcontext.VariantContext
import org.apache.spark.rdd.RDD
import org.apache.spark.{SparkConf, SparkContext}
import org.genomicsdb.reader.GenomicsDBQuery.{Interval, VariantCall}
import org.genomicsdb.spark.api.GenomicsDBQueryInputFormat

/**
  * This factory class exposes multiple ways of retrieving RDD of either variant contexts (htsjdk)
  * or List of VariantCalls from GenomicsDB
  */
object GenomicsDBScalaSparkFactory {

  def usingGenomicsRDD(args: Array[String]): Unit = {

    val master = args(0)
    val port = args(1)
    val loaderJsonFile = args(2)
    val queryJsonFile = args(3)
    val hostfile = args(4)

    val conf = new SparkConf()
      .setMaster("spark://" + master + ":" + port)
      .setAppName("GenomicsDBTest using GenomicsDBRDD")
    val sc = new SparkContext(conf)
    val hadoopConf = sc.hadoopConfiguration
    hadoopConf.set(GenomicsDBConfiguration.LOADERJSON, loaderJsonFile)
    hadoopConf.set(GenomicsDBConfiguration.QUERYJSON, queryJsonFile)
    hadoopConf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostfile)

    val gc = new GenomicsDBContext(hadoopConf, sc)
    val myrdd = gc.getVariantContexts
    System.out.println(myrdd.count())
    myrdd.collect().foreach(println)
  }

  def usingNewAPIHadoopRDD(args: Array[String]): Unit = {
    val master = args(0)
    val port = args(1)
    val loaderJsonFile = args(2)
    val queryJsonFile = args(3)
    val hostfile = args(4)
    val conf = new SparkConf()
      .setMaster("spark://" + master + ":" + port)
      .setAppName("GenomicsDBTest using newAPIHadoopRDD")

//    val classObjects: Array[Class[_]] = Array(Class.forName("org.apache.hadoop.io.LongWritable"),
//      Class.forName("org.apache.hadoop.io.Text"))
//    conf.registerKryoClasses(classObjects)

    val sc = new SparkContext(conf)
    val hadoopConf = sc.hadoopConfiguration
    hadoopConf.set(GenomicsDBConfiguration.LOADERJSON, loaderJsonFile)
    hadoopConf.set(GenomicsDBConfiguration.QUERYJSON, queryJsonFile)
    hadoopConf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostfile)

    val myrdd: RDD[(String, VariantContext)] =
      sc.newAPIHadoopRDD[String, VariantContext,
      GenomicsDBInputFormat[VariantContext, PositionalBufferedStream]](sc.hadoopConfiguration,
      classOf[GenomicsDBInputFormat[VariantContext, PositionalBufferedStream]],
      classOf[String], classOf[VariantContext])

    System.out.println(myrdd.count())
  }

  def usingGenomicsDBQuery(args: Array[String]): Unit = {
    val conf = new SparkConf().setAppName("GenomicsDBQuery API Test with Spark")

    if (args.length == 4) {
      val master = args(0)
      val port = args(1)
      conf.setMaster("spark://" + master + ":" + port)
    }
    val loaderJsonFile = if (args.length==4) args(2) else args(0)
    val queryJsonFile = if (args.length==4) args(2) else args(1)

    val sc = new SparkContext(conf)

    val hadoopConf = sc.hadoopConfiguration
    hadoopConf.set(GenomicsDBConfiguration.LOADERJSON, loaderJsonFile)
    hadoopConf.set(GenomicsDBConfiguration.QUERYJSON, queryJsonFile)

    val intervalRDDs: RDD[(Interval, List[VariantCall])] = sc.newAPIHadoopRDD(hadoopConf,
      classOf[GenomicsDBQueryInputFormat], classOf[Interval], classOf[List[VariantCall]])

    val count = intervalRDDs.count()
    System.out.println("Count=" + count)

    val intervals: Array[(Interval, List[VariantCall])] = intervalRDDs.collect()
    intervals.foreach(println)
  }

  def main(args: Array[String]): Unit = {
    //   usingNewAPIHadoopRDD(args)
    //   usingGenomicsRDD(args)
    usingGenomicsDBQuery(args)
  }
}
