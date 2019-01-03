/*
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

package org.genomicsdb.spark;

import org.genomicsdb.reader.GenomicsDBFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import org.apache.hadoop.conf.Configurable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.mapreduce.*;
import org.apache.log4j.Logger;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

public class GenomicsDBInputFormat<VCONTEXT extends Feature, SOURCE>
  extends InputFormat<String, VCONTEXT> implements Configurable {

  private Configuration configuration;
  private GenomicsDBInput input;

  Logger logger = Logger.getLogger(GenomicsDBInputFormat.class);

  /**
   * When this function is called, it is already assumed that configuration
   * object is set
   *
   * @param jobContext  Hadoop Job context passed from newAPIHadoopRDD
   *                    defined in SparkContext
   * @return  Returns a list of input splits
   * @throws FileNotFoundException  Thrown if creaing configuration object fails
   */
  public List<InputSplit> getSplits(JobContext jobContext) throws FileNotFoundException {

    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration(configuration);
    genomicsDBConfiguration.setLoaderJsonFile(
      configuration.get(GenomicsDBConfiguration.LOADERJSON));
    genomicsDBConfiguration.setQueryJsonFile(
      configuration.get(GenomicsDBConfiguration.QUERYJSON));
    genomicsDBConfiguration.setHostFile(
      configuration.get(GenomicsDBConfiguration.MPIHOSTFILE));

    input.setGenomicsDBConfiguration(genomicsDBConfiguration);
    return input.divideInput();
  }

  public RecordReader<String, VCONTEXT>
    createRecordReader(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
      throws IOException, InterruptedException {

    String loaderJson;
    String queryJson;

    GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader;
    GenomicsDBRecordReader<VCONTEXT, SOURCE> recordReader;
    GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)inputSplit;

    if (taskAttemptContext != null) {
      Configuration configuration = taskAttemptContext.getConfiguration();
      loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
      queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
    } else {
      // If control comes here, means this method is called from
      // GenomicsDBRDD. Hence, the configuration object must be
      // set by setConf method, else this will lead to
      // NullPointerException
      assert(configuration!=null);
      loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
      queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
    }

    // Need to amend query file being passed in based on inputSplit
    // so we'll create a temporary query file
    // only do this IF we are using hdfs compliant data store i.e., 
    // getPartitionInfo is not null
    String amendedQuery = queryJson;
    if (gSplit.getPartitionInfo() != null) {
      try {
        amendedQuery = GenomicsDBInput.createTmpQueryFile(queryJson, 
                                                          gSplit.getPartitionInfo(),
                                                          gSplit.getQueryInfo());
      }
      catch (ParseException e) {
        e.printStackTrace();
        return null;
      }
    }
    //GenomicsDBExportConfiguration.ExportConfiguration.Builder exportConfigurationBuilder = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
    //JsonFormat.merge(queryJson, exportConfigurationBuilder);
    //GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration = exportConfigurationBuilder
            //.setWorkspace("").setReferenceGenome("").build();

    featureReader = new GenomicsDBFeatureReader<>(amendedQuery,
            (FeatureCodec<VCONTEXT, SOURCE>) new BCF2Codec(), Optional.of(loaderJson));
    recordReader = new GenomicsDBRecordReader<>(featureReader);
    return recordReader;
  }

  /**
   * default constructor
   */
  public GenomicsDBInputFormat() {
    input = new GenomicsDBInput<GenomicsDBInputSplit>(null, null, 1, Long.MAX_VALUE, GenomicsDBInputSplit.class);
  }

  public GenomicsDBInputFormat(GenomicsDBConfiguration conf) {
    input = new GenomicsDBInput<GenomicsDBInputSplit>(conf, null, 1, Long.MAX_VALUE, GenomicsDBInputSplit.class);
  }

  /**
   * Set the loader JSON file path
   *
   * @param jsonFile  Full qualified path of the loader JSON file
   * @return  Returns the same object for forward function calls
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setLoaderJsonFile(String jsonFile) {
    input.getGenomicsDBConfiguration().setLoaderJsonFile(jsonFile);
    return this;
  }

  /**
   * Set the query JSON file path
   * @param jsonFile  Full qualified path of the query JSON file
   * @return  Returns the same object for forward function calls
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setQueryJsonFile(String jsonFile) {
    input.getGenomicsDBConfiguration().setQueryJsonFile(jsonFile);
    return this;
  }

  /**
   * Set the host file path
   * @param hostFile  Full qualified path of the hosts file
   * @return  Returns the same object for forward function calls
   * @throws FileNotFoundException thrown if the hosts file is not found
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setHostFile(String hostFile)
      throws FileNotFoundException {
    input.getGenomicsDBConfiguration().setHostFile(hostFile);
    return this;
  }

  @Override
  public void setConf(Configuration configuration) {
    this.configuration = configuration;
  }

  @Override
  public Configuration getConf() {
    return configuration;
  }
}
