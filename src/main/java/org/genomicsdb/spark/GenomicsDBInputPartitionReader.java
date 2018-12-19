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

import org.genomicsdb.reader.GenomicsDBFeatureReader;
import org.apache.spark.sql.sources.v2.reader.InputPartitionReader;
import org.apache.spark.sql.catalyst.InternalRow;
import org.apache.spark.sql.catalyst.encoders.*;
import org.apache.spark.sql.RowFactory;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.types.*;
import org.apache.spark.sql.types.DataTypes.*;
import org.apache.spark.sql.*;
import org.apache.spark.ml.linalg.SQLDataTypes;

import htsjdk.variant.variantcontext.*;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.tribble.readers.PositionalBufferedStream;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import scala.Predef;
import scala.collection.Seq;
import scala.collection.JavaConverters;
import java.util.HashMap;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.*;
import java.lang.RuntimeException;

public class GenomicsDBInputPartitionReader implements InputPartitionReader<InternalRow> {

  private GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> fReader;
  private GenomicsDBRecordReader<VariantContext, PositionalBufferedStream> rReader;
  private GenomicsDBInputPartition inputPartition;

  public GenomicsDBInputPartitionReader() {
  }

  public GenomicsDBInputPartitionReader(GenomicsDBInputPartition iPartition) {
    inputPartition = iPartition;
    String loaderJson = inputPartition.getLoaderFile();
    String queryJson = inputPartition.getQueryFile();
    String amendedQuery = queryJson;
    if (inputPartition.getPartitionInfo() != null) {
      try {
        amendedQuery = createTmpQueryFile(queryJson);
      }
      catch (ParseException | IOException e) {
        e.printStackTrace();
      }
    }

    try {
      fReader = new GenomicsDBFeatureReader<>(amendedQuery,
              (FeatureCodec<VariantContext, PositionalBufferedStream>) new BCF2Codec(), Optional.of(loaderJson));
    }
    catch (IOException e) {
      e.printStackTrace();
      fReader = null;
    }
    rReader = new GenomicsDBRecordReader<>(fReader);
    try {
      rReader.initialize(inputPartition);
    }
    catch (IOException e) {
      e.printStackTrace();
    }
  }

  public InternalRow get() {
    VariantContext val;
    try {
      val = (VariantContext)rReader.getCurrentValue();
    }
    catch (IOException | InterruptedException e) {
      e.printStackTrace();
      return null;
    }

    // get ref allele info
    String[] alleles = new String[val.getNAlleles()];
    List<Allele> vcAlleles = val.getAlleles();
    for(int i=0; i<val.getNAlleles(); i++) {
      alleles[i] = vcAlleles.get(i).getBaseString();
    }

    // get alt allele info
    List<Allele> vcAltAlleles = val.getAlternateAlleles();
    String[] altAlleles = null;
    if(vcAltAlleles != null && vcAltAlleles.size() > 0) {
      altAlleles = new String[vcAltAlleles.size()];
      for(int i=0; i<vcAltAlleles.size(); i++) {
        altAlleles[i] = vcAltAlleles.get(i).getBaseString();
      }
    }

    String[] sampleNames = new String[val.getNSamples()];
    String[] genotypes = new String[val.getNSamples()];

    // get sample name and genotype info
    int i=0;
    for (Genotype g: val.getGenotypesOrderedByName()) {
      sampleNames[i] = g.getSampleName();
      genotypes[i++] = g.getGenotypeString();
    }

    //ArrayList<Object> rowObjects = new ArrayList<>(inputPartition.getSchema().size());
    Object[] rowObjects = new Object[inputPartition.getSchema().size()];
    // put down the default schema fields
    rowObjects[0] = val.getContig();
    rowObjects[1] = val.getStart();
    rowObjects[2] = val.getType().toString();
    rowObjects[3] = alleles;
    rowObjects[4] = altAlleles;
    rowObjects[5] = sampleNames;
    rowObjects[6] = genotypes;
    i = 7;
    // go through the schema and try to extract items if we haven't
    // already extracted all that is specified there
    for (StructField sf: JavaConverters.asJavaIterableConverter(inputPartition.getSchema().toIterable()).asJava()) {
      if (sf.name().equals("contig") ||
          sf.name().equals("startPos") ||
          sf.name().equals("variantType") ||
          sf.name().equals("alleles") ||
          sf.name().equals("alternateAlleles") ||
          sf.name().equals("sampleNames") ||
          sf.name().equals("genotypes")) {
        // we've already added these fields
        continue;
      }
      else {
        //getDataFromVariantContext(sf, val, rowObjects);
        rowObjects[i] = getDataFromVariantContext(sf, val);
        i++;
      }
    }
    // in case this is needed...below is the needed transformation to send a Map/HashMap
    // JavaConverters.mapAsScalaMapConverter(genotypes).asScala().toMap(Predef.$conforms()));
    Row row = RowFactory.create(rowObjects);
    ExpressionEncoder<Row> encoder = RowEncoder.apply(inputPartition.getSchema());
    InternalRow iRow = encoder.toRow(row);
    return iRow;
  }

  public boolean next() {
    try {
      return rReader.nextKeyValue();
    }
    catch (IOException | InterruptedException e) {
      e.printStackTrace();
      return false;
    }
  }

  public void close() {
    try {
      rReader.close();
    }
    catch (IOException e) {
      e.printStackTrace();
    }
  }

  private Object getDataFromVariantContext(StructField sf, VariantContext vc) 
    throws RuntimeException {
    if (sf.name().startsWith("INFO_")) {
      String infoname = sf.name().substring(5);
      //rowObjects.add(vc.getAttribute(infoname));
      return vc.getAttribute(infoname);
    } else if (sf.name().startsWith("FORMAT_")) {
      String formatname = sf.name().substring(7);
      Object[] rowObjectArray = new Object[vc.getNSamples()];
      int i=0;
      for (Genotype g: vc.getGenotypesOrderedByName()) {
        rowObjectArray[i++] = g.getAnyAttribute(formatname);
        // TODO: special treatment here for fields/attributes that
        // have arrays? think they show up as arraylists here so might
        // need to do the conversion
      }
      //rowObjects.add(rowObjectArray);
      return rowObjectArray;
    }
    else {
      throw new RuntimeException("Unsupported StructField "+sf.name()+
                                 ". StructField names must start with INFO_ or FORMAT_");
    }
  }

  private String createTmpQueryFile(String queryJson) 
		  throws FileNotFoundException, IOException, ParseException {
    JSONObject amended = new JSONObject();
    amended.put("workspace",inputPartition.getPartitionInfo().getWorkspace());
    amended.put("array",inputPartition.getPartitionInfo().getArrayName());
    amended.put("query_column_ranges",inputPartition.getPartitionInfo().getArrayName());
    if (inputPartition.getQueryInfo().getBeginPosition() == inputPartition.getQueryInfo().getEndPosition()) {
      amended.put("query_column_ranges","[["+inputPartition.getQueryInfo().getBeginPosition()+"]]");
    }
    else {
      amended.put("query_column_ranges","[[["+inputPartition.getQueryInfo().getBeginPosition()+","+
                  inputPartition.getQueryInfo().getEndPosition()+"]]]");
    }

    try {

      JSONParser parser = new JSONParser();
      FileReader queryJsonReader = new FileReader(queryJson);
      JSONObject obj = null;
      try {
        obj = (JSONObject)parser.parse(queryJsonReader);
      }
      catch(ParseException | IOException e) {
        queryJsonReader.close();
        throw e;
      }
  
      obj.forEach((k, v) -> {
        if (!(k.equals("workspace") || k.equals("array"))) {
          amended.put(k, v);
        }
      });
      queryJsonReader.close();
    }
    catch (FileNotFoundException e) {
      e.printStackTrace();
      return null;
    }
    File tmpQueryFile = File.createTempFile("queryJson", ".json");
    tmpQueryFile.deleteOnExit();
    FileWriter fptr = new FileWriter(tmpQueryFile);
    try {
        fptr.write(amended.toJSONString());
    }
    catch(IOException e) {
        fptr.close();
        throw new IOException(e);
    }
    fptr.close();
    return tmpQueryFile.getAbsolutePath();
  }
}
