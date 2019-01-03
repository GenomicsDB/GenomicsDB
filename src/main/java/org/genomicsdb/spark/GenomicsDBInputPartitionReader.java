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
import org.apache.spark.sql.catalyst.util.ArrayData;
import org.apache.spark.sql.RowFactory;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.types.*;
import org.apache.spark.sql.types.DataTypes.*;
import org.apache.spark.sql.*;
import org.apache.spark.ml.linalg.SQLDataTypes;
import org.apache.spark.unsafe.types.UTF8String;

import htsjdk.variant.variantcontext.*;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.CloseableTribbleIterator;

import org.json.simple.JSONObject;
import org.json.simple.JSONArray;
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
import java.util.*;
import java.lang.RuntimeException;

public class GenomicsDBInputPartitionReader implements InputPartitionReader<InternalRow> {

  private GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> fReader;
  private CloseableTribbleIterator<VariantContext> iterator;
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
        amendedQuery = GenomicsDBInput.createTmpQueryFile(queryJson,
                                                          inputPartition.getPartitionInfo(),
                                                          inputPartition.getQueryInfo());
      }
      catch (ParseException | IOException e) {
        e.printStackTrace();
      }
    }

    try {
      fReader = new GenomicsDBFeatureReader<>(amendedQuery,
              (FeatureCodec<VariantContext, PositionalBufferedStream>) new BCF2Codec(), Optional.of(loaderJson));
      this.iterator = fReader.iterator();
    }
    catch (IOException e) {
      e.printStackTrace();
      fReader = null;
    }
  }

  public InternalRow get() {
    VariantContext val;
    val = (VariantContext)this.iterator.next();

    // get ref allele info
    UTF8String[] alleles = new UTF8String[val.getNAlleles()];
    List<Allele> vcAlleles = val.getAlleles();
    for(int i=0; i<val.getNAlleles(); i++) {
      alleles[i] = UTF8String.fromString(vcAlleles.get(i).getBaseString());
    }

    // get alt allele info
    List<Allele> vcAltAlleles = val.getAlternateAlleles();
    UTF8String[] altAlleles = null;
    if(vcAltAlleles != null && vcAltAlleles.size() > 0) {
      altAlleles = new UTF8String[vcAltAlleles.size()];
      for(int i=0; i<vcAltAlleles.size(); i++) {
        altAlleles[i] = UTF8String.fromString(vcAltAlleles.get(i).getBaseString());
      }
    }

    UTF8String[] sampleNames = new UTF8String[val.getNSamples()];
    UTF8String[] genotypes = new UTF8String[val.getNSamples()];

    // get sample name and genotype info
    int i=0;
    for (Genotype g: val.getGenotypesOrderedByName()) {
      sampleNames[i] = UTF8String.fromString(g.getSampleName());
      genotypes[i++] = UTF8String.fromString(g.getGenotypeString());
    }

    ArrayList<Object> rowObjects = new ArrayList<>(inputPartition.getSchema().size());
    //Object[] rowObjects = new Object[inputPartition.getSchema().size()];
    // put down the default schema fields
    rowObjects.add(UTF8String.fromString(val.getContig()));
    rowObjects.add(val.getStart());
    rowObjects.add(UTF8String.fromString(val.getType().toString()));
    rowObjects.add(ArrayData.toArrayData(alleles));
    rowObjects.add(ArrayData.toArrayData(altAlleles));
    rowObjects.add(ArrayData.toArrayData(sampleNames));
    rowObjects.add(ArrayData.toArrayData(genotypes));

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
        rowObjects.add(getDataFromVariantContext(sf, val));
      }
    }
    // in case this is needed...below is the needed transformation to send a Map/HashMap
    // JavaConverters.mapAsScalaMapConverter(genotypes).asScala().toMap(Predef.$conforms()));
    //Row row = RowFactory.create(rowObjects);
    //ExpressionEncoder<Row> encoder = RowEncoder.apply(inputPartition.getSchema());
    //InternalRow iRow = encoder.toRow(row);
    InternalRow iRow = InternalRow.fromSeq(JavaConverters.asScalaIteratorConverter(rowObjects.iterator()).asScala().toSeq());
    return iRow;
  }

  public boolean next() {
    return this.iterator.hasNext();
  }

  public void close() {
    this.iterator.close();
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
}
