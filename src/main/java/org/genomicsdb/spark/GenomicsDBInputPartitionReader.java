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
import java.lang.NumberFormatException;
import java.lang.NullPointerException;
import java.lang.reflect.Array;

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
    // put down the default schema fields
    rowObjects.add(UTF8String.fromString(val.getContig()));
    rowObjects.add(val.getStart());
    rowObjects.add(UTF8String.fromString(val.getID()));
    rowObjects.add(UTF8String.fromString(val.getType().toString()));
    rowObjects.add(UTF8String.fromString(val.getReference().getBaseString()));
    if(altAlleles == null) {
      rowObjects.add(null);
    }
    else {
      rowObjects.add(ArrayData.toArrayData(altAlleles));
    }
    rowObjects.add(ArrayData.toArrayData(sampleNames));
    rowObjects.add(ArrayData.toArrayData(genotypes));

    // go through the schema and try to extract items if we haven't
    // already extracted all that is specified there
    for (StructField sf: JavaConverters.asJavaIterableConverter(inputPartition.getSchema().toIterable()).asJava()) {
      if (sf.name().equals("contig") ||
          sf.name().equals("startPos") ||
          sf.name().equals("ID") ||
          sf.name().equals("variantType") ||
          sf.name().equals("refAllele") ||
          sf.name().equals("alternateAlleles") ||
          sf.name().equals("sampleNames") ||
          sf.name().equals("GT")) {
        // we've already added these fields
        continue;
      }
      else {
        Map<String, GenomicsDBVidSchema> vMap = inputPartition.getGenomicsDBVidSchema();
        GenomicsDBVidSchema field = vMap.get(sf.name());
        // TODO: add logging and mention couldn't find field....
        if (field == null) {
          rowObjects.add(null);
          continue;
        }

        // if INFO field and single object per VariantContext
        if (field.isInfo() && (field.getLength().equals("") || 
                               field.getLength().equals("1"))) {
          rowObjects.add(getDataFromVariantContext(sf.name(), val, field));
        }
        else if (field.isInfo() && field.getFieldClass().equals(String.class)) {
          rowObjects.add(UTF8String.fromString((String)getDataFromVariantContext(sf.name(), val, field)));
        }
        // FORMAT fields and INFO with multiple objects per VC
        // will use ArrayType
        else {
          System.out.println("SF name is:"+sf.name());
          Object obj = getDataFromVariantContext(sf.name(), val, field);
          if (obj!= null) {
            rowObjects.add(ArrayData.toArrayData(obj));
          }
          else {
            rowObjects.add(null);
          }
        }


/*
        if (sf.name().startsWith("INFO_")) {
          rowObjects.add(getDataFromVariantContext(sf, val));
        }
        else if (sf.name().startsWith("FORMAT_")) {
          rowObjects.add(ArrayData.toArrayData(getDataFromVariantContext(sf, val)));
        }
*/
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

  private int getFieldLength(VariantContext vc, GenomicsDBVidSchema field) 
    throws RuntimeException {
    String length = field.getLength();
    int l;
    try {
      l = Integer.parseInt(length);
    } catch (NumberFormatException | NullPointerException e) {
      switch (length) {
        case "A": l = vc.getAlternateAlleles().size(); break;
        case "R": l = vc.getAlleles().size(); break;
        case "G": l = vc.getAlleles().size()*(vc.getAlleles().size()+1)/2; break;
        default: throw new RuntimeException("Unsupported length field "+length+
                                            " in vid mapping");
      }
    }
    return l;
  }

  private Object[] createObjectArray(VariantContext vc, GenomicsDBVidSchema field) {
    int l = getFieldLength(vc, field);
    return (Object[])Array.newInstance(field.getFieldClass(), l);
  }

  private Object getDataFromVariantContext(String fName, VariantContext vc, GenomicsDBVidSchema field) 
    throws RuntimeException {
    if (field.isInfo() && (field.getLength().equals("1") ||
                           field.getFieldClass().equals(String.class))) {
      return vc.getAttribute(fName);
    }
    else {
      Object[] rowObjectArray;
      // these fields are just Arrays of elements
      // either INFO fields that are arrays...
      if (field.isInfo()) { 
        rowObjectArray = vc.getAttributeAsList(fName)
                           .toArray((Object[])Array.newInstance(field.getFieldClass(), 0));
/*
        rowObjectArray = createObjectArray(vc, field);
        int i=0;
        for (Genotype g: vc.getGenotypesOrderedByName()) {
          rowObjectArray[i++] = field.getFieldClass().cast(g.getAnyAttribute(fName));
        }
*/
      }
      // ...Or FORMAT/genotype fields that have a single object per sample
      else if (!field.isInfo() && (field.getLength().equals("1") ||
                              field.getFieldClass().equals(String.class))) {
        if (field.getFieldClass().equals(String.class)) {
          rowObjectArray = (Object[])Array.newInstance(UTF8String.class,
                                                       vc.getNSamples());
        }
        else {
          rowObjectArray = (Object[])Array.newInstance(field.getFieldClass(),
                                                       vc.getNSamples());
        }
        // sorta hacky but...in cases where we've added the _FORMAT
        // suffix, we need to remove it before querying for that attribute
        String noSuffix = fName;
        if (fName.endsWith("_FORMAT")) {
          noSuffix = fName.substring(0, fName.length()-7);
        }
        int i=0;
        for (Genotype g: vc.getGenotypesOrderedByName()) {
          if (g.hasAnyAttribute(noSuffix)) {
            if (field.getFieldClass().equals(String.class)) {
              rowObjectArray[i++] = UTF8String.fromString(
                                                    (String)((g.getAnyAttribute(noSuffix))));
            }
            else {
              rowObjectArray[i++] = field.getFieldClass().cast((g.getAnyAttribute(noSuffix)));
            }
          }
          else {
            rowObjectArray[i++] = null;
          }
        }
      }
      // these fields are Arrays of Arrays
      // basically FORMAT/genotype fields that are also arrays
      else {
        rowObjectArray = (Object[])Array.newInstance(ArrayData.class, 
                                                     vc.getNSamples());

        int i=0;
        for (Genotype g: vc.getGenotypesOrderedByName()) {
          if (g.hasAnyAttribute(fName)) {
/*
            int l = getFieldLength(vc, field);
            Object[] tempObj = (Object[])Array.newInstance(field.getFieldClass(), l);
            List attrArray = (List)(g.getAnyAttribute(fName));
            for(int j=0; j<l; j++) {
              tempObj[j] = (Object)attrArray.get(j);
            }
            rowObjectArray[i++] = ArrayData.toArrayData(tempObj);
*/
            rowObjectArray[i++] = ArrayData.toArrayData(((List)(g.getAnyAttribute(fName)))
                                         .toArray((Object[])Array.newInstance(field.getFieldClass(), 0)));
          }
          else {
            rowObjectArray[i++] = null;
          }
        }
      }
      return rowObjectArray;
    }

/*
    if (sf.name().startsWith("INFO_")) {
      String infoname = sf.name().substring(5);
      return vc.getAttribute(infoname);
    } 
    else if (sf.name().startsWith("FORMAT_")) {
      String formatname = sf.name().substring(7);
      if (!(sf.dataType() instanceof ArrayType)) {
        throw new RuntimeException("Unsupported dataype: schema field "+sf.name()+
             " has type "+sf.dataType().simpleString()+" but must be ArrayType");
      }

      Object[] rowObjectArray;
      DataType eType = ((ArrayType)sf.dataType()).elementType();
      if (eType instanceof IntegerType) {
        rowObjectArray = new Integer[vc.getNSamples()];
      }
      else if (eType instanceof LongType) {
        rowObjectArray = new Long[vc.getNSamples()];
      }
      else if (eType instanceof DoubleType) {
        rowObjectArray = new Double[vc.getNSamples()];
      }
      else if (eType instanceof FloatType) {
        rowObjectArray = new Float[vc.getNSamples()];
      }
      else if (eType instanceof StringType) {
        rowObjectArray = new String[vc.getNSamples()];
      }
      else {
        throw new RuntimeException("Unsupported element type "+
              ((ArrayType)sf.dataType()).elementType().simpleString()+
              " for schema field "+sf.name());
      }
      int i=0;
      
      for (Genotype g: vc.getGenotypesOrderedByName()) {
        Class<?> cType = rowObjectArray.getClass().getComponentType();
        rowObjectArray[i++] = cType.cast(g.getAnyAttribute(formatname));
        // TODO: special treatment here for fields/attributes that
        // have arrays? think they show up as arraylists here so might
        // need to do the conversion
        // Not sure we'll have any fields like that so ignore for now
      }
      return rowObjectArray;
    }
    else {
      throw new RuntimeException("Unsupported StructField "+sf.name()+
                                 ". StructField names must start with INFO_ or FORMAT_");
    }
*/
  }
}
