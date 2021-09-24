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

package org.genomicsdb.spark;

import org.apache.spark.sql.types.StructType;
import org.genomicsdb.model.GenomicsDBVidMapProto;
import org.genomicsdb.model.GenomicsDBVidMapProto.VidMappingPB;
import org.apache.spark.sql.types.*;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import scala.collection.JavaConverters;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Collection;
import java.util.Set;
import java.util.HashSet;

/**
 * Defines the default schema and provides methods to extend upon the default schema. 
 * Can build a vid map that is used for GenomicsDBInput, and will also ensure that types 
 * are validated using DataType DDL.
 *
 * When a vid map is needed, an instance should be defined with a constructor. When 
 * the vid map is not needed, such as from the scala api - then this is treated as a 
 * utility class.
 **/
public class GenomicsDBSchemaFactory {

  private final Map<String, GenomicsDBVidSchema> vidMap;

  public GenomicsDBSchemaFactory(String loader){
    this.vidMap = buildVidSchemaMap(loader);
  }

  public GenomicsDBSchemaFactory(VidMappingPB vidmapPB) {
    this.vidMap = new HashMap<String, GenomicsDBVidSchema>();
    for (GenomicsDBVidMapProto.GenomicsDBFieldInfo f : vidmapPB.getFieldsList()) {
      String key = f.getName();
      Class<?> clazz = getFieldType(f.getType(0));
      String length = getFieldLength(f);
      if (f.getVcfFieldClassCount() == 1) {
        this.vidMap.put((String) key, new GenomicsDBVidSchema(f.getVcfFieldClass(0).equals("INFO"), clazz, length));
      }
      // if field is both, add entries for INFO under it's name
      // and then add an entry for FORMAT by adding the suffix _FORMAT
      // mainly this is for DP
      else {
        this.vidMap.put((String) key, new GenomicsDBVidSchema(true, clazz, length));
        this.vidMap.put((String) key + "_FORMAT", new GenomicsDBVidSchema(false, clazz, length));
      }
    }
  }

  private Class<?> getFieldType(String fieldtype) {
    // TODO: we don't currently support multiple types for field
    switch (fieldtype) {
      case "int":
      case "integer":
      case "Integer":
        return Integer.class;
      case "char":
      case "String":
        return String.class;
      case "float":
      case "Float":
      case "double":
      case "Double":
        // TODO: switching to double below because that is what VC
        // seems to be giving us. Otherwise spark will complain later
        // when we try to create the dataframe about unboxing from double to float
        return Double.class;
      case "flag":
        return Boolean.class;
      default:
        throw new RuntimeException("Unsupported type " + fieldtype + " in vid mapping");
    }
  }

  private String getFieldLength(GenomicsDBVidMapProto.GenomicsDBFieldInfo field) {
    // TODO: we don't currently support multiple lengths for field
    if (field.getLengthCount() == 0) {
      return "1";
    }
    if (field.getLength(0).hasFixedLength()) {
      return Integer.toString(field.getLength(0).getFixedLength());
    }
    else {
      return field.getLength(0).getVariableLengthDescriptor();
    }
  }
   
  public Map<String, GenomicsDBVidSchema> getVidMap(){
    return vidMap;
  }
 
  public static List<StructField> defaultFields(){
    List<StructField> fields = new ArrayList<>(); 
    fields.add(DataTypes.createStructField("contig", DataTypes.StringType, false));
    fields.add(DataTypes.createStructField("startPos", DataTypes.IntegerType, false));
    fields.add(DataTypes.createStructField("ID", DataTypes.StringType, true));
    fields.add(DataTypes.createStructField("variantType", DataTypes.StringType, true));
    fields.add(DataTypes.createStructField("refAllele", DataTypes.StringType, true));
    fields.add(DataTypes.createStructField("alternateAlleles", 
          DataTypes.createArrayType(DataTypes.StringType, true), true));
    fields.add(DataTypes.createStructField("sampleNames", 
          DataTypes.createArrayType(DataTypes.StringType, false), false));
    fields.add(DataTypes.createStructField("GT", 
          DataTypes.createArrayType(DataTypes.StringType, true), true));
    return fields;
  }

  public static List<String> defaultFieldNames() {
    List<String> names = new ArrayList<>();
    for (StructField f: defaultFields()) {
      names.add(f.name());
    }

    return names;
  }
 
  public static StructType defaultSchema(){
    return DataTypes.createStructType(defaultFields());
  }

  private List<StructField> addFieldsWithVid(StructField[] addFields){
    List<StructField> fields = new ArrayList<>();
    List<String> defaults = defaultFieldNames();
    fields.addAll(defaultFields());
    for (StructField f: addFields){
      if (defaults.contains(f.name())) {
        continue;
      }
      GenomicsDBVidSchema field = this.vidMap.get(f.name());
      if (field == null) {
        throw new RuntimeException("Schema field " + f.name() + " not available in vid");      
      }
      fields.add(DataTypes.createStructField(f.name(), DataType.fromDDL(field.getDataType()), f.nullable()));
    }
    return fields;
  }

  public StructType buildSchemaWithVid(StructField[] addFields){
    return DataTypes.createStructType(addFieldsWithVid(addFields));
  }

  private Map<String, GenomicsDBVidSchema> buildVidSchemaMap(String loader){
    Map<String, GenomicsDBVidSchema> vMap = new HashMap<String, GenomicsDBVidSchema>();
    String vidMapping = "";
    try {
      FileReader loaderReader = new FileReader(loader);
      JSONParser parser = new JSONParser();
      JSONObject objLoad = (JSONObject) parser.parse(loaderReader);

      vidMapping = (String) objLoad.get("vid_mapping_file");
      loaderReader.close();

      FileReader vidReader = new FileReader(vidMapping);
      JSONObject objVid = (JSONObject) parser.parse(vidReader);

      HashMap<?,?> fields = (HashMap<?,?>) objVid.get("fields");
      fields.forEach(
        (k, vObj) -> {
          // ignore fields without vcf_field_class
          HashMap<?,?> v = (HashMap<?,?>) vObj;
          JSONArray fieldClass = (JSONArray) v.get("vcf_field_class");
          if (fieldClass != null) {
            String vType = (String) v.get("type");
            Class<?> clazz = getFieldType(vType);
            
            String length;
            if(v.get("length") == null) {
              length = "1";
            } else {
              length = v.get("length").toString();
            }
            // if field is INFO or FORMAT
            if (fieldClass.size() == 1) {
              vMap.put((String) k, new GenomicsDBVidSchema(fieldClass.get(0).equals("INFO"), clazz, length));
            }
            // if field is both, add entries for INFO under it's name
            // and then add an entry for FORMAT by adding the suffix _FORMAT
            // mainly this is for DP
            else {
              vMap.put((String) k, new GenomicsDBVidSchema(true, clazz, length));
              vMap.put((String) k + "_FORMAT", new GenomicsDBVidSchema(false, clazz, length));
            }
          }
        });
      vidReader.close();
    } catch (ParseException | IOException e) {
      e.printStackTrace();
      return null;
    }
    return vMap;
  }

}
