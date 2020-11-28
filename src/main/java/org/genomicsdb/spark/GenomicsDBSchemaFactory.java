package org.genomicsdb.spark;

import org.apache.spark.sql.types.StructType;
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
 
  public static StructType defaultSchema(){
    return DataTypes.createStructType(defaultFields());
  }

  public static List<StructField> extendSchema(StructField[] addFields){
    List<StructField> defaults = defaultFields();
    StructType additional = DataTypes.createStructType(addFields);
    Set<String> names = new HashSet<String>(Arrays.asList(additional.fieldNames()));
    List<StructField> fields = new ArrayList<>();
    for (StructField f: defaults){
      if (names.contains(f.name())){
        continue;
      }else{
        fields.add(f);
      }
    }
    fields.addAll(Arrays.asList(addFields));
    return fields;
  }

  private List<StructField> extendSchemaWithVid(StructField[] addFields){
    List<StructField> defaults = defaultFields();
    List<StructField> fields = new ArrayList<>();
    StructType additional = DataTypes.createStructType(fields);
    Set<String> names =  new HashSet<String>(Arrays.asList(additional.fieldNames())); 
    for (StructField f: defaults){
      if (names.contains(f.name())){
        continue;
      }else{
        fields.add(f);
      }
    }
    fields.addAll(addFieldsWithVid(addFields));
    return fields;
  }

  private List<StructField> addFieldsWithVid(StructField[] addFields){
    List<StructField> fields = new ArrayList<>();
    for (StructField f: addFields){
      GenomicsDBVidSchema field = this.vidMap.get(f.name());
      if (field == null) {
        throw new RuntimeException("Schema field " + f.name() + " not available in vid");      
      }
      fields.add(DataTypes.createStructField(f.name(), DataType.fromDDL(field.getDataType()), f.nullable()));
    }
    return fields;
  }

  public static StructType buildSchema(StructField[] addFields){
    return DataTypes.createStructType(extendSchema(addFields));
  }

  public StructType buildSchemaWithVid(StructField[] addFields){
    return DataTypes.createStructType(extendSchemaWithVid(addFields));
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
            Class<?> clazz;
            String vType = (String) v.get("type");
            switch (vType) {
              case "int":
                clazz = Integer.class;
                break;
              case "char":
                clazz = String.class;
                break;
                // TODO: switching to double below because that is what VC
                // seems to be giving us. Otherwise spark will complain later
                // when we try to create the dataframe about unboxing from double to float
              case "float":
                clazz = Double.class;
                break;
              case "flag":
                clazz = Boolean.class;
                break;
              default:
                throw new RuntimeException("Unsupported type " + vType + " in vid mapping");
            }
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
