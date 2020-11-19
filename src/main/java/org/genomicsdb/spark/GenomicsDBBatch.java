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

import org.apache.spark.sql.catalyst.InternalRow;
//import org.apache.spark.sql.sources.v2.DataSourceOptions;
//import org.apache.spark.sql.sources.v2.reader.DataSourceReader;
//import org.apache.spark.sql.sources.v2.reader.InputPartition;
import org.apache.spark.sql.connector.read.Batch;
import org.apache.spark.sql.connector.read.InputPartition;
import org.apache.spark.sql.connector.read.PartitionReaderFactory;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

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
import java.util.List;
import java.util.Map;

//public class GenomicsDBDataSourceReader implements DataSourceReader {

public class GenomicsDBBatch implements Batch {

  private GenomicsDBInput<GenomicsDBInputPartition> input;
  private final StructType schema;
  private final Map<String, String> properties;
  private final CaseInsensitiveStringMap options;

  //public GenomicsDBDataSourceReader() {}

  public GenomicsDBBatch(StructType schema, Map<String, String> properties, 
      CaseInsensitiveStringMap options) {
    this.schema = schema;
    this.properties = properties;
    this.options = options;
    this.input = setSchemaOptions(options, schema);
  }


  /**
   * TODO documentation
   * @param options todo 
   * @param schema todo
   * @return GenomcsDBInput
   **/
  private GenomicsDBInput<GenomicsDBInputPartition> setSchemaOptions(CaseInsensitiveStringMap options, StructType schema)
      throws RuntimeException {
    StructField[] fields = new StructField[] {
          new StructField("contig", DataTypes.StringType, false, Metadata.empty()),
          new StructField("startPos", DataTypes.IntegerType, false, Metadata.empty()),
          new StructField("ID", DataTypes.StringType, true, Metadata.empty()),
          new StructField("variantType", DataTypes.StringType, true, Metadata.empty()),
          new StructField("refAllele", DataTypes.StringType, true, Metadata.empty()),
          new StructField(
              "alternateAlleles",
              new ArrayType(DataTypes.StringType, true),
              true,
              Metadata.empty()),
          new StructField(
              "sampleNames", new ArrayType(DataTypes.StringType, false), false, Metadata.empty()),
          new StructField("GT", new ArrayType(DataTypes.StringType, true), true, Metadata.empty())
        };
    StructType defaultSchema = new StructType(fields);
    // comment out below for now...revisit if needed
    //                             StructField("encodedGenotypes", VectorType, false))

    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration();
    if(options.containsKey(GenomicsDBConfiguration.LOADERJSON)) {
      genomicsDBConfiguration.setLoaderJsonFile(options.get(GenomicsDBConfiguration.LOADERJSON));
    }
    else {
      throw new RuntimeException("Must specify "+GenomicsDBConfiguration.LOADERJSON);
    }
    if(options.containsKey(GenomicsDBConfiguration.QUERYPB)) {
      genomicsDBConfiguration.setQueryPB(options.get(GenomicsDBConfiguration.QUERYPB));
    }
    else if(options.containsKey(GenomicsDBConfiguration.QUERYJSON)) {
      genomicsDBConfiguration.setQueryJsonFile(options.get(GenomicsDBConfiguration.QUERYJSON));
    }
    else {
      throw new RuntimeException("Must specify either "+GenomicsDBConfiguration.QUERYJSON+
              " or "+GenomicsDBConfiguration.QUERYPB);
    }
    if(options.containsKey(GenomicsDBConfiguration.MPIHOSTFILE)) {
      try {
        genomicsDBConfiguration.setHostFile(options.get(GenomicsDBConfiguration.MPIHOSTFILE));
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
    }
    Map<String, GenomicsDBVidSchema> vMap = buildVidSchemaMap(options.get(GenomicsDBConfiguration.LOADERJSON));

    StructType finalSchema = defaultSchema;
    if (schema != null) {
      for (StructField sf : JavaConverters.asJavaIterableConverter(schema.toIterable()).asJava()) {
        GenomicsDBVidSchema field = vMap.get(sf.name());
        if (field == null) {
          throw new RuntimeException("Schema field " + sf.name() + " not available in vid");
        }
        String typeName;
        // could use some sort of reflection to go from java class to scala
        // type I guess...but we only care about limited types so this
        // should do fine
        if (field.getFieldClass().equals(String.class)) {
          typeName = "String";
        } else if (field.getFieldClass().equals(Integer.class)) {
          typeName = "Integer";
        }
        // TODO: switching to double below because that is what VC
        // seems to be giving us. Otherwise spark will complain later
        // when we try to create the dataframe about unboxing from double to float
        else if (field.getFieldClass().equals(Double.class)) {
          typeName = "Double";
        } else if (field.getFieldClass().equals(Boolean.class)) {
          typeName = "Boolean";
        } else {
          throw new RuntimeException("Unsupported field type in vid!");
        }
        String dtypeDDL;
        // single item INFO fields
        if (field.isInfo()
            && (field.getLength().equals("1") || field.getFieldClass().equals(String.class))) {
          dtypeDDL = typeName;
        }
        // multi-item INFO fields OR
        // single item FORMAT fields
        else if (field.isInfo()
            || field.getLength().equals("1")
            || field.getFieldClass().equals(String.class)) {
          dtypeDDL = "ARRAY<" + typeName + ">";
        }
        // multi-item FORMAT fields
        else {
          dtypeDDL = "ARRAY<ARRAY<" + typeName + ">>";
        }

        finalSchema = finalSchema.add(
                new StructField(
                    sf.name(), DataType.fromDDL(dtypeDDL), sf.nullable(), sf.metadata()));
      }
    }
    Long blocksize = new Long(1);
    Long maxblock = Long.MAX_VALUE;
    if (options.get("genomicsdb.minqueryblocksize") != null){
      blocksize = Long.valueOf(options.get("genomicsdb.minqueryblocksize"));
    }
    if (options.get("genomicsdb.maxqueryblocksize") != null){
      maxblock = Long.valueOf(options.get("genomicsdb.maxqueryblocksize"));
    }
    input = new GenomicsDBInput<>(
            genomicsDBConfiguration,
            finalSchema,
            vMap,
            blocksize,
            maxblock,
            GenomicsDBInputPartition.class);
    return input;
  }

  private Map<String, GenomicsDBVidSchema> buildVidSchemaMap(String loader) {
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

  @Override
  @SuppressWarnings("unchecked")
  public InputPartition[] planInputPartitions(){
    InputPartition[] ipartitions = new InputPartition[0];
    return input.divideInput().toArray(ipartitions);
  }
  
  @Override
  public PartitionReaderFactory createReaderFactory(){
    return new GenomicsDBPartitionReaderFactory();
  }  

}
