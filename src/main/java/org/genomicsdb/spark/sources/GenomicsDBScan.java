package org.genomicsdb.spark.sources;

import org.apache.spark.sql.connector.read.Batch;
import org.apache.spark.sql.connector.read.Scan;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

import java.util.Map;

/**
 * Scan is a logical plan of a GenomicsDB Scan.
 **/
public class GenomicsDBScan implements Scan {

  private final StructType schema;
  private final Map<String, String> properties;
  private final CaseInsensitiveStringMap options;

  public GenomicsDBScan(StructType schema, Map<String, String> properties,
    CaseInsensitiveStringMap options){

    this.schema = schema;
    this.properties = properties;
    this.options = options;

  }

  @Override 
  public StructType readSchema(){
    return schema;
  }

  @Override
  public String description(){
    return "GenomicsDBSource Scan";
  }

  @Override 
  public Batch toBatch(){
    return new GenomicsDBBatch(schema, properties, options);
  }

}
