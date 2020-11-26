package org.genomicsdb.spark;

import org.apache.spark.sql.connector.catalog.Table;
import org.apache.spark.sql.connector.catalog.TableProvider;
import org.apache.spark.sql.connector.expressions.Transform;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.types.*;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

import java.util.Map;

/** 
 * This replaces the primary data source interface in the DatasourceV2 API, 
 * which was removed as of Spark 3.0.x and replaced with TableProvider interface. 
 * This is the base interface for all custom datasource. 
 * 
 **/

public class GenomicsDBSource implements TableProvider {

  public GenomicsDBSource(){}

  /** 
   * This sets the default schema, which will contain all the default fields coming 
   * from the backend. This was previously in the setSchemaOptions in the GenomicsDBBatch.
   *
   **/
  public StructType inferSchema(CaseInsensitiveStringMap options){
    return GenomicsDBSchemaFactory.defaultSchema();
  }

  @Override
  public Table getTable(StructType schema, Transform[] partitioning, Map<String, String> properties){
    return new GenomicsDBTable(schema, properties);
  }

  @Override
  public boolean supportsExternalMetadata(){
    return true;
  }

}
