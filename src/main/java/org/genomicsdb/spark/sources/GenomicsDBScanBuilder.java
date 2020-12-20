package org.genomicsdb.spark.sources;

import org.apache.spark.sql.connector.read.Scan;
import org.apache.spark.sql.connector.read.ScanBuilder;
import org.apache.spark.sql.types.StructType;
import org.apache.spark.sql.util.CaseInsensitiveStringMap;

import java.util.Map;

/**
 * Returns a Scan instance to create a logical plan.
 **/
public class GenomicsDBScanBuilder implements ScanBuilder {

  private final StructType schema;
  private final Map<String, String> properties;
  private final CaseInsensitiveStringMap options;

  public GenomicsDBScanBuilder(StructType schema, Map<String, String> properties,
    CaseInsensitiveStringMap options) {

    this.schema = schema;
    this.properties = properties;
    this.options = options;
  }

  @Override
  public Scan build() {
    return new GenomicsDBScan(schema,properties,options);
  }


}
