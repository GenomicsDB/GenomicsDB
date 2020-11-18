package org.genomicsdb.spark;

import org.apache.spark.sql.catalyst.InternalRow;
import org.apache.spark.sql.connector.read.InputPartition;
import org.apache.spark.sql.connector.read.PartitionReader;
import org.apache.spark.sql.connector.read.PartitionReaderFactory;
import org.apache.spark.sql.types.StructType;

import java.io.FileNotFoundException;
import java.net.URISyntaxException;

public class GenomicsDBPartitionReaderFactory implements PartitionReaderFactory{
  
  public GenomicsDBPartitionReaderFactory(){}

  @Override
  public PartitionReader<InternalRow> createReader(InputPartition partition){
    return new GenomicsDBInputPartitionReader((GenomicsDBInputPartition) partition);
  }

}
