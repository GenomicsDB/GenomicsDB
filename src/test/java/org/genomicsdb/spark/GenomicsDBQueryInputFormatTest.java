package org.genomicsdb.spark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.mapreduce.Job;
import org.genomicsdb.spark.api.GenomicsDBQueryInputFormat;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;

public class GenomicsDBQueryInputFormatTest {

  @Test
  public void testBasicQueryFormatClass() throws IOException, InterruptedException {
    GenomicsDBQueryInputFormat inputFormat = new GenomicsDBQueryInputFormat();
    Assert.assertNull(inputFormat.getConf());
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, "loader.json");
    conf.set(GenomicsDBConfiguration.QUERYJSON, "query.json");
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, "hostfile");

    inputFormat.setConf(conf);
    Assert.assertEquals(inputFormat.getConf().get(GenomicsDBConfiguration.LOADERJSON), "loader.json");
    Assert.assertEquals(inputFormat.getConf().get(GenomicsDBConfiguration.QUERYJSON), "query.json");
    Assert.assertEquals(inputFormat.getConf().get(GenomicsDBConfiguration.MPIHOSTFILE), "hostfile");

    inputFormat = new GenomicsDBQueryInputFormat(new GenomicsDBConfiguration(conf));
    Assert.assertEquals(inputFormat.getConf().get(GenomicsDBConfiguration.LOADERJSON), "loader.json");
    Assert.assertEquals(inputFormat.getConf().get(GenomicsDBConfiguration.QUERYJSON), "query.json");
    Assert.assertEquals(inputFormat.getConf().get(GenomicsDBConfiguration.MPIHOSTFILE), "hostfile");
  }

  @Test
  public void testNullBasicInputSplits() throws IOException, InterruptedException {
    final GenomicsDBQueryInputFormat tryNullInputFormat = new GenomicsDBQueryInputFormat(new GenomicsDBConfiguration());
    tryNullInputFormat.getSplits(null);

    final GenomicsDBConfiguration tryConfiguration = new GenomicsDBConfiguration();
    tryConfiguration.set(GenomicsDBConfiguration.LOADERJSON, "xxx");
    final GenomicsDBQueryInputFormat tryInputFormat = new GenomicsDBQueryInputFormat(tryConfiguration);
    tryInputFormat.getSplits(null);
  }

}
