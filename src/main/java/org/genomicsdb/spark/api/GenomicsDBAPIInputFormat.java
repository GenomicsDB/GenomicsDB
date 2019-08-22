package org.genomicsdb.spark;

import org.apache.hadoop.conf.Configurable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.mapreduce.*;
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.reader.GenomicsDBQuery;
import org.genomicsdb.reader.GenomicsDBQuery.Interval;
import org.genomicsdb.reader.GenomicsDBQuery.Pair;
import org.genomicsdb.reader.GenomicsDBQuery.VariantCall;
import org.genomicsdb.spark.GenomicsDBConfiguration;
import org.genomicsdb.spark.GenomicsDBInput;
import org.genomicsdb.spark.GenomicsDBInputSplit;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class GenomicsDBAPIInputFormat extends InputFormat<Interval, List<VariantCall>> implements Configurable {

  private Configuration configuration;
  private GenomicsDBInput<GenomicsDBInputSplit> input;

  @Override
  public List<InputSplit> getSplits(JobContext jobContext) throws IOException, InterruptedException {
    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration(configuration);
    genomicsDBConfiguration.setLoaderJsonFile(
            configuration.get(GenomicsDBConfiguration.LOADERJSON));
    if (configuration.get(GenomicsDBConfiguration.QUERYPB) != null) {
      genomicsDBConfiguration.setQueryJsonFile(
              configuration.get(GenomicsDBConfiguration.QUERYPB));
    }
    else {
      genomicsDBConfiguration.setQueryJsonFile(
              configuration.get(GenomicsDBConfiguration.QUERYJSON));
    }
    if (configuration.get(GenomicsDBConfiguration.MPIHOSTFILE) != null) {
      genomicsDBConfiguration.setHostFile(
              configuration.get(GenomicsDBConfiguration.MPIHOSTFILE));
    }

    input.setGenomicsDBConfiguration(genomicsDBConfiguration);
    return (List)input.divideInput();
  }

  @Override
  public RecordReader<Interval, List<VariantCall>> createRecordReader(
      InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
      throws IOException, InterruptedException {
    String loaderJson;
    String queryJson;

    GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)inputSplit;

    boolean isPB;
    if (taskAttemptContext != null) {
      Configuration configuration = taskAttemptContext.getConfiguration();
      loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
      if (configuration.get(GenomicsDBConfiguration.QUERYPB) != null) {
        queryJson = configuration.get(GenomicsDBConfiguration.QUERYPB);
        isPB = true;
      }
      else {
        queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
        isPB = false;
      }
    } else {
      // If control comes here, means this method is called from
      // GenomicsDBRDD. Hence, the configuration object must be
      // set by setConf method, else this will lead to
      // NullPointerException
      assert(configuration!=null);
      loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
      if (configuration.get(GenomicsDBConfiguration.QUERYPB) != null) {
        queryJson = configuration.get(GenomicsDBConfiguration.QUERYPB);
        isPB = true;
      }
      else {
        queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
        isPB = false;
      }
    }

    // Need to amend query file being passed in based on inputSplit
    // so we'll create an appropriate protobuf object
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration;
    try {
      exportConfiguration =
              GenomicsDBInput.createTargetExportConfigurationPB(queryJson,
                      gSplit.getPartitionInfo(),
                      gSplit.getQueryInfoList(), isPB);
    }
    catch (ParseException e) {
      e.printStackTrace();
      return null;
    }

    return new RecordReader<Interval, List<VariantCall>>() {
      List<Interval> intervals;
      Iterator<Interval> intervalIterator;
      Interval currentInterval;
      int numProcessedIntervals = 0;

      List<Pair> ToColumnRangePairs(
          List<Coordinates.GenomicsDBColumnOrInterval> intervalLists) {
        List<Pair> intervalPairs = new ArrayList<>();
        for (Coordinates.GenomicsDBColumnOrInterval interval : intervalLists) {
          assert (interval.getColumnInterval().hasTiledbColumnInterval());
          Coordinates.TileDBColumnInterval tileDBColumnInterval =
              interval.getColumnInterval().getTiledbColumnInterval();
          assert (tileDBColumnInterval.hasBegin() && tileDBColumnInterval.hasEnd());
          intervalPairs.add(
              new Pair(tileDBColumnInterval.getBegin(), tileDBColumnInterval.getEnd()));
        }
        return intervalPairs;
      }

      List<Pair> ToRowRangePairs(
              List<GenomicsDBExportConfiguration.RowRange> rowRanges) {
        List<Pair> rangePairs = new ArrayList<>();
        for (GenomicsDBExportConfiguration.RowRange range : rowRanges) {
          assert(range.hasLow() && range.hasLow());
          rangePairs.add(
                  new Pair(range.getLow(), range.getHigh()));
        }
        return rangePairs;
      }

      @Override
      public void initialize(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
          throws IOException, InterruptedException {
        GenomicsDBQuery query = new GenomicsDBQuery();
        if (exportConfiguration.hasWorkspace()
            && exportConfiguration.hasVidMappingFile()
            && exportConfiguration.hasCallsetMappingFile()
            && exportConfiguration.hasReferenceGenome()) {
          long queryHandle = query.connect(
                  exportConfiguration.getWorkspace(),
                  exportConfiguration.getVidMappingFile(),
                  exportConfiguration.getCallsetMappingFile(),
                  exportConfiguration.getReferenceGenome(),
                  exportConfiguration.getAttributesList(),
                  exportConfiguration.getSegmentSize());
          assert(exportConfiguration.hasArrayName());
          if (exportConfiguration.getQueryRowRangesCount() > 0) {
            if (exportConfiguration.getQueryColumnRangesCount() > 0) {
              assert(exportConfiguration.getQueryColumnRangesCount() == 0);
              intervals = query.queryVariantCalls(
                      queryHandle,
                      exportConfiguration.getArrayName(),
                      ToColumnRangePairs(exportConfiguration.getQueryColumnRanges(0).getColumnOrIntervalListList()));
            }
            assert(exportConfiguration.getQueryRowRangesCount() == 0);
            intervals = query.queryVariantCalls(
                    queryHandle,
                    exportConfiguration.getArrayName(),
                    ToColumnRangePairs(exportConfiguration.getQueryColumnRanges(0).getColumnOrIntervalListList()),
                    ToRowRangePairs(exportConfiguration.getQueryRowRanges(0).getRangeListList()));
          } else {
            intervals = query.queryVariantCalls(queryHandle, exportConfiguration.getArrayName());
          }
          intervalIterator = intervals.iterator();
        }
      }

      @Override
      public boolean nextKeyValue() throws IOException, InterruptedException {
        if (intervalIterator.hasNext()) {
          currentInterval = intervalIterator.next();
          numProcessedIntervals++;
          return true;
        } else {
          currentInterval = null;
          return false;
        }
      }

      @Override
      public Interval getCurrentKey() throws IOException, InterruptedException {
        return currentInterval;
      }

      @Override
      public List<VariantCall> getCurrentValue() throws IOException, InterruptedException {
        return getCurrentKey().getCalls();
      }

      @Override
      public float getProgress() throws IOException, InterruptedException {
        return numProcessedIntervals / intervals.size();
      }

      @Override
      public void close() throws IOException {}
    };
  }

  /** default constructor */
  public GenomicsDBAPIInputFormat() {
    input = new GenomicsDBInput<>(null, null, null, 1, Long.MAX_VALUE, GenomicsDBInputSplit.class);
  }

  public GenomicsDBAPIInputFormat(GenomicsDBConfiguration conf) {
    input = new GenomicsDBInput<>(conf, null, null, 1, Long.MAX_VALUE, GenomicsDBInputSplit.class);
  }

  @Override
  public void setConf(Configuration configuration) {
    this.configuration = configuration;
  }

  @Override
  public Configuration getConf() {
    return configuration;
  }
}
