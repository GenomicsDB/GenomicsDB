package org.genomicsdb.spark.api;

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
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class GenomicsDBQueryInputFormat extends InputFormat<Interval, List<VariantCall>> implements Configurable {

  private Configuration configuration;
  private GenomicsDBInput<GenomicsDBInputSplit> input;

  @Override
  public List<InputSplit> getSplits(JobContext jobContext) throws IOException, InterruptedException {
    GenomicsDBConfiguration genomicsDBConfiguration = new GenomicsDBConfiguration(configuration);

    // At least a query in the form of json or protobuf should be passed
    if (configuration.get(GenomicsDBConfiguration.QUERYPB) == null && configuration.get(GenomicsDBConfiguration.QUERYJSON) == null) {
      throw new IOException("Query json or query protobuf has to be specified.");
    }

    if (configuration.get(GenomicsDBConfiguration.LOADERJSON) != null) {
      genomicsDBConfiguration.setLoaderJsonFile(
          configuration.get(GenomicsDBConfiguration.LOADERJSON));
    }
    if (configuration.get(GenomicsDBConfiguration.QUERYPB) != null) {
      genomicsDBConfiguration.setQueryPB(
              configuration.get(GenomicsDBConfiguration.QUERYPB));
    } else if (configuration.get(GenomicsDBConfiguration.QUERYJSON) != null){
      genomicsDBConfiguration.setQueryJsonFile(
              configuration.get(GenomicsDBConfiguration.QUERYJSON));
    }
    if (configuration.get(GenomicsDBConfiguration.MPIHOSTFILE) != null) {
      genomicsDBConfiguration.setHostFile(
              configuration.get(GenomicsDBConfiguration.MPIHOSTFILE));
    }
    setConf(genomicsDBConfiguration);
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
      configuration = taskAttemptContext.getConfiguration();
    } else {
      assert(configuration!=null);
    }

    loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
    if (configuration.get(GenomicsDBConfiguration.QUERYPB) != null) {
      queryJson = configuration.get(GenomicsDBConfiguration.QUERYPB);
      isPB = true;
    } else {
      queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
      isPB = false;
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
      boolean initialized = false;
      List<Interval> intervals;
      Iterator<Interval> intervalIterator;
      Interval currentInterval;
      int numProcessedIntervals = 0;

      GenomicsDBQuery query = new GenomicsDBQuery();

      List<Pair> ToColumnRangePairs(List<Coordinates.GenomicsDBColumnOrInterval> intervalLists) {
        List<Pair> intervalPairs = new ArrayList<>();
        for (Coordinates.GenomicsDBColumnOrInterval interval : intervalLists) {
          assert (interval.getColumnInterval().hasTiledbColumnInterval());
          Coordinates.TileDBColumnInterval tileDBColumnInterval =
              interval.getColumnInterval().getTiledbColumnInterval();
          interval.getColumnInterval().getTiledbColumnInterval();
          assert (tileDBColumnInterval.hasBegin() && tileDBColumnInterval.hasEnd());
          intervalPairs.add(
              new Pair(tileDBColumnInterval.getBegin(), tileDBColumnInterval.getEnd()));
        }
        return intervalPairs;
      }

      List<Pair> ToRowRangePairs(List<GenomicsDBExportConfiguration.RowRange> rowRanges) {
        List<Pair> rangePairs = new ArrayList<>();
        for (GenomicsDBExportConfiguration.RowRange range : rowRanges) {
          assert (range.hasLow() && range.hasLow());
          rangePairs.add(new Pair(range.getLow(), range.getHigh()));
        }
        return rangePairs;
      }

      private boolean check_configuration(String key, String value) {
        if (value == null || value.isEmpty()) {
          System.err.println("GenomicsDB Configuration does not contain value for key=" + key);
          return false;
        } else {
          return true;
        }
      }

      @Override
      public void initialize(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
              throws IOException, InterruptedException {
        String workspace = null;
        String vidMappingFile = null;
        String callsetMappingFile = null;
        String referenceGenome = null;
        Long segmentSize = 0L;

        JSONObject jsonObject = null;
        try {
          JSONParser parser = new JSONParser();
          jsonObject = (JSONObject)parser.parse(new FileReader(configuration.get(GenomicsDBConfiguration.LOADERJSON)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
          e.printStackTrace();
        } catch (ParseException e) {
          e.printStackTrace();
        }

        assert (exportConfiguration.hasArrayName());

        long queryHandle;
        if (jsonObject != null) {
          workspace = (exportConfiguration.hasWorkspace())?exportConfiguration.getWorkspace():(String)jsonObject.get("workspace");
          vidMappingFile = (exportConfiguration.hasVidMappingFile())?exportConfiguration.getVidMappingFile():(String)jsonObject.get("vid_mapping_file");
          callsetMappingFile = (exportConfiguration.hasCallsetMappingFile())?exportConfiguration.getCallsetMappingFile():(String)jsonObject.get("callset_mapping_file");
          referenceGenome = (exportConfiguration.hasReferenceGenome())?exportConfiguration.getReferenceGenome():(String)jsonObject.get("reference_genome");
          segmentSize = (exportConfiguration.hasSegmentSize())?exportConfiguration.getSegmentSize():(Long)jsonObject.get("segment_size");
          if (!check_configuration("workspace", workspace) ||
                  !check_configuration("vid_mapping_file", vidMappingFile) ||
                  !check_configuration("callset_mapping_file", callsetMappingFile) ||
                  !check_configuration("reference_genome", referenceGenome)) {
            throw new RuntimeException("GenomicsDBConfiguration is incomplete. Add required configuration values and restart the operation");
          }
          List<String> attributesList = exportConfiguration.getAttributesList();
          if (segmentSize > 0) {
            queryHandle = query.connect(workspace, vidMappingFile, callsetMappingFile, referenceGenome, attributesList, segmentSize.longValue());
          } else {
            queryHandle = query.connect(workspace, vidMappingFile, callsetMappingFile, referenceGenome, attributesList);
          }
          intervals = query.queryVariantCalls(queryHandle, exportConfiguration.getArrayName(),
                  ToColumnRangePairs(exportConfiguration.getQueryColumnRanges(0).getColumnOrIntervalListList()),
                  ToRowRangePairs(exportConfiguration.getQueryRowRanges(0).getRangeListList()));
        } else {
          queryHandle = query.connectExportConfiguration(exportConfiguration);
          intervals = query.queryVariantCalls(queryHandle, exportConfiguration.getArrayName());
        }
        query.disconnect(queryHandle);

        intervalIterator = intervals.iterator();
        initialized = true;
      }

      @Override
      public boolean nextKeyValue() throws IOException, InterruptedException {
        if (initialized && intervalIterator.hasNext()) {
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
  public GenomicsDBQueryInputFormat() {
    input = new GenomicsDBInput<>(null, null, null, 1, Long.MAX_VALUE, GenomicsDBInputSplit.class);
  }

  public GenomicsDBQueryInputFormat(GenomicsDBConfiguration conf) {
    this.configuration = conf;
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
