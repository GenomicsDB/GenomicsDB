/*
 * The MIT License (MIT)
 * Copyright (c) 2023 dātma, inc™
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

import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.reader.GenomicsDBQuery;

import java.util.List;

public class TestGenomicsDBDemo {

  public static void main(String[] args) {
    String ws = System.getenv("GENOMICSDB_DEMO_WS");
    if (ws == null || ws.isEmpty()) {
      System.err.println("Env GENOMICSDB_DEMO_WS not set");
      System.exit(1);
    }
    String arrayName = "allcontigs$1$3095677412";
    Coordinates.ContigInterval interval = Coordinates.ContigInterval.newBuilder()
            .setContig("17").setBegin(7571719).setEnd(7590868).build();
    GenomicsDBExportConfiguration.RowRangeList rowRangeList = GenomicsDBExportConfiguration.RowRangeList.newBuilder().addRangeList(GenomicsDBExportConfiguration.RowRange.newBuilder()
            .setLow(0l).setHigh(200000l).build()).build();

    String filters[] = {"", "REF==\"A\"", "REF==\"A\" && ALT|=\"T\"", "REF==\"A\" && ALT|=\"T\" && GT&=\"1/1\""};

    for (String filter: filters) {
      long startTime = System.currentTimeMillis();
      GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
              .setWorkspace(ws)
              .setVidMappingFile(ws + "/vidmap.json")
              .setCallsetMappingFile(ws + "/callset.json")
              .setArrayName(arrayName)
              .setEnableSharedPosixfsOptimizations(true)
              .setBypassIntersectingIntervalsPhase(true)
              .setQueryFilter(filter)
              .addAttributes("REF").addAttributes("ALT").addAttributes("GT")
              .addQueryContigIntervals(interval)
              .addQueryRowRanges(rowRangeList)
              .build();

      GenomicsDBQuery query = new GenomicsDBQuery();
      long genomicsDBHandle = query.connectExportConfiguration(exportConfiguration);

      List<GenomicsDBQuery.Interval> intervals = query.queryVariantCalls(genomicsDBHandle);
    /* for (GenomicsDBQuery.Interval query_interval: intervals) {
      System.out.println("Interval: " + query_interval.getInterval().toString());
      List<GenomicsDBQuery.VariantCall> calls = query_interval.getCalls();
      for (GenomicsDBQuery.VariantCall call: calls) {
        System.out.println("\t call: " + call.toString());
      }
    } */
      Long elapsedTime = (System.currentTimeMillis()-startTime)/1000;
      System.out.print("Summary"+(filter.isEmpty()?"":" for ")+filter+": ");
      System.out.println("Elapsed time: " + elapsedTime + "s");
      for (GenomicsDBQuery.Interval query_interval : intervals) {
        System.out.print("Interval: " + query_interval.getInterval().toString());
        System.out.print(" calls=" + query_interval.getCalls().size());
        System.out.println();
      }
      System.out.println();
      query.disconnect(genomicsDBHandle);
    }
  }
}