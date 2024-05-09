/**
 * src/test/cpp/src/test_genomicsdb_demo.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2024 dātma, inc™
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
 *
 * Test the GenomicsDB query api with the genomicsdb demo workspace
 */

#include "genomicsdb.h"
#include "genomicsdb_config_base.h"
#include "genomicsdb_logger.h"
#include "memory_measure.h"
#include "tiledb_utils.h"
#ifdef __linux__
#  include <malloc.h>
#endif

#ifdef USE_GPERFTOOLS_HEAP
#include "gperftools/heap-profiler.h"
#endif

#include "genomicsdb_export_config.pb.h"

/**
 * Instructions to run this test...
 * Define GENOMICSDB_DEMO_WS(test does not run otherwise) and NUM_ITERATIONS(default=1) env variables
 * GENOMICSDB_DEMO_WS=<genomicsdb_demo_ws_path> NUM_ITERATIONS=10 ./test_genomicsdb_demo
 * with valgrind...
 * GENOMICSDB_DEMO_WS=<genomicsdb_demo_ws_path> valgrind --leak-check=full --suppressions=<GenomicsDB>/src/test/inputs/valgrind.supp ./test_genomicsdb_demo
 * To attach a degugger with valgrind, invoke valgrind --vgdb=yes --vgdb-error=0 --suppressions=<GenomiscDB>/src/test/inputs/valgrind.supp ./test_genomicsdb_demo
 */
class CountCellsProcessor : public GenomicsDBVariantCallProcessor {
 public:
  CountCellsProcessor() {
  };

  void process(const interval_t& interval) {
    m_intervals++;
  };

  void process(const std::string& sample_name,
               const int64_t* coordinates,
               const genomic_interval_t& genomic_interval,
               const std::vector<genomic_field_t>& genomic_fields) {
    m_count++;
    m_coordinates[0] = coordinates[0];
    m_coordinates[1] = coordinates[1];
  };

  int m_intervals = 0;
  int m_count = 0;
  int64_t m_coordinates[2];
};

int main(int argc, char** argv) {
  using namespace genomicsdb_pb;

  char *genomicsdb_demo_workspace = getenv("GENOMICSDB_DEMO_WS");
  if (!genomicsdb_demo_workspace) return 0;

  unsigned num_iterations = 1;
  char *num_iterations_str = getenv("NUM_ITERATIONS");
  if (num_iterations_str) num_iterations = std::stoul(num_iterations_str);

  printf("num_iterations=%u\n", num_iterations);

  bool use_single_handle = false;
  char *use_single_handle_str = getenv("USE_SINGLE_HANDLE");
  if (use_single_handle_str && strncmp(use_single_handle_str, "true", 4) == 0) {
    use_single_handle = true;
  }

  ExportConfiguration config;

  std::string ws(genomicsdb_demo_workspace);
  config.set_workspace(ws);
  config.set_array_name("allcontigs$1$3095677412");
  config.set_callset_mapping_file(ws+"/callset.json");
  config.set_vid_mapping_file(ws+"/vidmap.json");

  // query_contig_intervals
  auto* contig_interval = config.add_query_contig_intervals();
  contig_interval->set_contig("17");
  contig_interval->set_begin(7571719);
  contig_interval->set_end(7590868);

  // query_row_ranges
  auto* row_range = config.add_query_row_ranges()->add_range_list();
  row_range->set_low(0);
  row_range->set_high(200000);

  // query_attributes
  config.add_attributes()->assign("REF");
  config.add_attributes()->assign("ALT");
  config.add_attributes()->assign("GT");

  // other
  config.set_bypass_intersecting_intervals_phase(true);
  config.set_enable_shared_posixfs_optimizations(true);

  std::string config_string;
  config.SerializeToString(&config_string);

  print_memory_usage("before loop");

#ifdef USE_GPERFTOOLS_HEAP
    HeapProfilerStart("test_genomicsdb_demo.gperf.heap");
#endif

  if (use_single_handle) {
    GenomicsDB gdb(config_string, GenomicsDB::PROTOBUF_BINARY_STRING);
    for (auto i=0u; i<num_iterations; i++) {
      CountCellsProcessor count_cells_processor;
      print_memory_usage("before query");
#ifdef USE_GPERFTOOLS_HEAP
      HeapProfilerDump("before_query");
#endif
      gdb.query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
#ifdef USE_GPERFTOOLS_HEAP
      HeapProfilerDump("after_query_before_trim");
#endif
      print_memory_usage("after query before trim");
#ifdef __linux__
//      malloc_trim(0);
#endif
#ifdef USE_GPERFTOOLS_HEAP
      HeapProfilerDump("after_query");
#endif
      print_memory_usage("after query");
      printf("count=%d\n", count_cells_processor.m_count);
    }
  } else {
    for (auto i=0u; i<num_iterations; i++) {
      print_memory_usage("before gdb and query");
      GenomicsDB *gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING);
      CountCellsProcessor count_cells_processor;
      print_memory_usage("before query");
      gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
      delete gdb;
      print_memory_usage("after query");
      printf("count=%d\n", count_cells_processor.m_count);
   }
  }

#ifdef USE_GPERFTOOLS_HEAP
    HeapProfilerStop();
#endif
}
