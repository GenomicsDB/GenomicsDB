/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 * Copyright (c) 2019-2022 Omics Data Automation, Inc.
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
 **/

#include <iostream>
#include <string>
#include <getopt.h>
#include <mpi.h>
#include "json_config.h"
#include "timer.h"
#include "query_variants.h"
#include "broad_combined_gvcf.h"
#include "vid_mapper_pb.h"
#include "genomicsdb.h"
#include "tiledb.h"
#include "tiledb_utils.h"

#ifdef USE_BIGMPI
#include "bigmpi.h"
#endif

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

#ifdef USE_GPERFTOOLS_HEAP
#include "gperftools/heap-profiler.h"
#endif

//#define VERBOSE 1

enum ArgsEnum {
  ARGS_IDX_SKIP_QUERY_ON_ROOT=1000,
  ARGS_IDX_PRODUCE_BROAD_GVCF,
  ARGS_IDX_PRODUCE_HISTOGRAM,
  ARGS_IDX_PRINT_CALLS,
  ARGS_IDX_PRINT_CSV,
  ARGS_IDX_BGEN,
  ARGS_IDX_BGEN_COMPRESSION,
  ARGS_IDX_BED,
  ARGS_IDX_TPED,
  ARGS_IDX_PLINK,
  ARGS_IDX_PLINK_PREFIX,
  ARGS_IDX_PLINK_ONE_PASS,
  ARGS_IDX_FAM_LIST,
  ARGS_IDX_PROG,
  ARGS_IDX_VERSION,
  ARGS_IDX_PRODUCE_INTERESTING_POSITIONS,
  ARGS_IDX_PRINT_ALT_ALLELE_COUNTS,
  ARGS_IDX_COLUMNAR_GVCF,
  ARGS_IDX_GVCF_PROFILE
};

enum CommandsEnum {
  COMMAND_RANGE_QUERY=0,
  COMMAND_PRODUCE_BROAD_GVCF,
  COMMAND_PRODUCE_HISTOGRAM,
  COMMAND_PRINT_CALLS,
  COMMAND_PRINT_CSV,
  COMMAND_PLINK,
  COMMAND_PRINT_ALT_ALLELE_COUNTS,
  COMMAND_COLUMNAR_GVCF
};

enum ProduceBroadGVCFSubOperation {
  PRODUCE_BROAD_GVCF_PRODUCE_GVCF=0,
  PRODUCE_BROAD_GVCF_PRODUCE_INTERESTING_POSITIONS,
  PRODUCE_BROAD_GVCF_COUNT_LINES,
  PRODUCE_BROAD_GVCF_UNKNOWN
};

#define MegaByte (1024*1024)

#ifdef DO_PROFILING
enum TimerTypesEnum {
  TIMER_TILEDB_QUERY_RANGE_IDX=0u,
  TIMER_BINARY_SERIALIZATION_IDX,
  TIMER_MPI_GATHER_IDX,
  TIMER_ROOT_BINARY_DESERIALIZATION_IDX,
  TIMER_JSON_PRINTING_IDX,
  TIMER_NUM_TIMERS
};

auto g_timer_names = std::vector<std::string> {
  "TileDB-query",
  "Binary-serialization",
  "MPI-gather",
  "Binary-deserialization",
  "JSON-printing"
};
#endif

#ifdef NDEBUG
#define ASSERT(X) if(!(X)) { std::cerr << "Assertion failed - exiting\n"; exit(-1); }
#else
#define ASSERT(X) assert(X)
#endif

//id_mapper could be NULL - use for contig/callset name mapping only if non-NULL
void run_range_query(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config, const VidMapper& id_mapper,
                     const std::string& output_format, const bool is_partitioned_by_column, int num_mpi_processes, int my_world_mpi_rank, bool skip_query_on_root) {
  //Check if id_mapper is initialized before using it
  //if(id_mapper.is_initialized())
  GTProfileStats* stats_ptr = 0;
#ifdef DO_PROFILING
  //Performance measurements
  Timer timer;
  std::vector<double> timings(2u*TIMER_NUM_TIMERS, 0);  //[cpu-time,wall-clock time]
  ASSERT(g_timer_names.size() == TIMER_NUM_TIMERS);
  GTProfileStats stats;
  stats_ptr = &stats;
#endif
  //Variants vector
  std::vector<Variant> variants;
  uint64_t num_column_intervals = query_config.get_num_column_intervals();
  std::vector<uint64_t> queried_column_positions(num_column_intervals * 2, 0ull);
  std::vector<uint64_t> query_column_lengths(num_column_intervals, 0ull);
  //Perform query if not root or !skip_query_on_root
  if (my_world_mpi_rank != 0 || !skip_query_on_root) {
#ifdef DO_PROFILING
    timer.start();
#endif
    for (auto i=0u; i<query_config.get_num_column_intervals(); ++i) {
      qp.gt_get_column_interval(qp.get_array_descriptor(), query_config, i, variants, 0, stats_ptr);
      query_column_lengths[i] = variants.size();
      queried_column_positions[i * 2] = query_config.get_column_begin(i);
      queried_column_positions[i * 2 + 1] = query_config.get_column_end(i);
    }

#ifdef DO_PROFILING
    timer.stop();
    timer.get_last_interval_times(timings, TIMER_TILEDB_QUERY_RANGE_IDX);
#endif
  }
#ifdef DO_PROFILING
  timer.start();
#endif
#if VERBOSE>0
  std::cerr << "[Rank "<< my_world_mpi_rank << " ]: Completed query, obtained "<<variants.size()<<" variants\n";
#endif
  //serialized variant data
  std::vector<uint8_t> serialized_buffer;
  serialized_buffer.resize(1000000u);       //1MB, arbitrary value - will be resized if necessary by serialization functions
  uint64_t serialized_length = 0ull;
  for (const auto& variant : variants)
    variant.binary_serialize(serialized_buffer, serialized_length);
#if VERBOSE>0
  std::cerr << "[Rank "<< my_world_mpi_rank << " ]: Completed serialization, serialized data size "
            << std::fixed << std::setprecision(3) << ((double)serialized_length)/MegaByte  << " MBs\n";
#endif
#ifdef DO_PROFILING
  timer.stop();
  timer.get_last_interval_times(timings, TIMER_BINARY_SERIALIZATION_IDX);
  //Gather profiling data at root
  auto num_timing_values_per_mpi_process = 2*(TIMER_BINARY_SERIALIZATION_IDX+1);
  std::vector<double> gathered_timings(num_mpi_processes*num_timing_values_per_mpi_process);
  ASSERT(MPI_Gather(&(timings[0]), num_timing_values_per_mpi_process, MPI_DOUBLE,
                    &(gathered_timings[0]), num_timing_values_per_mpi_process, MPI_DOUBLE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  timer.start();
#endif
  //Gather all serialized lengths (in bytes) at root
  std::vector<uint64_t> lengths_vector(num_mpi_processes, 0ull);
  ASSERT(MPI_Gather(&serialized_length, 1, MPI_UNSIGNED_LONG_LONG, &(lengths_vector[0]), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  //Total size of gathered data (valid at root only)
  auto total_serialized_size = 0ull;
  //Buffer to receive all gathered data, will be resized at root
  std::vector<uint8_t> receive_buffer(1u);
#ifdef USE_BIGMPI
  std::vector<MPI_Count> recvcounts(num_mpi_processes);
  std::vector<MPI_Aint> displs(num_mpi_processes);
#else
  //MPI uses int for counts and displacements
  std::vector<int> recvcounts(num_mpi_processes);
  std::vector<int> displs(num_mpi_processes);
#endif
  //root
  if (my_world_mpi_rank == 0) {
    for (auto val : lengths_vector)
      total_serialized_size += val;
#if VERBOSE>0
    std::cerr << "[Rank "<< my_world_mpi_rank << " ]: Gathered lengths, total size "
              << std::fixed << std::setprecision(3) << (((double)total_serialized_size)/MegaByte)<<" MBs\n";
#endif
#ifndef USE_BIGMPI
    if (total_serialized_size >= static_cast<uint64_t>(INT_MAX)) { //max 32 bit signed int
      std::cerr << "Serialized size beyond 32-bit int limit - exiting.\n";
      std::cerr <<  "Use the BigMPI library: https://github.com/jeffhammond/BigMPI and recompile the TileDB library with USE_BIGMPI=<path>\n";
      exit(-1);
    }
#endif
    receive_buffer.resize(total_serialized_size);
    auto curr_displ = 0ull;
    //Fill in recvcounts and displs vectors
    for (auto i=0u; i<recvcounts.size(); ++i) {
      auto curr_length = lengths_vector[i];
      recvcounts[i] = curr_length;
      displs[i] = curr_displ;
      curr_displ += curr_length;
    }
  }
  //Gather serialized variant data
#if VERBOSE>0
  if (my_world_mpi_rank == 0)
    std::cerr << "Starting MPIGather at root into buffer of size "<<((double)receive_buffer.size())/MegaByte<<"\n";
#endif
#ifdef USE_BIGMPI
  ASSERT(MPIX_Gatherv_x(&(serialized_buffer[0]), serialized_length, MPI_UNSIGNED_CHAR, &(receive_buffer[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
#else
  ASSERT(MPI_Gatherv(&(serialized_buffer[0]), serialized_length, MPI_UNSIGNED_CHAR, &(receive_buffer[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif
  std::vector<uint64_t> gathered_num_column_intervals(num_mpi_processes);
  ASSERT(MPI_Gather(&num_column_intervals, 1, MPI_UNSIGNED_LONG_LONG,
                    &(gathered_num_column_intervals[0]), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD) == MPI_SUCCESS);

  uint64_t total_columns = 0ull;
  if (my_world_mpi_rank == 0) {
    for (auto i = 0ull; i < recvcounts.size(); ++i) {
      // Reuse the displs and recvcounts vectors declared for variants serialization
      displs[i] = total_columns;
      total_columns += gathered_num_column_intervals[i];
      recvcounts[i] = gathered_num_column_intervals[i];
    }
  }

  std::vector<uint64_t> gathered_query_column_lengths(total_columns);
#ifdef USE_BIGMPI
  ASSERT(MPIX_Gatherv_x(&(query_column_lengths[0]), num_column_intervals, MPI_UNSIGNED_LONG_LONG,
                        &(gathered_query_column_lengths[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
                        0, MPI_COMM_WORLD) == MPI_SUCCESS);
#else
  ASSERT(MPI_Gatherv(&(query_column_lengths[0]), num_column_intervals, MPI_UNSIGNED_LONG_LONG,
                     &(gathered_query_column_lengths[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
                     0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif

  // Update the recvcounts and displs but multiply by 2 because
  // each column has the start and end query position
  total_columns = 0ull;
  if (my_world_mpi_rank == 0) {
    for (auto i = 0ull; i < recvcounts.size(); ++i) {
      // Reuse the displs and recvcounts vectors declared for variants serialization
      displs[i] = total_columns;
      total_columns += gathered_num_column_intervals[i] * 2;
      recvcounts[i] = gathered_num_column_intervals[i] * 2;
    }
  }

  std::vector<uint64_t> gathered_queried_column_positions(total_columns);
#ifdef USE_BIGMPI
  ASSERT(MPIX_Gatherv_x(&(queried_column_positions[0]), num_column_intervals * 2, MPI_UNSIGNED_LONG_LONG,
                        &(gathered_queried_column_positions[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
                        0, MPI_COMM_WORLD) == MPI_SUCCESS);
#else
  ASSERT(MPI_Gatherv(&(queried_column_positions[0]), num_column_intervals * 2, MPI_UNSIGNED_LONG_LONG,
                     &(gathered_queried_column_positions[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
                     0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif

#if VERBOSE>0
  if (my_world_mpi_rank == 0)
    std::cerr << "Completed MPI_Gather\n";
#endif
#ifdef DO_PROFILING
  timer.stop();
  timer.get_last_interval_times(timings, TIMER_MPI_GATHER_IDX);
  timer.start();
#endif
  //Deserialize at root
  if (my_world_mpi_rank == 0) {
    variants.clear();
    uint64_t offset = 0ull;
    while (offset < total_serialized_size) {
      variants.emplace_back();
      auto& variant = variants.back();
      qp.binary_deserialize(variant, query_config, receive_buffer, offset);
    }
#if VERBOSE>0
    std::cerr << "Completed binary deserialization at root\n";
#endif
#ifdef DO_PROFILING
    timer.stop();
    timer.get_last_interval_times(timings, TIMER_ROOT_BINARY_DESERIALIZATION_IDX);
    timer.start();
#endif
    print_variants(variants, output_format, query_config, std::cout, is_partitioned_by_column, &id_mapper,
                   gathered_query_column_lengths, gathered_num_column_intervals, gathered_queried_column_positions);
#ifdef DO_PROFILING
    timer.stop();
    timer.get_last_interval_times(timings, TIMER_JSON_PRINTING_IDX);
    std::cerr << "Root received "<< std::fixed << std::setprecision(3) << (((double)total_serialized_size)/MegaByte)
              << " MBs of variant data in binary format\n";
    for (auto i=0u; i<TIMER_NUM_TIMERS; ++i) {
      std::cerr << g_timer_names[i];
      if (i >= TIMER_MPI_GATHER_IDX) //only root info
        std::cerr << std::fixed << std::setprecision(3) << "," << timings[2*i] << ",," << timings[2*i+1u] << "\n";
      else {
        assert(2*i+1 < static_cast<unsigned>(num_timing_values_per_mpi_process));
        for (auto j=0u; j<gathered_timings.size(); j+=num_timing_values_per_mpi_process)
          std::cerr << std::fixed << std::setprecision(3) << "," << gathered_timings[j + 2*i];
        std::cerr << ",";
        for (auto j=0u; j<gathered_timings.size(); j+=num_timing_values_per_mpi_process)
          std::cerr << std::fixed << std::setprecision(3) << "," << gathered_timings[j + 2*i + 1];
        std::cerr << "\n";
      }
    }
#endif
  }
}

void scan_and_produce_Broad_GVCF(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config,
                                 VCFAdapter& vcf_adapter, const VidMapper& id_mapper,
                                 const ProduceBroadGVCFSubOperation sub_operation_type, int my_world_mpi_rank, bool skip_query_on_root) {
  //Read output in batches if required
  //Must initialize buffer before constructing gvcf_op
  RWBuffer rw_buffer;
  auto serialized_vcf_adapter_ptr = dynamic_cast<VCFSerializedBufferAdapter*>(&vcf_adapter);
  if (serialized_vcf_adapter_ptr)
    serialized_vcf_adapter_ptr->set_buffer(rw_buffer);
  SingleVariantOperatorBase* op_ptr = 0;
  switch (sub_operation_type) {
  case ProduceBroadGVCFSubOperation::PRODUCE_BROAD_GVCF_PRODUCE_GVCF:
    op_ptr = new BroadCombinedGVCFOperator(vcf_adapter, id_mapper, query_config);
    break;
  case ProduceBroadGVCFSubOperation::PRODUCE_BROAD_GVCF_PRODUCE_INTERESTING_POSITIONS:
    op_ptr = new InterestingLocationsPrinter(std::cout, query_config);
    break;
  case ProduceBroadGVCFSubOperation::PRODUCE_BROAD_GVCF_COUNT_LINES:
    op_ptr = new ProfilerOperator(&id_mapper, &query_config);
    break;
  default:
    throw VariantOperationException(std::string("Unknown gvcf sub-operation type: ")
                                    + std::to_string(sub_operation_type) + "n");
  }
  Timer timer;
  timer.start();
  //At least 1 iteration
  VariantQueryProcessorScanState scan_state;
  for (auto i=0u; i<std::max(1u, query_config.get_num_column_intervals()); ++i) {
    while (!scan_state.end()) {
      qp.scan_and_operate(qp.get_array_descriptor(), query_config, *op_ptr, i, true, &scan_state);
      if (serialized_vcf_adapter_ptr) {
        serialized_vcf_adapter_ptr->do_output();
        rw_buffer.m_num_valid_bytes = 0u;
      }
    }
    scan_state.reset();
  }
  timer.stop();
  timer.print(std::string("Total scan_and_produce_Broad_GVCF time")+" for rank "+std::to_string(my_world_mpi_rank), std::cerr);
  if(sub_operation_type == ProduceBroadGVCFSubOperation::PRODUCE_BROAD_GVCF_COUNT_LINES)
    std::cerr << "Count "<< dynamic_cast<ProfilerOperator*>(op_ptr)->get_value() << "\n";
  delete op_ptr;
}

void iterate_columnar_gvcf(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config, int command_idx,
    int sub_operation_type, VCFAdapter& vcf_adapter) {
  switch(sub_operation_type) {
    case ProduceBroadGVCFSubOperation::PRODUCE_BROAD_GVCF_PRODUCE_GVCF:
      {
        BroadCombinedGVCFOperator op(vcf_adapter, query_config.get_vid_mapper(), query_config, false, true, false);
        qp.iterate_over_gvcf_entries(qp.get_array_descriptor(), query_config, op, true);
        break;
      }
    default:
      {
        ProfilerOperator counter;
        qp.iterate_over_gvcf_entries(qp.get_array_descriptor(), query_config, counter, true);
        std::cerr << "Counter "<<counter.get_value() << "\n";
        break;
      }
  }
}

void print_calls(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config, int command_idx, const VidMapper& id_mapper) {
  switch (command_idx) {
  case COMMAND_PRINT_CALLS: {
    std::string indent_prefix = "    ";
    std::cout << "{\n";
    //variant_calls is an array of dictionaries
    std::cout << indent_prefix << "\"variant_calls\": [\n";
    VariantCallPrintOperator printer(std::cout, indent_prefix+indent_prefix, &id_mapper);
    qp.iterate_over_cells(qp.get_array_descriptor(), query_config, printer, true);
    std::cout << "\n" << indent_prefix << "]\n";
    std::cout << "}\n";
    break;
  }
  case COMMAND_PRINT_CSV: {
    VariantCallPrintCSVOperator printer(std::cout);
    qp.iterate_over_cells(qp.get_array_descriptor(), query_config, printer, true);
    break;
  }
  case COMMAND_PRINT_ALT_ALLELE_COUNTS: {
    AlleleCountOperator AC_counter(id_mapper, query_config);
    qp.iterate_over_cells(qp.get_array_descriptor(), query_config, AC_counter, true);
    AC_counter.print_allele_counts();
    break;
  }
  default:
    std::cerr << "Unknown print_calls command "<<command_idx<<"\n";
    exit(-1);
  }
}

void produce_column_histogram(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config, uint64_t bin_size,
                              const std::vector<uint64_t>& num_equi_load_bins) {
  ColumnHistogramOperator histogram_op(0, 4000000000ull, bin_size);
  qp.iterate_over_cells(qp.get_array_descriptor(), query_config, histogram_op, true);
  for (auto val : num_equi_load_bins)
    histogram_op.equi_partition_and_print_bins(val);
}
  
void print_usage() {
  std::cout << "FIXME remove" << std::endl;
  char message[] = "aaaaaasdfsdfkmmmmmmeeeeeaaaaaasdfasdfasdfasdfasdfasdfbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbnml\0";
  std::string str(message);
  std::cout << "Original size: " << str.length() << std::endl;
  char* data;
  size_t data_size;

  void *codec;
  TileDBUtils::create_codec(&codec, TILEDB_ZSTD, Z_DEFAULT_COMPRESSION);
  TileDBUtils::compress(codec, (unsigned char*)message, str.length(), (void**)&data, data_size);
  TileDBUtils::finalize_codec(codec);

  std::cout << "Compressed size " << data_size << std::endl;

  std::cout << "Usage: gt_mpi_gather [options]\n"
            << "where options include:\n"
            << "\t \e[1m--help\e[0m, \e[1m-h\e[0m Print a usage message summarizing options available and exit\n"
            << "\t \e[1m--json-config\e[0m=<query json file>, \e[1m-j\e[0m <query json file>\n"
            << "\t\t Can specify workspace, array, query_column_ranges, query_row_ranges, vid_mapping_file,\n"
            << "\t\t callset_mapping_file, query_attributes, query_filter, reference_genome, etc. as fields in the json file e.g.\n"
            << "\t\t\t { \"workspace\" : \"/tmp/ws\",\n"
            << "\t\t\t   \"array\" : \"t0_1_2\",\n"
            << "\t\t\t   \"query_column_ranges\" : [ [ [0, 100 ], 500 ] ],\n"
            << "\t\t\t   \"query_row_ranges\" : [ [ [0, 2 ] ],\n"
            << "\t\t\t   \"vid_mapping_file\" : \"/tests/inputs/vid.json\",\n"
	    << "\t\t\t   \"callset_mapping_file\": \"/tests/inputs/callset_mapping.json\",\n"
            << "\t\t\t   \"query_attributes\" : [ \"REF\", \"ALT\", \"BaseQRankSum\", \"MQ\", \"MQ0\", \"ClippingRankSum\", \"MQRankSum\", \"ReadPosRankSum\", \"DP\", \"GT\", \"GQ\", \"SB\", \"AD\", \"PL\", \"DP_FORMAT\", \"MIN_DP\" ] }\n"
            << "\t \e[1m--loader-json-config\e[0m=<loader json file>, \e[1m-l\e[0m <loader json file>\n"
            << "\t\t Optional, if vid_mapping_file and callset_mapping_file fields are specified in the query json file\n"
            << "\t \e[1m--workspace\e[0m=<workspace dir>, \e[1m-w\e[0m <GenomicsDB workspace dir>\n"
            << "\t\t Optional, if workspace is specified in any of the json config files\n"
            << "\t \e[1m--array\e[0m=<array dir>, \e[1m-A\e[0m <GenomicsDB array dir>\n"
            << "\t\t Optional, if array is specified in any of the json config files\n"
            << "\t \e[1m--print-calls\e[0m\n"
            << "\t\t Optional, prints VariantCalls in a JSON format\n"
            << "\t \e[1m--print-csv\e[0m\n"
            << "\t\t Optional, outputs CSV with the fields and the order of CSV lines determined by the query attributes\n"
            << "\t \e[1m--produce-bgen\e[0m\n"
            << "\t\t Optional, outputs .bgen v1.2\n"
            << "\t\t \e[1m--bgen-compression=<type>\e[0m\n"
            << "\t\t\t Optional, used with \e[1m--produce-bgen\e[0m, specify the compression method used for the probability data storage\n"
            << "\t\t\t 0 - no compression, 1 - gzip (default), 2 - z standard\n"
            << "\t \e[1m--produce-bed\e[0m\n"
            << "\t\t Optional, outputs plink1.9 .bed, .bim, and .fam file\n"
            << "\t \e[1m--produce-tped\e[0m\n"
            << "\t\t Optional, outputs plink1.9 .tped and .fam files\n"
            << "\t \e[1m--produce-plink\e[0m\n"
            << "\t\t Optional, equivalent to \e[1m--produce-bgen\e[0m, \e[1m--produce-bed\e[0m, and \e[1m--produce-tped\e[0m\n"
            << "\t\t \e[1m--plink-prefix=<prefix>\e[0m\n"
            << "\t\t\t Optional, used with \e[1m--produce-bgen\e[0m, \e[1m--produce-bed\e[0m, \e[1m--produce-tped\e[0m, and \e[1m--produce-plink\e[0m\n"
            << "\t\t\t When specified, all files generated by these options will be named <prefix>.filetype. When this option is omitted <prefix> defaults to \"output\"\n"
            << "\t\t \e[1m--fam-list=<filename>\e[0m\n"
            << "\t\t\t Optional, used with \e[1m--produce-bed\e[0m and \e[1m--produce-tped\e[0m\n"
            << "\t\t\t Takes a file containing the names of .fam files to include this information in the generated .fam file. The output .fam file will only contain samples found by the query.\n"
            << "\t\t\t The within-family id (second column) values of the .fam files should match the sample names in GenomicsDB to be correctly associated\n"
            << "\t\t \e[1m--progress\e[0m[=interval]\n"
            << "\t\t\t Optional, Show query progress (currently implemented for: \e[1m--produce-bgen\e[0m, \e[1m--produce-bed\e[0m, \e[1m--produce-tped\e[0m, and \e[1m--produce-plink\e[0m\n"
            << "\t\t\t specify minimum amount of time between progress messages\n"
            << "\t\t\t where <interval> is a floating point number. Default units are seconds, explicitly specify seconds, minutes, or hours by appending s, m, or h to the end of the number\n"
            << "\t\t \e[1m--plink-one-pass\e[0m\n"
            << "\t\t\t Optional, produce plink files with only one pass, using samples from callset. Improves runtime, but will include samples in query range without any data.\n"
            << "\t \e[1m--produce-Broad-GVCF\e[0m\n"
            << "\t\t Optional, produces combined gVCF from the GenomicsDB data constrained by the query configuration\n"
            << "\t\t \e[1m--output-format\e[0m=<output_format>, \e[1m-O\e[0m <output_format>\n"
            << "\t\t\t used with \e[1m--produce-Broad-GVCF\e[0m\n"
            << "\t\t\t Output format can be one of the following strings: \"z[0-9]\" (compressed VCF),\"b[0-9]\" (compressed BCF)\n"
            << "\t\t\t or \"bu\" (uncompressed BCF). Default is uncompressed VCF if not specified.\n"
            << "\t \e[1m--produce-histogram\e[0m\n"
            << "\t\t Optional\n"
            << "\t \e[1m--produce-interesting-positions\e[0m\n"
            << "\t\t Optional\n"
            << "\t \e[1m--version\e[0m Print version and exit\n" 
            << "\tIf none of the print/produce arguments are specified, the tool prints all the Variants constrained by the query configuration in a JSON format\n\n"
            << "\t\e[1mParallel Querying\e[0m\n"
            << "\t\t MPI could be used for parallel querying, e.g.\n"
            << "\t\t mpirun -n <num_processes> -hostfile <hostfile> ./bin/gt_mpi_gather -j <query.json> -l <loader.json> [<other_args>]\n";
}

void initialize_MPI_env(int *my_world_mpi_rank, int *num_mpi_processes) {
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD, num_mpi_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, my_world_mpi_rank);
#ifdef DEBUG
  //Print host, rank and LD_LIBRARY_PATH
  std::vector<char> hostname;
  hostname.resize(500u);
  gethostname(&(hostname[0]), hostname.size());
  if (my_world_mpi_rank == 0)
    std::cerr << "#processes "<<num_mpi_processes<<"\n";
  auto* ld_library_path_cstr = getenv("LD_LIBRARY_PATH");
  std::cerr << "Host : "<< &(hostname[0]) << " rank "<< my_world_mpi_rank << " LD_LIBRARY_PATH= "<<(ld_library_path_cstr ? ld_library_path_cstr : "")<< "\n";
#endif
}

int main(int argc, char *argv[]) {
  int my_world_mpi_rank = 0;
  int num_mpi_processes = 0;
  initialize_MPI_env(&my_world_mpi_rank, &num_mpi_processes);

  // Define long options
  static struct option long_options[] = {
    {"page-size",1,0,'p'},
    {"rank",1,0,'r'},
    {"output-format",1,0,'O'},
    {"workspace",1,0,'w'},
    {"json-config",1,0,'j'},
    {"loader-json-config",1,0,'l'},
    {"segment-size",1,0,'s'},
    {"skip-query-on-root",0,0,ARGS_IDX_SKIP_QUERY_ON_ROOT},
    {"produce-Broad-GVCF",0,0,ARGS_IDX_PRODUCE_BROAD_GVCF},
    {"produce-interesting-positions",0,0,ARGS_IDX_PRODUCE_INTERESTING_POSITIONS},
    {"produce-histogram",0,0,ARGS_IDX_PRODUCE_HISTOGRAM},
    {"print-calls",0,0,ARGS_IDX_PRINT_CALLS},
    {"print-csv",0,0,ARGS_IDX_PRINT_CSV},
    {"produce-bgen", 0, 0, ARGS_IDX_BGEN},
    {"bgen-compression", 1, 0, ARGS_IDX_BGEN_COMPRESSION},
    {"produce-bed", 0, 0, ARGS_IDX_BED},
    {"produce-tped", 0, 0, ARGS_IDX_TPED},
    {"produce-plink", 0, 0, ARGS_IDX_PLINK},
    {"plink-prefix",1,0,ARGS_IDX_PLINK_PREFIX},
    {"plink-one-pass", 0, 0, ARGS_IDX_PLINK_ONE_PASS},
    {"progress",2,0, ARGS_IDX_PROG},
    {"fam-list", 1, 0, ARGS_IDX_FAM_LIST},
    {"print-AC",0,0,ARGS_IDX_PRINT_ALT_ALLELE_COUNTS},
    {"array",1,0,'A'},
    {"version",0,0,ARGS_IDX_VERSION},
    {"columnar-gvcf",0,0,ARGS_IDX_COLUMNAR_GVCF},
    {"gvcf-profile",0,0,ARGS_IDX_GVCF_PROFILE},
    {"help",0,0,'h'},
    {0,0,0,0},
  };

  int c;
  uint64_t page_size = 0u;
  std::string output_format = "";
  std::string workspace = "";
  std::string array_name = "";
  std::string json_config_file = "";
  std::string loader_json_config_file = "";
  bool skip_query_on_root = false;
  unsigned command_idx = COMMAND_RANGE_QUERY;
  size_t segment_size = 10u*1024u*1024u; //in bytes = 10MB
  auto segment_size_set_in_command_line = false;
  auto sub_operation_type = ProduceBroadGVCFSubOperation::PRODUCE_BROAD_GVCF_UNKNOWN;
  auto use_columnar_iterator = false;
  std::string plink_prefix = "output";
  unsigned char plink_formats = 0;
  bool one_pass = false;
  std::string fam_list = "";
  double progress_interval = -1;
  int bgen_compression = 1;
  while ((c=getopt_long(argc, argv, "j:l:w:A:p:O:s:r:h", long_options, NULL)) >= 0) {
    switch (c) {
    case 'p':
      page_size = strtoull(optarg, 0, 10);
      std::cerr << "WARNING: page size is ignored except for scan now\n";
      break;
    case 'r':
      my_world_mpi_rank = strtoull(optarg, 0, 10);
      break;
    case 'O':
      output_format = std::move(std::string(optarg));
      break;
    case 'w':
      workspace = std::move(std::string(optarg));
      break;
    case 'A':
      array_name = std::move(std::string(optarg));
      break;
    case 's':
      segment_size = strtoull(optarg, 0, 10);
      segment_size_set_in_command_line = true;
      break;
    case ARGS_IDX_SKIP_QUERY_ON_ROOT:
      skip_query_on_root = true;
      break;
    case ARGS_IDX_PRODUCE_BROAD_GVCF:
      command_idx = COMMAND_PRODUCE_BROAD_GVCF;
      sub_operation_type = PRODUCE_BROAD_GVCF_PRODUCE_GVCF;
      break;
    case ARGS_IDX_COLUMNAR_GVCF:
      use_columnar_iterator = true;
      break;
    case ARGS_IDX_PRODUCE_INTERESTING_POSITIONS:
      command_idx = COMMAND_PRODUCE_BROAD_GVCF;
      sub_operation_type = PRODUCE_BROAD_GVCF_PRODUCE_INTERESTING_POSITIONS;
      break;
    case ARGS_IDX_GVCF_PROFILE:
      command_idx = COMMAND_PRODUCE_BROAD_GVCF;
      sub_operation_type = PRODUCE_BROAD_GVCF_COUNT_LINES;
      break;
    case ARGS_IDX_PRODUCE_HISTOGRAM:
      command_idx = COMMAND_PRODUCE_HISTOGRAM;
      break;
    case 'j':
      json_config_file = std::move(std::string(optarg));
      break;
    case ARGS_IDX_PRINT_CALLS:
      command_idx = COMMAND_PRINT_CALLS;
      break;
    case ARGS_IDX_PRINT_CSV:
      command_idx = COMMAND_PRINT_CSV;
      break;
    case ARGS_IDX_BGEN:
      command_idx = COMMAND_PLINK;
      plink_formats = plink_formats | 1;
      break;
    case ARGS_IDX_BGEN_COMPRESSION:
      if (optarg) {
        try {
          bgen_compression = std::stoi(std::string(optarg));
          if(bgen_compression < 0 || bgen_compression > 2) {
            bgen_compression = 1;
          }
        }
        catch ( std::exception& e ) {
        }
      }
    case ARGS_IDX_BED:
      command_idx = COMMAND_PLINK;
      plink_formats = plink_formats | 2;
      break;
    case ARGS_IDX_TPED:
      command_idx = COMMAND_PLINK;
      plink_formats = plink_formats | 4;
      break;
    case ARGS_IDX_PLINK:
      command_idx = COMMAND_PLINK;
      plink_formats = 7;
      break;
    case ARGS_IDX_PLINK_PREFIX:
      if(optarg) {
        plink_prefix = std::string(optarg);
      }
      break;
    case ARGS_IDX_PLINK_ONE_PASS:
      one_pass = true;
      break;
    case ARGS_IDX_PROG:
      if (optarg) {
        try {
          int unit_multiplier = 1;
          std::string optstring(optarg);
          switch(optstring.back()){
            case 's': optstring.pop_back(); break;
            case 'm': unit_multiplier=60; optstring.pop_back(); break;
            case 'h': unit_multiplier=3600; optstring.pop_back(); break;
          }
          progress_interval = (int)(std::stod(std::string(optarg)) * 1000 * unit_multiplier);
        }
        catch ( std::exception& e ) {
          progress_interval = 5000;
        }
      }
      else {
        progress_interval = 5000;
      }
      break;
    case ARGS_IDX_FAM_LIST:
      if(optarg) {
        fam_list = std::string(optarg);
      }
      break;
    case ARGS_IDX_PRINT_ALT_ALLELE_COUNTS:
      command_idx = COMMAND_PRINT_ALT_ALLELE_COUNTS;
      break;
    case 'l':
      loader_json_config_file = std::move(std::string(optarg));
      break;
    case ARGS_IDX_VERSION:
      std::cout << GENOMICSDB_VERSION <<"\n";
      return 0;
    case 'h':
      print_usage();
      return 0;
    default:
      std::cerr << "Unknown command line argument\n";
      print_usage();
      return -1;
    }
  }

  if (json_config_file.empty()) {
    std::cerr << "Query JSON file (-j) is a mandatory argument - unspecified\n";
    print_usage();
    return -1;
  }

  if(use_columnar_iterator) {
    command_idx = COMMAND_COLUMNAR_GVCF;
  }

  int rc=0;
  try {
    //Use VariantQueryConfig to setup query info
    VariantQueryConfig query_config;
    VCFAdapter vcf_adapter_base;
    VCFSerializedBufferAdapter serialized_vcf_adapter(true, true);
    auto& vcf_adapter = (page_size > 0u) ? dynamic_cast<VCFAdapter&>(serialized_vcf_adapter) : vcf_adapter_base;

    //Loader configuration - optional
    GenomicsDBImportConfig loader_config;
    if (!loader_json_config_file.empty()) {
      loader_config.read_from_file(loader_json_config_file, my_world_mpi_rank);
      query_config.update_from_loader(loader_config, my_world_mpi_rank);
    }
    //Info from loader (if specified) is obtained before reading query JSON
    query_config.read_from_file(json_config_file, my_world_mpi_rank);
    //Discard intervals not part of this partition
    if (!loader_json_config_file.empty()) {
      query_config.subset_query_column_ranges_based_on_partition(loader_config, my_world_mpi_rank);
    }
    //Command line overrides
    if (page_size > 0u) {
      query_config.set_combined_vcf_records_buffer_size_limit(page_size);
    }
    if (command_idx == COMMAND_PRODUCE_BROAD_GVCF) {
      if (query_config.get_reference_genome().empty()) {
        throw GenomicsDBConfigException("No reference genome specified in query config");
      }
      query_config.set_vcf_output_format(output_format);
    }
    vcf_adapter.initialize(query_config);
    workspace = query_config.get_workspace(my_world_mpi_rank);
    array_name = query_config.get_array_name(my_world_mpi_rank);
    assert(!workspace.empty() && !array_name.empty());

#ifdef USE_GPERFTOOLS
    ProfilerStart("gt_mpi_gather.gperf.prof");
#endif
#ifdef USE_GPERFTOOLS_HEAP
    HeapProfilerStart("gt_mpi_gather.gperf.heap");
#endif

    segment_size = segment_size_set_in_command_line ? segment_size
        : query_config.get_segment_size();
#if VERBOSE>0
    std::cerr << "Segment size: "<<segment_size<<" bytes\n";
#endif

    /*Create storage manager*/
    VariantStorageManager sm(workspace, segment_size);
    /*Create query processor*/
    VariantQueryProcessor qp(&sm, array_name, query_config.get_vid_mapper());
    auto require_alleles = ((command_idx == COMMAND_RANGE_QUERY)
                            || (command_idx == COMMAND_PRODUCE_BROAD_GVCF)
			    || (command_idx == COMMAND_COLUMNAR_GVCF));
    qp.do_query_bookkeeping(qp.get_array_schema(), query_config, query_config.get_vid_mapper(), require_alleles);
    switch (command_idx) {
      case COMMAND_RANGE_QUERY:
        run_range_query(qp, query_config, query_config.get_vid_mapper(), output_format,
                        (loader_json_config_file.empty() || loader_config.is_partitioned_by_column()),
                        num_mpi_processes, my_world_mpi_rank, skip_query_on_root);
        break;
      case COMMAND_PRODUCE_BROAD_GVCF:
        scan_and_produce_Broad_GVCF(qp, query_config, vcf_adapter, query_config.get_vid_mapper(),
                                    sub_operation_type, my_world_mpi_rank, skip_query_on_root);
        break;
      case COMMAND_COLUMNAR_GVCF:
	iterate_columnar_gvcf(qp, query_config, command_idx, sub_operation_type, vcf_adapter);
	break;
      case COMMAND_PRODUCE_HISTOGRAM:
        produce_column_histogram(qp, query_config, 100, std::vector<uint64_t>({ 128, 64, 32, 16, 8, 4, 2 }));
        break;
      case COMMAND_PRINT_CALLS:
      case COMMAND_PRINT_CSV:
      case COMMAND_PRINT_ALT_ALLELE_COUNTS:
        print_calls(qp, query_config, command_idx, query_config.get_vid_mapper());
        break;
      case COMMAND_PLINK:
        GenomicsDB gdb(query_config.get_workspace(my_world_mpi_rank),
                       query_config.get_callset_mapping_file(),
                       query_config.get_vid_mapping_file(),
                       query_config.get_reference_genome());
        gdb.generate_plink(array_name, &query_config, plink_formats, bgen_compression, one_pass, progress_interval, plink_prefix, fam_list);
        break;
    }
#ifdef USE_GPERFTOOLS_HEAP
    HeapProfilerStop();
#endif
#ifdef USE_GPERFTOOLS
    ProfilerStop();
#endif
    sm.close_array(qp.get_array_descriptor());
    GenomicsDBProtoBufInitAndCleanup::shutdown_protobuf_library();
  } catch(const GenomicsDBConfigException& genomicsdb_ex) {
    std::cerr << genomicsdb_ex.what() << "\n";
    std::cerr << "Do the config files specified to gt_mpi_gather exist? Are they parseable as JSON?\n";
    rc = -1;
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << "\n";
    std::cerr << "Try running gt_mpi_gather --help for usage" << std::endl;
    rc = -1;
  }

  MPI_Finalize();
  return rc;
}
