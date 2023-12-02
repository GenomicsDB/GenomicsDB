/**
 * src/test/cpp/src/test_genomicsdb_api.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2019-2023 Omics Data Automation, Inc.
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
 *
 * Test the GenomicsDB query api
 *
 */

#include <catch2/catch.hpp>

#include "genomicsdb.h"
#include "genomicsdb_config_base.h"
#include "tiledb_utils.h"

#include "genomicsdb_export_config.pb.h"

#include "test_base.h"

#include <iostream>
#include <string>
#include <thread>
#include <utility>

TEST_CASE_METHOD(TempDir, "utils", "[genomicsdb_utils]") {
  REQUIRE(genomicsdb::version().size() > 0);
  REQUIRE_THAT(genomicsdb::version(), Catch::Equals(GENOMICSDB_VERSION));

  std::string hello_world = "Hello World";
  std::string filename =  append("hello.txt");
  CHECK(TileDBUtils::write_file(filename, hello_world.data(), hello_world.length(), true) == TILEDB_OK);

  CHECK(genomicsdb::is_file(filename));
  CHECK(!genomicsdb::is_file(filename+".nonexistent"));
  
  CHECK(genomicsdb::file_size(filename) == 11);
  CHECK(genomicsdb::file_size(filename+".new") == -1);
  
  char *txt;
  size_t length;
  CHECK(genomicsdb::read_entire_file(filename, (void **)&txt, &length) == GENOMICSDB_OK);
  CHECK(length == 11);
  CHECK(hello_world == txt);
  CHECK(genomicsdb::read_entire_file(filename+".another", (void **)txt, &length) == TILEDB_ERR);
}

TEST_CASE("api empty_args", "[empty_args]") {
  std::string empty_string;

  // Constructor that allows workspace/callset/vidmapping/reference genome specification
  CHECK_THROWS_AS(new GenomicsDB(empty_string, empty_string, empty_string), std::exception);
  CHECK_THROWS_AS(new GenomicsDB("ws", empty_string, empty_string), std::exception);
  CHECK_THROWS_AS(new GenomicsDB("ws", "callset", empty_string), std::exception);

  // Constructor that allows json files for configuration
  CHECK_THROWS_AS(new GenomicsDB(empty_string), std::exception);
  CHECK_THROWS_AS(new GenomicsDB(empty_string, GenomicsDB::JSON_FILE), std::exception);
  CHECK_THROWS_AS(new GenomicsDB(empty_string, GenomicsDB::JSON_STRING, "loader_config.json"), std::exception);
  CHECK_THROWS_AS(new GenomicsDB(empty_string, GenomicsDB::PROTOBUF_BINARY_STRING, "loader_config.json"), std::exception);

  // Check if there is an error message
  try {
    new GenomicsDB(empty_string);
    FAIL();
  } catch (GenomicsDBException& e) {
    REQUIRE(e.what());
    CHECK(strlen(e.what()) > 0);
  }
}

static std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

static std::string workspace(ctests_input_dir+"ws");
static std::string callset_mapping(ctests_input_dir+"callset_t0_1_2.json");
static std::string vid_mapping(ctests_input_dir+"vid.json");
static std::string reference_genome(ctests_input_dir+"chr1_10MB.fasta.gz");

static std::string query_json(ctests_input_dir+"query.json");
static std::string loader_json(ctests_input_dir+"loader.json");

// Test Phased Ploidy
static std::string workspace_PP(ctests_input_dir+"ws_phased_ploidy");
static std::string vid_mapping_PP(ctests_input_dir+"vid_phased_ploidy.json");

// Test workspace created with vcf2genomicsdb_init
static std::string workspace_new(ctests_input_dir+"ws_new");

static std::string array("t0_1_2");

// Test Shared PosixFS Optimizations
static std::string query_with_shared_posixfs_optimizations(ctests_input_dir+"query_with_shared_posixfs_optimizations.json");
static std::string consolidation_lock_file(workspace+"/"+array+"/"+".__consolidation_lock");

void check_query_variants_results(GenomicsDB* gdb, const std::string& array, GenomicsDBVariants variants) {
   REQUIRE(variants.size() == 4);
   CHECK(variants.at(5) == nullptr);
   auto variant = variants.next();
   auto variant1 = variants.at(0);
   REQUIRE(variant == variant1);
   auto interval = gdb->get_interval(variant);
   CHECK(interval.first == 12140);
   CHECK(interval.second == 12294);
   auto genomic_interval = gdb->get_genomic_interval(variant);
   CHECK(genomic_interval.contig_name == "1");
   CHECK(genomic_interval.interval.first == 12141);
   CHECK(genomic_interval.interval.second == 12295);
   auto genomic_fields = gdb->get_genomic_fields(array, variant);
   CHECK(genomic_fields.size() == 0);

   // Check variant calls
   auto variant_calls = gdb->get_variant_calls(array, variant);
   REQUIRE(variant_calls.size() == 1);
   CHECK(variant_calls.at(1) == nullptr);
   auto variant_call = variant_calls.next();
   auto variant_call1 = variant_calls.at(0);
   REQUIRE(variant_call == variant_call1);
   CHECK(gdb->get_row(variant_call) == 0);
   auto call_interval = gdb->get_interval(variant_call);
   CHECK(call_interval.first == 12140);
   CHECK(call_interval.second == 12294);
   auto call_genomic_interval = gdb->get_genomic_interval(variant_call);
   CHECK(call_genomic_interval.contig_name == "1");
   CHECK(call_genomic_interval.interval.first == 12141);
   CHECK(call_genomic_interval.interval.second == 12295);
   auto call_genomic_fields = gdb->get_genomic_fields(array, variant_call);
   CHECK(call_genomic_fields.size() >= 2);
   CHECK(call_genomic_fields[0].name == "REF");
   CHECK(variant_calls.get_genomic_field_type(call_genomic_fields[0].name).is_string());
   CHECK(call_genomic_fields[0].str_value() == "C");
   CHECK(call_genomic_fields[1].name == "ALT");
   CHECK(variant_calls.get_genomic_field_type(call_genomic_fields[1].name).is_cppstring());
   CHECK(call_genomic_fields[1].to_string(variant_calls.get_genomic_field_type(call_genomic_fields[1].name)) == "&");
}

TEST_CASE("api query_variants direct DP", "[query_variants]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}, {{0,3}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array));
  delete gdb;
}

TEST_CASE("api query_variants direct DP with huge row range", "[query_variants_huge_range]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}, {{0,10000000000}}));
  delete gdb;
}

TEST_CASE("api query_variants direct DP with huge disjoint row ranges", "[query_variants_huge_disjoint_range]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,10000000000}}, {{0, 1}, {2,1000000000}, {2000000000, 3000000000}}));
  delete gdb;
}

TEST_CASE("api query_variants direct DP with range composed of multiple points", "[query_variants_points]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);
  REQUIRE(gdb->query_variants(array, {{0,10000000000}}, {{0, 0}, {2,2}}).size() == 2);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,10000000000}}, {{0, 0}, {2,2}, {1,1}}));
  delete gdb;
}

TEST_CASE("api query_variants direct DP with range composed of overlapping ranges", "[query_variants_overlap]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);
  REQUIRE(gdb->query_variants(array, {{0,10000000000}}, {{1000,2000}, {1, 5}, {2,2}}).size() == 3);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,10000000000}}, {{0, 10}, {1,5}, {2,3}}));
  delete gdb;
}

TEST_CASE("api query_variants direct DP with invalid row range", "[query_variants_invalid_range]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);
  REQUIRE_THROWS(check_query_variants_results(gdb, array, gdb->query_variants(array, {{-100, 1}}, {{0,10000000000}})));
  delete gdb;
}

TEST_CASE("api query_variants direct DP and GT", "[query_variants_DP_GT]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP", "GT"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}, {{0,3}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array));
  delete gdb;
}

TEST_CASE("api query_variants direct DP and GT with PP", "[query_variants_DP_GT_with_PP]") {
  GenomicsDB* gdb = new GenomicsDB(workspace_PP, callset_mapping, vid_mapping_PP, {"DP", "GT"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}, {{0,3}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array));
  delete gdb;
}

TEST_CASE("api query_variants with json", "[query_variants_json]") {
  GenomicsDB* gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json);
  check_query_variants_results(gdb, array, gdb->query_variants());
  delete gdb;

  gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json, 0);
  check_query_variants_results(gdb, array, gdb->query_variants());
  delete gdb;

  CHECK_THROWS_AS(new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json, 1), GenomicsDBConfigException);
}

TEST_CASE("api query variants with json string", "[query_variants_json_string]") {
  char *query_json_buffer=NULL;
  size_t query_json_buffer_length;
  CHECK(TileDBUtils::read_entire_file(query_json, (void **)&query_json_buffer, &query_json_buffer_length) == TILEDB_OK);
  CHECK(query_json_buffer);
  CHECK(query_json_buffer_length > 0);

  GenomicsDB* gdb = new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json);
  check_query_variants_results(gdb, array, gdb->query_variants());
  // Consolidation lock file created as shared posixfs optimization is not set by default
  std::cerr << consolidation_lock_file << std::endl;
  CHECK(TileDBUtils::is_file(consolidation_lock_file));
  delete gdb;

  gdb = new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json, 0);
  check_query_variants_results(gdb, array, gdb->query_variants());
  delete gdb;

  CHECK_THROWS_AS(new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json, 1), GenomicsDBConfigException);

  free(query_json_buffer);
}

TEST_CASE("api query variants with enabled shared posixfs optimizations", "[query_variants_json_string_posixfs_optimizations]") {
  char *query_json_buffer=NULL;
  size_t query_json_buffer_length;
  CHECK(TileDBUtils::read_entire_file(query_with_shared_posixfs_optimizations, (void **)&query_json_buffer, &query_json_buffer_length) == TILEDB_OK);
  CHECK(query_json_buffer);
  CHECK(query_json_buffer_length > 0);

  if (TileDBUtils::is_file(consolidation_lock_file)) {
    if (TileDBUtils::delete_file(consolidation_lock_file) == TILEDB_OK) {
      GenomicsDB* gdb = new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json);
      check_query_variants_results(gdb, array, gdb->query_variants());
      // A consolidation lock file should not be created when querying with shared posixfs
      // optimization is true
      CHECK(!(TileDBUtils::is_file(consolidation_lock_file)));
      delete gdb;
    }
  }

  free(query_json_buffer);
}

class NullVariantCallProcessor : public GenomicsDBVariantCallProcessor {
 public:
  NullVariantCallProcessor() {
  };

  void process(const interval_t& interval) {
  };

  void process(const std::string& sample_name,
               const int64_t* coordinates,
               const genomic_interval_t& genomic_interval,
               const std::vector<genomic_field_t>& genomic_fields) {
  };
};

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

  void check_single_match() {
    CHECK(m_intervals == 1);
    CHECK(m_count == 1);
    CHECK(m_coordinates[0] == 1);
    CHECK(m_coordinates[1] == 17384);
  }

  int m_intervals = 0;
  int m_count = 0;
  int64_t m_coordinates[2];
};

class OneQueryIntervalProcessor : public GenomicsDBVariantCallProcessor {
 public:
  OneQueryIntervalProcessor(const std::vector<std::string> attributes = {}, bool is_PP=false) {
    m_attributes_size = attributes.size();
    m_is_PP = is_PP;
  }

  ~OneQueryIntervalProcessor() {
    CHECK(m_num_query_intervals == 1);
    CHECK(m_processed_rows == 5);
    CHECK(m_sample_found);
  }

  void reinitialize(int max_int_check =false) {
    m_num_query_intervals = 0;
    m_processed_rows = 0;
    m_sample_found = false;
    m_max_int_check = max_int_check;
  }

  void process(const interval_t& interval) {
    m_num_query_intervals++;
    CHECK(interval.first == 0);
    if (m_max_int_check) {
      CHECK(interval.second == INT64_MAX - 1);
    } else {
      CHECK(interval.second == 1000000000);
    }
  };

  void process(const std::string& sample_name,
               const int64_t* coordinates,
               const genomic_interval_t& genomic_interval,
               const std::vector<genomic_field_t>& genomic_fields) {
    if (sample_name == "HG01530") {
      m_sample_found = true;
    }
    m_processed_rows++;
    auto row = coordinates[0];
    CHECK(row <= 2);
    if (row == 2) {
      CHECK(sample_name == "HG01530");
      CHECK(genomic_interval.contig_name == "1");
      CHECK(genomic_interval.interval.first == 17385);
      CHECK(genomic_interval.interval.second == 17385);
      if (m_attributes_size) {
        CHECK(genomic_fields.size() == m_attributes_size + 2); // Add '2' for REF and ALT fields
      }
      size_t found_fields = 0;
      for (auto i=0u; i<genomic_fields.size(); i++) {
        if (genomic_fields[i].name == "REF") {
          found_fields++;
          CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
          CHECK(genomic_fields[i].str_value() == "G");
        } else if (genomic_fields[i].name == "ALT") {
          found_fields++;
          CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
          CHECK(genomic_fields[i].str_value() == "A|&");
          CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[1].name)) == "[A, <NON_REF>]");
        } else if (genomic_fields[i].name == "DP") {
          found_fields++;
          CHECK(get_genomic_field_type(genomic_fields[i].name).is_int());
          CHECK(genomic_fields[i].int_value_at(0) == 76);
          CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[2].name)) == "76");
        } else if (genomic_fields[i].name == "GT") {
          found_fields++;
          CHECK(get_genomic_field_type(genomic_fields[i].name).is_int());
          if (m_is_PP) {
            CHECK(genomic_fields[i].get_num_elements() == 3);
          } else {
            CHECK(genomic_fields[i].get_num_elements() == 2);
            CHECK(genomic_fields[i].int_value_at(0) == 0);
            CHECK(genomic_fields[i].int_value_at(1) == 1);
            CHECK_THROWS_AS(genomic_fields[i].int_value_at(2), std::exception);
          }
          CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "0/1");
        }
      }
      if (m_attributes_size) {
        CHECK(found_fields == m_attributes_size + 2);
      }
    }
    if (sample_name == "HG01958" && coordinates[1] == 12144) {
      for (auto i=0u; i<genomic_fields.size(); i++) {
        if (genomic_fields[i].name == "GT") {
          if (m_is_PP) {
            CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "0|0");
          } else {
            CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "0/0");
          }
        }
      }
    }
  };

  size_t m_attributes_size;
  int m_num_query_intervals = 0;
  int m_processed_rows = 0;
  int m_sample_found = false;
  int m_is_PP;
  bool m_max_int_check = false;
};

class TwoQueryIntervalsProcessor : public GenomicsDBVariantCallProcessor {
 public:
  TwoQueryIntervalsProcessor() {};

  ~TwoQueryIntervalsProcessor() {
    CHECK(m_num_query_intervals == 2);
    CHECK(m_processed_rows == 5);
  }

  void reinitialize() {
    m_num_query_intervals = 0;
    m_processed_rows = 0;
  }

  void process(const interval_t& interval) {
    m_num_query_intervals++;
    if (m_num_query_intervals == 1) {
      CHECK(interval.first == 0);
      CHECK(interval.second == 17000);
    } else {
      CHECK(interval.first == 17000);
      CHECK(interval.second == 18000);
    }
  };

  void process(const std::string& sample_name,
               const int64_t* coordinates,
               const genomic_interval_t& genomic_interval,
               const std::vector<genomic_field_t>& genomic_fields) {
    m_processed_rows++;
  };

  int m_num_query_intervals = 0;
  int m_processed_rows = 0;
};

TEST_CASE("api query_variant_calls direct", "[query_variant_calls_direct]") {
  const std::vector<std::string> attributes = {"DP"};

  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, attributes, 40);

  gdb->query_variant_calls(array);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor(attributes);
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}}, {{0,3}});
  one_query_interval_processor.reinitialize();
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}});
  one_query_interval_processor.reinitialize(true);
  gdb->query_variant_calls(one_query_interval_processor, array);

  TwoQueryIntervalsProcessor two_query_intervals_processor;
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}}, {{0,3}});
  two_query_intervals_processor.reinitialize();
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}});

  delete gdb;
}

TEST_CASE("api query_variant_calls direct DP and GT", "[query_variant_calls_direct_DP_GT]") {
  const std::vector<std::string> attributes = {"GT", "DP"};

  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, attributes, 40);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor(attributes);
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}}, {{0,3}});
  one_query_interval_processor.reinitialize();
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}});
  one_query_interval_processor.reinitialize(true);
  gdb->query_variant_calls(one_query_interval_processor, array);

  TwoQueryIntervalsProcessor two_query_intervals_processor;
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}}, {{0,3}});
  two_query_intervals_processor.reinitialize();
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}});

  delete gdb;
}

TEST_CASE("api query_variant_calls direct DP and GT with PP", "[query_variant_calls_direct_DP_GT_with_PP]") {
  const std::vector<std::string> attributes = {"GT", "DP"};

  GenomicsDB* gdb = new GenomicsDB(workspace_PP, callset_mapping, vid_mapping_PP, attributes, 40);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor(attributes, true);
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}}, {{0,3}});
  one_query_interval_processor.reinitialize();
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}});
  one_query_interval_processor.reinitialize(true);
  gdb->query_variant_calls(one_query_interval_processor, array);

  TwoQueryIntervalsProcessor two_query_intervals_processor;
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}}, {{0,3}});
  two_query_intervals_processor.reinitialize();
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}});

  delete gdb;
}

TEST_CASE("api query_variant_calls with json", "[query_variant_calls_with_json]") {
  GenomicsDB* gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json);
  gdb->query_variant_calls();

  OneQueryIntervalProcessor one_query_interval_processor;
  gdb->query_variant_calls(one_query_interval_processor, "", GenomicsDB::NONE);
  delete gdb;

  OneQueryIntervalProcessor another_query_interval_processor;
  gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json, 0);
  gdb->query_variant_calls(another_query_interval_processor, "", GenomicsDB::NONE);
  delete gdb;

  CHECK_THROWS_AS(new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json, 1), GenomicsDBConfigException);
}

TEST_CASE("api query_variant_calls with protobuf", "[query_variant_calls_with_protobuf]") {
  using namespace genomicsdb_pb;

  ExportConfiguration *config = new ExportConfiguration();

  config->set_workspace(workspace);
  config->set_callset_mapping_file(callset_mapping);
  config->set_vid_mapping_file(vid_mapping);

  // query_row_ranges
  RowRangeList* row_ranges = config->add_query_row_ranges();
  RowRange* row_range = row_ranges->add_range_list();
  row_range->set_low(0);
  row_range->set_high(3);

  // query_column_ranges
  GenomicsDBColumnOrIntervalList* column_ranges = config->add_query_column_ranges();
  GenomicsDBColumnOrInterval* column_range = column_ranges->add_column_or_interval_list();

  TileDBColumnInterval* tiledb_column_interval = new TileDBColumnInterval();
  tiledb_column_interval->set_begin(0);
  tiledb_column_interval->set_end(1000000000);

  GenomicsDBColumnInterval* column_interval = new GenomicsDBColumnInterval();
  column_interval->set_allocated_tiledb_column_interval(tiledb_column_interval);

  column_range->set_allocated_column_interval(column_interval);

  // query_attributes
  config->add_attributes()->assign("GT");
  config->add_attributes()->assign("DP");

  config->set_reference_genome("inputs/chr1_10MB.fasta.gz");
  config->set_segment_size(40);
  config->set_vcf_header_filename("inputs/template_vcf_header.vcf");

  std::string config_string;
  CHECK(config->SerializeToString(&config_string));

  OneQueryIntervalProcessor one_query_interval_processor;
  
  // Try with no array name set, GenomicsDB should find the single array t0_1_2 during query
  GenomicsDB* gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  gdb->query_variant_calls();
  gdb->query_variant_calls(one_query_interval_processor, "", GenomicsDB::NONE);
  delete gdb;
  
  config->set_array_name(array);
  CHECK(config->SerializeToString(&config_string));
  gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  gdb->query_variant_calls();
  one_query_interval_processor.reinitialize();
  gdb->query_variant_calls(one_query_interval_processor, "", GenomicsDB::NONE);
  delete gdb;

  // Filter expressions are supported only for workspaces/arrays with ploidy information
  config->set_workspace(workspace_PP);

  SECTION("Try query with contig intervals instead of tiledb column intervals") {
    ContigInterval* contig_interval = new ContigInterval();
    contig_interval->set_contig("1");
    contig_interval->set_begin(1);
    contig_interval->set_end(249250621);
    column_interval->Clear();
    column_interval->set_allocated_contig_interval(contig_interval);
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 5);
    delete gdb;
  }

  config->set_query_filter("REF == \"G\" && resolve(GT, REF, ALT) &= \"T/T\" && ALT |= \"T\"");
  SECTION("Try query with a filter") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    count_cells_processor.check_single_match();
    delete gdb;
  }

  config->set_workspace(workspace_PP);
  SECTION("Try query with a filter and workspace with GT ploidy") {
    config->set_bypass_intersecting_intervals_phase(true);
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    count_cells_processor.check_single_match();
    delete gdb;
  }

  SECTION("Try query with a filter and a small segment size to force TileDB buffers overflow") {
    config->set_segment_size(40);
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    count_cells_processor.check_single_match();
    delete gdb;
  }

  config->set_query_filter("POS==99999999 && ROW==1");
  SECTION("Try query with POS and ROW with no match") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 0);
    CHECK(count_cells_processor.m_count == 0);
    delete gdb;
  }


  config->set_query_filter("POS==17384 && ROW==1 && ISHOMALT && resolve(GT, REF, ALT) &= \"T/T\"");
  SECTION("Try query with POS and ROW") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    count_cells_processor.check_single_match();
    delete gdb;
  }
  
  config->set_query_filter("ISHOMREF");
  SECTION("Try query with ISHOMREF") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 2);
    delete gdb;
  }

  config->set_query_filter("!ISHOMREF");
  SECTION("Try query with ISHOMREF") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 3);
    delete gdb;
  }

  config->set_query_filter("ISHOMALT");
  SECTION("Try query with ISHOMREF") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 1);
    delete gdb;
  }

  config->set_query_filter("ISHET");
  SECTION("Try query with ISHOMREF") {
    CHECK(config->SerializeToString(&config_string));
    gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 2);
    delete gdb;
  }

}

TEST_CASE("api query_variant_calls with protobuf and config", "[query_variant_calls_with_protobuf_and_explicit_configuration]") {
  using namespace genomicsdb_pb;

  ExportConfiguration *config = new ExportConfiguration();

  config->set_workspace(workspace_new);
  config->set_array_name("t0_1_2");
  config->set_callset_mapping_file(workspace_new+"/callset.json");
  config->set_vid_mapping_file(workspace_new+"/vidmap.json");
  config->set_enable_shared_posixfs_optimizations(true);
  config->set_bypass_intersecting_intervals_phase(true);
  // query_attributes
  config->add_attributes()->assign("REF");
  config->add_attributes()->assign("ALT");
  config->add_attributes()->assign("GT");
  config->set_segment_size(128);

  std::string config_string;
  CHECK(config->SerializeToString(&config_string));
  GenomicsDB *gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING);

  SECTION("Base configuration") {
    CountCellsProcessor count_cells_processor;
    gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 5);
  }

  SECTION("Subset Array with pb string that is not Query or ExportConfiguration") {
    CountCellsProcessor count_cells_processor;
    CHECK_THROWS_AS(gdb->query_variant_calls(count_cells_processor, "NotAProtobufConfiguration",
                                             GenomicsDB::PROTOBUF_BINARY_STRING), std::exception);
  }

  SECTION("Subset Array using QueryConfiguration") {
    CountCellsProcessor count_cells_processor;
    QueryConfiguration query_config;
    query_config.set_array_name("t0_1_2");
    ContigInterval* contig_interval = query_config.add_query_contig_intervals();
    contig_interval->set_contig("1");
    contig_interval->set_begin(1);
    contig_interval->set_end(100000);
    RowRangeList* row_ranges = query_config.add_query_row_ranges();
    RowRange* row_range = row_ranges->add_range_list();
    row_range->set_low(0);
    row_range->set_high(3);
    std::string query_string;
    CHECK(query_config.SerializeToString(&query_string));
    gdb->query_variant_calls(count_cells_processor, query_string, GenomicsDB::PROTOBUF_BINARY_STRING);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 5);
  }

  SECTION("Subset Array using ExportConfiguration") {
    CountCellsProcessor count_cells_processor;
    ExportConfiguration query_config;
    ContigInterval* contig_interval = query_config.add_query_contig_intervals();
    contig_interval->set_contig("1");
    contig_interval->set_begin(1);
    contig_interval->set_end(100000);
    std::string query_string;
    CHECK(query_config.SerializeToString(&query_string));
    gdb->query_variant_calls(count_cells_processor, query_string, GenomicsDB::PROTOBUF_BINARY_STRING);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 5);
  }

  SECTION("Subset Array using QueryConfiguration - two intervals") {
    CountCellsProcessor count_cells_processor;
    QueryConfiguration query_config;
    ContigInterval* contig_interval = query_config.add_query_contig_intervals();
    contig_interval->set_contig("1");
    contig_interval->set_begin(1);
    contig_interval->set_end(17000);
    contig_interval = query_config.add_query_contig_intervals();
    contig_interval->set_contig("1");
    contig_interval->set_begin(17000);
    contig_interval->set_end(18000);
    std::string query_string;
    CHECK(query_config.SerializeToString(&query_string));
    gdb->query_variant_calls(count_cells_processor, query_string, GenomicsDB::PROTOBUF_BINARY_STRING);
    CHECK(count_cells_processor.m_intervals == 2);
    CHECK(count_cells_processor.m_count == 5);
  }

   SECTION("Subset Array using QueryConfiguration - multiple queries") {
    QueryConfiguration query_config;

    // First query
    CountCellsProcessor count_cells_processor;
    ContigInterval* contig_interval = query_config.add_query_contig_intervals();
    contig_interval->set_contig("1");
    contig_interval->set_begin(1);
    contig_interval->set_end(17000);
    std::string query_string;
    CHECK(query_config.SerializeToString(&query_string));
    gdb->query_variant_calls(count_cells_processor, query_string, GenomicsDB::PROTOBUF_BINARY_STRING);
    CHECK(count_cells_processor.m_intervals == 1);
    CHECK(count_cells_processor.m_count == 2);

    // Second query
    CountCellsProcessor another_count_cells_processor;
    contig_interval->set_contig("1");
    contig_interval->set_begin(17000);
    contig_interval->set_end(18000);
    CHECK(query_config.SerializeToString(&query_string));
    gdb->query_variant_calls(another_count_cells_processor, query_string, GenomicsDB::PROTOBUF_BINARY_STRING);
    CHECK(another_count_cells_processor.m_intervals == 1);
    CHECK(another_count_cells_processor.m_count == 3);
  }

  delete gdb;
}

TEST_CASE("api query_variant_calls with protobuf new", "[query_variant_calls_with_protobuf_new]") {
  using namespace genomicsdb_pb;

  ExportConfiguration *config = new ExportConfiguration();

  config->set_workspace(workspace_new);
  config->set_array_name("t0_1_2");
  config->set_callset_mapping_file(workspace_new+"/callset.json");
  config->set_vid_mapping_file(workspace_new+"/vidmap.json");

  // query_row_ranges
  RowRangeList* row_ranges = config->add_query_row_ranges();
  RowRange* row_range = row_ranges->add_range_list();
  row_range->set_low(0);
  row_range->set_high(3);

  // query contig_interval
  ContigInterval* contig_interval = config->add_query_contig_intervals();
  contig_interval->set_contig("1");
  contig_interval->set_begin(1);
  contig_interval->set_end(249250621);

  // query_attributes
  config->add_attributes()->assign("GT");
  config->add_attributes()->assign("DP");

  config->set_bypass_intersecting_intervals_phase(true);
  config->set_query_filter("REF == \"G\" && resolve(GT, REF, ALT) &= \"T/T\" && ALT |= \"T\"");

  size_t segment_sizes[4] = { 16, 36, 96, 128 };

  for (auto i=0; i<4; i++) {
    SECTION("Test filter test for " + std::to_string(i)) {
      config->set_segment_size(segment_sizes[i]);

      std::string config_string;
      CHECK(config->SerializeToString(&config_string));
      GenomicsDB *gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING);

      CountCellsProcessor count_cells_processor;
      gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);

      CHECK(count_cells_processor.m_intervals == 1);
      CHECK(count_cells_processor.m_count == 1);
      delete gdb;
    }
  }
}

TEST_CASE("api query_variant_calls with JSONVariantCallProcessor", "[query_variant_calls_with_json_processor]") {
  using namespace genomicsdb_pb;

  ExportConfiguration *config = new ExportConfiguration();

  config->set_workspace(workspace_new);
  config->set_array_name("t0_1_2");
  config->set_callset_mapping_file(workspace_new+"/callset.json");
  config->set_vid_mapping_file(workspace_new+"/vidmap.json");

  // query_row_ranges
  RowRangeList* row_ranges = config->add_query_row_ranges();
  RowRange* row_range = row_ranges->add_range_list();
  row_range->set_low(0);
  row_range->set_high(3);

  // query contig_interval
  ContigInterval* contig_interval = config->add_query_contig_intervals();
  contig_interval->set_contig("1");
  contig_interval->set_begin(1);
  contig_interval->set_end(249250621);

  // query_attributes
  config->add_attributes()->assign("GT");
  config->add_attributes()->assign("DP");

  config->set_bypass_intersecting_intervals_phase(true);

  config->set_segment_size(128);

  std::string config_string;
  CHECK(config->SerializeToString(&config_string));
  GenomicsDB *gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING);

  JSONVariantCallProcessor json_processor;
  gdb->query_variant_calls(json_processor, "", GenomicsDB::NONE);

  auto output = json_processor.construct_json_output();
  // {"HG00141":{"CHROM":["1","1"],"POS":[12141,17385],"DP":[".","."],"GT":["C/C","G/A"]},"HG01530":{"CHROM":["1"],"POS":[17385],"DP":[76],"GT":["G/A"]},"HG01958":{"CHROM":["1","1"],"POS":[12145,17385],"DP":[".",120],"GT":["C/C","T/T"]}}
  CHECK(strlen(output) == 232);
  printf("%lu %s\n\n", strlen(output), output);

  JSONVariantCallProcessor json_processor1(JSONVariantCallProcessor::samples_with_ncalls);
  gdb->query_variant_calls(json_processor1, "", GenomicsDB::NONE);
  output = json_processor1.construct_json_output();
  // {"HG00141":2,"HG01530":1,"HG01958":2}
  CHECK(strlen(output) == 37);
  printf("%lu %s\n\n", strlen(output), output);

  JSONVariantCallProcessor json_processor2(JSONVariantCallProcessor::just_ncalls);
  gdb->query_variant_calls(json_processor2, "", GenomicsDB::NONE);
  output = json_processor2.construct_json_output();
  // {"num_calls":5}
  CHECK(strlen(output)== 15);
  printf("%lu %s\n\n", strlen(output), output);

  JSONVariantCallProcessor json_processor3;
  json_processor3.set_payload_mode(JSONVariantCallProcessor::just_ncalls);
  gdb->query_variant_calls(json_processor3, "", GenomicsDB::NONE);
  output = json_processor3.construct_json_output();
  // {"num_calls":5}
  CHECK(strlen(output)== 15);
  printf("%lu %s\n\n", strlen(output), output);

  delete gdb;
}


TEST_CASE("Test genomicsdb demo test case", "[genomicsdb_demo]") {
  using namespace genomicsdb_pb;
  
  char *genomicsdb_demo_workspace = getenv("GENOMICSDB_DEMO_WS");
  if (!genomicsdb_demo_workspace) return;

  ExportConfiguration *config = new ExportConfiguration();

  std::string ws(genomicsdb_demo_workspace);
  config->set_workspace(ws);
  config->set_array_name("allcontigs$1$3095677412");
  config->set_callset_mapping_file(ws+"/callset.json");
  config->set_vid_mapping_file(ws+"/vidmap.json");

  // query_contig_intervals
  auto* contig_interval = config->add_query_contig_intervals();
  contig_interval->set_contig("17");
  contig_interval->set_begin(7571719);
  contig_interval->set_end(7590868);

  // query_row_ranges
  auto* row_range = config->add_query_row_ranges()->add_range_list();
  row_range->set_low(0);
  row_range->set_high(200000);

  // query_attributes
  config->add_attributes()->assign("REF");
  config->add_attributes()->assign("ALT");
  config->add_attributes()->assign("GT");

  // other
  config->set_bypass_intersecting_intervals_phase(true);
  config->set_enable_shared_posixfs_optimizations(true);

  // filters
  std::vector<std::string> filters = {""/*, // zlib arm64 - 1s
    "REF==\"A\"", // 2s
    "REF==\"A\" && ALT|=\"T\"", // 2s
    "REF==\"A\" &&  ALT|=\"T\" && ISHOMALT", // 3s
    "REF==\"A\" && ALT|=\"T\" && resolve(GT, REF, ALT) &= \"T|T\"" // 3s*/
  };

  // sizes
  std::vector<size_t> segment_sizes = { 0/*, 10240, 20480, 40960*/ };

  // results
  std::vector<int64_t> counts = { 2962039, 400432, 82245, 82245, 69548 };

  for (auto i=0u; i<filters.size(); i++) {
    for (auto j=0u; j<segment_sizes.size(); j++) {
      SECTION("Demo test for " + std::to_string(i) + std::to_string(j)) {
        if (segment_sizes[j] > 0) config->set_segment_size(segment_sizes[j]);
        config->set_query_filter(filters[i]);

        std::string config_string;
        CHECK(config->SerializeToString(&config_string));

        Catch::Timer t;
        t.start();

        GenomicsDB* gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING);
        CountCellsProcessor count_cells_processor;
        gdb->query_variant_calls(count_cells_processor, "", GenomicsDB::NONE);

        CHECK(count_cells_processor.m_intervals == 1);
        CHECK(count_cells_processor.m_count == counts[i]);
        printf("Elapsed Time=%us for filter=%s segment_size=%zu\n", t.getElapsedMilliseconds()/1000,
               filters[i].c_str(), segment_sizes[j]);
  
        delete gdb;
      }
    }
  }
}

TEST_CASE("api generate_vcf direct", "[query_generate_vcf_direct]") {
  TempDir temp_dir;
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, {"DP"}, 40);

  const std::string vcf_file1 = temp_dir.append("1.vcf.gz");
  gdb->generate_vcf(array, {{0,1000000000}}, {{0,3}}, reference_genome, "", vcf_file1);
  CHECK(TileDBUtils::is_file(vcf_file1));
  CHECK(!TileDBUtils::is_file(vcf_file1+".tbi"));

  gdb->generate_vcf(array, {{0,1000000000}}, {{0,3}}, reference_genome, "", vcf_file1, "z", true);
  CHECK(TileDBUtils::is_file(vcf_file1));
  CHECK(TileDBUtils::is_file(vcf_file1+".tbi"));

  CHECK_THROWS_AS(gdb->generate_vcf(array, {{0,1000000000}}, {{0,3}}, reference_genome, "", vcf_file1, "z", false), std::exception);

  const std::string vcf_file2 = temp_dir.append("3.vcf.gz");
  gdb->generate_vcf(array, {{0,13000}, {13000, 1000000000}}, {{0,3}}, reference_genome, "", vcf_file2, "z", true);
  CHECK(TileDBUtils::is_file(vcf_file2));
  CHECK(TileDBUtils::is_file(vcf_file2+".tbi"));

  delete gdb;
}

TEST_CASE("api generate_vcf with json", "[query_generate_with_json]") {
  TempDir temp_dir;
  GenomicsDB* gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json);

  const std::string vcf_file = temp_dir.append("1.vcf.gz");
  gdb->generate_vcf(vcf_file, "z", true);
  CHECK(TileDBUtils::is_file(vcf_file));
  CHECK(TileDBUtils::is_file(vcf_file+".tbi"));

  delete gdb;
}


TEST_CASE("api generate_vcf with json multiple threads", "[query_generate_with_json_multiple_threads]") {
  // Define a lambda expression
  auto test_genomicsdb_fn = [](int i) {
    TempDir temp_dir;
    GenomicsDB* gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json);
    const std::string vcf_file = temp_dir.append(std::to_string(i)+".vcf.gz");
    gdb->generate_vcf(vcf_file, "z", true);
    CHECK(TileDBUtils::is_file(vcf_file));
    CHECK(TileDBUtils::is_file(vcf_file+".tbi"));
    delete gdb;
  };

  size_t num_threads = 8;
  std::vector<std::thread> threads;
  for (auto i=0ul; i<num_threads; i++) {
    std::thread thread_object(test_genomicsdb_fn, i);
    threads.push_back(std::move(thread_object));
  }

  CHECK(num_threads == threads.size());

  for (auto i=0ul; i<num_threads; i++) {
    threads[i].join();
  }
}

/**
  Call processor used for testing genomic field annotations
 */
class VariantAnnotationCallProcessor : public GenomicsDBVariantCallProcessor {
  void process(const std::string& sample_name,
               const int64_t* coordinates,
               const genomic_interval_t& genomic_interval,
               const std::vector<genomic_field_t>& genomic_fields) {

     if(sample_name == "HG00141" && coordinates[1] == 12140) {
       CHECK(genomic_fields.size() == 3);
     } else if (sample_name == "HG01958" && coordinates[1] == 12144) {
       CHECK(genomic_fields.size() == 3);
     } else if (sample_name == "HG00141" && coordinates[1] == 17384) {
       CHECK(genomic_interval.contig_name == "1");
       CHECK(genomic_interval.interval.first == 17385);
       CHECK(genomic_interval.interval.second == 17385);
       CHECK(genomic_fields.size() == 9);

       int countExpectedGenomicFields=0;
       CHECK(genomic_fields.size() == 9);
       for (auto i=0u; i<genomic_fields.size(); i++) {
         if (genomic_fields[i].name == "REF") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].str_value() == "G");
         } else if(genomic_fields[i].name == "ALT") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].str_value() == "A|&");
         } else if(genomic_fields[i].name == "GT") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].to_string(get_genomic_field_type("GT")) == "0/1");
         } else if(genomic_fields[i].name == "dataSourceZero_field0") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_int());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "2345");
         } else if(genomic_fields[i].name == "dataSourceZero_field1") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_int());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "5678");
         } else if(genomic_fields[i].name == "dataSourceZero_field2") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_float());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "3.141500");
         } else if(genomic_fields[i].name == "dataSourceZero_field3") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "noot");
         } else if(genomic_fields[i].name == "dataSourceZero_field7") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "");
         } else if(genomic_fields[i].name == "dataSourceZero_ID") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "id001");
         } else {
           countExpectedGenomicFields = -1;
         }
       }
       // Check not just the genomic_fields but that the expected fields were found.
       CHECK(countExpectedGenomicFields == 9);
     } else if (sample_name == "HG01958") {
       CHECK(genomic_interval.contig_name == "1");
       CHECK(genomic_interval.interval.first == 17385);
       CHECK(genomic_interval.interval.second == 17385);
       CHECK(genomic_fields.size() == 9);

       int countExpectedGenomicFields=0;
       for (auto i=0u; i<genomic_fields.size(); i++) {
         if (genomic_fields[i].name == "REF") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].str_value() == "G");
         } else if(genomic_fields[i].name == "ALT") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].str_value() == "T|&");
         } else if(genomic_fields[i].name == "GT") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].to_string(get_genomic_field_type("GT")) == "1/1");
         } else if(genomic_fields[i].name == "DP") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "120");
         } else if(genomic_fields[i].name == "dataSourceZero_field3") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "waldo");
         } else if(genomic_fields[i].name == "dataSourceZero_field4") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_int());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "5679");
         } else if(genomic_fields[i].name == "dataSourceZero_field5") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_float());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "6.660000");
         } else if(genomic_fields[i].name == "dataSourceZero_field6") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "biz");
         } else if(genomic_fields[i].name == "dataSourceZero_ID") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "id002");
         } else {
           countExpectedGenomicFields = -1;
         }
       }

       // Check not just the genomic_fields but that the expected fields were found.
       CHECK(countExpectedGenomicFields == 9);

     } else if (sample_name == "HG01530") {
       CHECK(genomic_interval.contig_name == "1");
       CHECK(genomic_interval.interval.first == 17385);
       CHECK(genomic_interval.interval.second == 17385);
       // For this sample i don't bother with the actual genomic_field values because because they are the same as HG00141 above.
       CHECK(genomic_fields.size() == 10);
     }
  };
};

TEST_CASE("api annotate query_variant_calls with test datasource 0", "[annotate_variant_calls_with_tds0]") {
  using namespace genomicsdb_pb;

  ExportConfiguration *config = new ExportConfiguration();

  config->set_workspace(workspace);
  config->set_callset_mapping_file(callset_mapping);
  config->set_vid_mapping_file(vid_mapping);

  // query_row_ranges
  RowRangeList* row_ranges = config->add_query_row_ranges();
  RowRange* row_range = row_ranges->add_range_list();
  row_range->set_low(0);
  row_range->set_high(3);

  // query_column_ranges
  GenomicsDBColumnOrIntervalList* column_ranges = config->add_query_column_ranges();
  GenomicsDBColumnOrInterval* column_range = column_ranges->add_column_or_interval_list();

  TileDBColumnInterval* tiledb_column_interval = new TileDBColumnInterval();
  tiledb_column_interval->set_begin(0);
  tiledb_column_interval->set_end(1000000000);

  GenomicsDBColumnInterval* column_interval = new GenomicsDBColumnInterval();
  column_interval->set_allocated_tiledb_column_interval(tiledb_column_interval);

  column_range->set_allocated_column_interval(column_interval);

  // query_attributes
  config->add_attributes()->assign("GT");
  config->add_attributes()->assign("DP");

  config->set_reference_genome("inputs/chr1_10MB.fasta.gz");
  config->set_segment_size(40);
  config->set_vcf_header_filename("inputs/template_vcf_header.vcf");

  // Add an annotation dataSource
  AnnotationSource* annotation_source0 = config->add_annotation_source();
  const std::string vcf_file0(ctests_input_dir+"test_datasource0.vcf.bgz");
  const std::string data_source0("dataSourceZero");
  annotation_source0->set_is_vcf(true);
  annotation_source0->set_filename(vcf_file0);
  annotation_source0->set_data_source(data_source0);
  annotation_source0->add_attributes()->assign("field0");
  annotation_source0->add_attributes()->assign("field1");
  annotation_source0->add_attributes()->assign("field2");
  annotation_source0->add_attributes()->assign("field3");
  annotation_source0->add_attributes()->assign("field4");
  annotation_source0->add_attributes()->assign("field5");
  annotation_source0->add_attributes()->assign("field6");
  annotation_source0->add_attributes()->assign("field7");
  annotation_source0->add_attributes()->assign("ID");

  // Add a chromosome specific data-source that should not be loaded.
  // The file is not a valid vcf, so an error will be thrown if the file is loaded,
  // but it won't be because the annotation service knows that the file doesn't have
  // any matches since it only has variants on the 99th chromosome.
  AnnotationSource* annotation_source1 = config->add_annotation_source();
  const std::string vcf_file1(ctests_input_dir+"test_datasource_invalid.vcf.bgz");
  const std::string data_source1("dataSource1");
  annotation_source1->set_is_vcf(true);
  annotation_source1->set_filename(vcf_file1);
  annotation_source1->set_data_source(data_source1);
  annotation_source1->add_attributes()->assign("chrField");
  annotation_source1->add_file_chromosomes()->assign("99");

  std::string config_string;
  CHECK(config->SerializeToString(&config_string));

  config->set_array_name(array);
  CHECK(config->SerializeToString(&config_string));

  GenomicsDB* gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);

  VariantAnnotationCallProcessor variant_annotation_processor;
  gdb->query_variant_calls(variant_annotation_processor, "", GenomicsDB::NONE);
  delete gdb;

  // Check if exception thrown when annotation buffer size is too small to hold the annotated field values
  config->set_annotation_buffer_size(4);
  CHECK(config->SerializeToString(&config_string));
  gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  CHECK_THROWS_AS(gdb->query_variant_calls(variant_annotation_processor, "", GenomicsDB::NONE), GenomicsDBException);
  delete gdb;

  // Add an info field that doesn't exist to the configuration; expect custom exception
  config->set_annotation_buffer_size(10240);
  annotation_source0->add_attributes()->assign("no_exist_field0");
  CHECK(config->SerializeToString(&config_string));
  gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  CHECK_THROWS_AS(gdb->query_variant_calls(variant_annotation_processor, "", GenomicsDB::NONE), GenomicsDBException);
  delete gdb;
}
