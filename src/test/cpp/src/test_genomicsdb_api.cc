/**
 * src/test/cpp/src/test_genomicsdb_api.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2019-2021 Omics Data Automation, Inc.
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

TEST_CASE("api get_version", "[get_version]") {
  REQUIRE(genomicsdb_version().size() > 0);
  REQUIRE_THAT(genomicsdb_version(), Catch::Equals(GENOMICSDB_VERSION));
}

TEST_CASE("api empty_args", "[empty_args]") {
  std::string empty_string;

  // Constructor that allows workspace/callset/vidmapping/reference genome specification
  CHECK_THROWS_AS(new GenomicsDB(empty_string, empty_string, empty_string, empty_string), std::exception);
  CHECK_THROWS_AS(new GenomicsDB("ws", empty_string, empty_string, empty_string), std::exception);
  CHECK_THROWS_AS(new GenomicsDB("ws", "callset", empty_string, empty_string), std::exception);
  CHECK_THROWS_AS(new GenomicsDB("ws", "callset", "vid", empty_string), std::exception);

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
   CHECK(variant_calls.get_genomic_field_type(call_genomic_fields[1].name).is_string());
   CHECK(call_genomic_fields[1].str_value() == "&");
}

TEST_CASE("api query_variants direct DP", "[query_variants]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, {"DP"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}, {{0,3}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array));
  delete gdb;
}

TEST_CASE("api query_variants direct DP and GT", "[query_variants_DP_GT]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, {"DP", "GT"}, 40);
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}, {{0,3}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array, {{0,1000000000}}));
  check_query_variants_results(gdb, array, gdb->query_variants(array));
  delete gdb;
}

TEST_CASE("api query_variants direct DP and GT with PP", "[query_variants_DP_GT_with_PP]") {
  GenomicsDB* gdb = new GenomicsDB(workspace_PP, callset_mapping, vid_mapping_PP, reference_genome, {"DP", "GT"}, 40);
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

  void process(const interval_t& interval) {
    m_num_query_intervals++;
    CHECK(interval.first == 0);
    CHECK(interval.second == 1000000000);
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
      int found_fields = 0;
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

  int m_attributes_size;
  int m_num_query_intervals = 0;
  int m_processed_rows = 0;
  int m_sample_found = false;
  int m_is_PP;
};

class TwoQueryIntervalsProcessor : public GenomicsDBVariantCallProcessor {
 public:
  TwoQueryIntervalsProcessor() {};

  ~TwoQueryIntervalsProcessor() {
    CHECK(m_num_query_intervals == 2);
    CHECK(m_processed_rows == 5);
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

  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, attributes, 40);

  gdb->query_variant_calls(array);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor(attributes);
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}}, {{0,3}});
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}});
  gdb->query_variant_calls(one_query_interval_processor, array);

  TwoQueryIntervalsProcessor two_query_intervals_processor;
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}}, {{0,3}});
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}});

  delete gdb;
}

TEST_CASE("api query_variant_calls direct DP and GT", "[query_variant_calls_direct_DP_GT]") {
  const std::vector<std::string> attributes = {"GT", "DP"};

  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, attributes, 40);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor(attributes);
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}}, {{0,3}});
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}});
  gdb->query_variant_calls(one_query_interval_processor, array);

  TwoQueryIntervalsProcessor two_query_intervals_processor;
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}}, {{0,3}});
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}});

  delete gdb;
}

TEST_CASE("api query_variant_calls direct DP and GT with PP", "[query_variant_calls_direct_DP_GT_with_PP]") {
  const std::vector<std::string> attributes = {"GT", "DP"};

  GenomicsDB* gdb = new GenomicsDB(workspace_PP, callset_mapping, vid_mapping_PP, reference_genome, attributes, 40);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor(attributes, true);
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}}, {{0,3}});
  gdb->query_variant_calls(one_query_interval_processor, array, {{0,1000000000}});
  gdb->query_variant_calls(one_query_interval_processor, array);

  TwoQueryIntervalsProcessor two_query_intervals_processor;
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}}, {{0,3}});
  gdb->query_variant_calls(two_query_intervals_processor, array, {{0,17000},{17000,18000}});

  delete gdb;
}

TEST_CASE("api query_variant_calls with json", "[query_variant_calls_with_json]") {
  GenomicsDB* gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json);
  gdb->query_variant_calls();

  OneQueryIntervalProcessor one_query_interval_processor;
  gdb->query_variant_calls(one_query_interval_processor);
  delete gdb;

  OneQueryIntervalProcessor another_query_interval_processor;
  gdb = new GenomicsDB(query_json, GenomicsDB::JSON_FILE, loader_json, 0);
  gdb->query_variant_calls(another_query_interval_processor);
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

  // no array name set, should throw exception
  CHECK_THROWS_AS(new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0), GenomicsDBConfigException);

  config->set_array_name(array);
  CHECK(config->SerializeToString(&config_string));

  GenomicsDB* gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  gdb->query_variant_calls();
  OneQueryIntervalProcessor one_query_interval_processor;
  gdb->query_variant_calls(one_query_interval_processor);
  delete gdb;

  // try query with contig intervals instead of tiledb column intervals
  ContigInterval* contig_interval = new ContigInterval();
  contig_interval->set_contig("1");
  contig_interval->set_begin(1);
  contig_interval->set_end(249250621);
  column_interval->Clear();
  column_interval->set_allocated_contig_interval(contig_interval);
  CHECK(config->SerializeToString(&config_string));
  gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  gdb->query_variant_calls();
  delete gdb;
}

TEST_CASE("api generate_vcf direct", "[query_generate_vcf_direct]") {
  TempDir temp_dir;
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, {"DP"}, 40);

  const std::string vcf_file1 = temp_dir.append("1.vcf.gz");
  gdb->generate_vcf(array, {{0,1000000000}}, {{0,3}}, vcf_file1);
  CHECK(TileDBUtils::is_file(vcf_file1));
  CHECK(!TileDBUtils::is_file(vcf_file1+".tbi"));

  gdb->generate_vcf(array, {{0,1000000000}}, {{0,3}}, vcf_file1, "z", true);
  CHECK(TileDBUtils::is_file(vcf_file1));
  CHECK(TileDBUtils::is_file(vcf_file1+".tbi"));

  CHECK_THROWS_AS(gdb->generate_vcf(array, {{0,1000000000}}, {{0,3}}, vcf_file1, "z", false), std::exception);

  const std::string vcf_file2 = temp_dir.append("3.vcf.gz");
  gdb->generate_vcf(array, {{0,13000}, {13000, 1000000000}}, {{0,3}}, vcf_file2, "z", true);
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

  int num_threads = 8;
  std::vector<std::thread> threads;
  for (auto i=0; i<num_threads; i++) {
    std::thread thread_object(test_genomicsdb_fn, i);
    threads.push_back(std::move(thread_object));
  }

  CHECK(num_threads == threads.size());

  for (auto i=0; i<num_threads; i++) {
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
           CHECK(genomic_fields[i].str_value() == "");
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
           CHECK(genomic_fields[i].str_value() == "noot");
         } else if(genomic_fields[i].name == "dataSourceZero_field7") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
           CHECK(genomic_fields[i].to_string(get_genomic_field_type(genomic_fields[i].name)) == "");
         } else if(genomic_fields[i].name == "dataSourceZero_ID") {
           ++countExpectedGenomicFields;
           CHECK(get_genomic_field_type(genomic_fields[i].name).is_string());
           CHECK(genomic_fields[i].str_value() == "id001");
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
           CHECK(genomic_fields[i].str_value().length() == 1);
         } else if(genomic_fields[i].name == "DP") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].str_value() == "x");
         } else if(genomic_fields[i].name == "dataSourceZero_field3") {
           ++countExpectedGenomicFields;
           CHECK(genomic_fields[i].str_value() == "waldo");
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
           CHECK(genomic_fields[i].str_value() == "id002");
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
  // The file is not a valid vcf, so if an error will be thrown if the file is loaded,
  // but it won't be because the annotation service knows that the file doesn't have
  // any matches since it only has chromosome 7 variants.
  AnnotationSource* annotation_source7 = config->add_annotation_source();
  const std::string vcf_file7(ctests_input_dir+"test_datasource_chr7.vcf.bgz");
  const std::string data_source7("dataSourceChr7");
  annotation_source7->set_is_vcf(true);
  annotation_source7->set_filename(vcf_file7);
  annotation_source7->set_data_source(data_source7);
  annotation_source7->add_attributes()->assign("chrField");
  annotation_source7->add_file_chromosomes()->assign("7");

  std::string config_string;
  CHECK(config->SerializeToString(&config_string));

  config->set_array_name(array);
  CHECK(config->SerializeToString(&config_string));

  GenomicsDB* gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);

  VariantAnnotationCallProcessor variant_annotation_processor;
  gdb->query_variant_calls(variant_annotation_processor);
  delete gdb;

  // Add an info field that doesn't exist to the configuration; expect custom exception
  annotation_source0->add_attributes()->assign("no_exist_field0");
  CHECK(config->SerializeToString(&config_string));
  config->set_array_name(array);
  CHECK(config->SerializeToString(&config_string));
  gdb = new GenomicsDB(config_string, GenomicsDB::PROTOBUF_BINARY_STRING, loader_json, 0);
  CHECK_THROWS_AS(gdb->query_variant_calls(variant_annotation_processor), GenomicsDBException);
  delete gdb;
}
