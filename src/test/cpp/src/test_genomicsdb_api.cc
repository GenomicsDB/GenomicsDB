/**
 * src/test/cpp/src/test_genomicsdb_api.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2019-2020 Omics Data Automation, Inc.
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

static std::string array("t0_1_2");

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
   CHECK(call_genomic_fields.size() == 2);
   CHECK(call_genomic_fields[0].name == "REF");
   CHECK(variant_calls.get_genomic_field_type(call_genomic_fields[0].name).is_string());
   CHECK(call_genomic_fields[0].str_value() == "C");
   CHECK(call_genomic_fields[1].name == "ALT");
   CHECK(variant_calls.get_genomic_field_type(call_genomic_fields[1].name).is_string());
   CHECK(call_genomic_fields[1].str_value() == "&");
}

TEST_CASE("api query_variants direct", "[query_variants]") {
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, {"DP"}, 40);
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
  char *query_json_buffer=NULL, *loader_json_buffer=NULL;
  size_t query_json_buffer_length, loader_json_buffer_length;
  CHECK(TileDBUtils::read_entire_file(query_json, (void **)&query_json_buffer, &query_json_buffer_length) == TILEDB_OK);
  CHECK(query_json_buffer);
  CHECK(query_json_buffer_length > 0);

  GenomicsDB* gdb = new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json);
  check_query_variants_results(gdb, array, gdb->query_variants());
  delete gdb;

  gdb = new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json, 0);
  check_query_variants_results(gdb, array, gdb->query_variants());
  delete gdb;
  
  CHECK_THROWS_AS(new GenomicsDB(query_json_buffer, GenomicsDB::JSON_STRING, loader_json, 1), GenomicsDBConfigException);
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
  OneQueryIntervalProcessor() {};

  ~OneQueryIntervalProcessor() {
    CHECK(m_num_query_intervals == 1);
    CHECK(m_processed_rows == 5);
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
    m_processed_rows++;
    auto row = coordinates[0];
    CHECK(row <= 2);
    if (row == 2) {
      CHECK(sample_name == "HG01530");
      CHECK(genomic_interval.contig_name == "1");
      CHECK(genomic_interval.interval.first == 17385);
      CHECK(genomic_interval.interval.second == 17385);
      CHECK(genomic_fields.size() == 3);
      CHECK(genomic_fields[0].name == "REF");
      CHECK(get_genomic_field_type(genomic_fields[0].name).is_string());
      CHECK(genomic_fields[0].str_value() == "G");
      CHECK(genomic_fields[1].name == "ALT");
      CHECK(get_genomic_field_type(genomic_fields[1].name).is_string());
      CHECK(genomic_fields[1].str_value() == "A|&");
      CHECK(genomic_fields[1].to_string(get_genomic_field_type(genomic_fields[1].name)) == "[A, <NON_REF>]");
      CHECK(genomic_fields[2].name == "DP");
      CHECK(get_genomic_field_type(genomic_fields[2].name).is_int());
      CHECK(genomic_fields[2].int_value_at(0) == 76);
      CHECK(genomic_fields[2].to_string(get_genomic_field_type(genomic_fields[2].name)) == "76");
    }
  };

  int m_num_query_intervals = 0;
  int m_processed_rows = 0;
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
  GenomicsDB* gdb = new GenomicsDB(workspace, callset_mapping, vid_mapping, reference_genome, {"DP"}, 40);

  // Default query variant call without processor, should just print the calls
  gdb->query_variant_calls(array, {{0,1000000000}}, {{0,3}});

  NullVariantCallProcessor null_processor;
  gdb->query_variant_calls(null_processor, array, {{0,1000000000}}, {{0,3}});

  OneQueryIntervalProcessor one_query_interval_processor;
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
