/**
 * src/test/cpp/src/test_config.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation, Inc.
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
 * Test GenomicsDBConfig/VariantQueryConfig/... classes
 *
 */

#include <catch2/catch.hpp>

#include "genomicsdb_config_base.h"
#include "tiledb_utils.h"
#include "variant_query_config.h"

static std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

static std::string query_json(ctests_input_dir+"query.json");
static std::string loader_json(ctests_input_dir+"loader.json");

TEST_CASE("Check configuration with json file", "[basic_config_check]") {
  GenomicsDBImportConfig loader_config;
  loader_config.read_from_file(loader_json);
  
  VariantQueryConfig query_config_from_file;
  query_config_from_file.update_from_loader(loader_config);
  query_config_from_file.read_from_file(query_json);
  query_config_from_file.subset_query_column_ranges_based_on_partition(loader_config);

  CHECK(query_config_from_file.get_workspace(0) == "inputs/ws");
  CHECK(query_config_from_file.get_array_name(0) == "t0_1_2");
  CHECK(query_config_from_file.get_query_column_ranges(0).size() ==  1);
  CHECK(query_config_from_file.get_query_column_ranges(0)[0].first ==  0);
  CHECK(query_config_from_file.get_query_column_ranges(0)[0].second ==  1000000000);
  CHECK(query_config_from_file.get_query_row_ranges(0).size() ==  1);
  CHECK(query_config_from_file.get_query_row_ranges(0)[0].first ==  0);
  CHECK(query_config_from_file.get_query_row_ranges(0)[0].second ==  3);
  CHECK(query_config_from_file.get_query_filter().empty());
  CHECK(query_config_from_file.get_segment_size() == 40);
  CHECK(query_config_from_file.get_determine_sites_with_max_alleles() == 0);
  CHECK(query_config_from_file.get_vcf_output_filename() == "-");
  CHECK(query_config_from_file.get_reference_genome() ==  "inputs/chr1_10MB.fasta.gz");
  CHECK(query_config_from_file.produce_GT_field() == false);
  CHECK(query_config_from_file.produce_FILTER_field() ==  false);
  CHECK(query_config_from_file.sites_only_query() == false);
  CHECK(query_config_from_file.index_output_VCF() ==  false);
  CHECK(query_config_from_file.produce_GT_with_min_PL_value_for_spanning_deletions() == false);
  CHECK(query_config_from_file.get_vid_mapper().is_initialized() ==  true);
  CHECK(query_config_from_file.get_vid_mapper().is_callset_mapping_initialized() == true);
}

TEST_CASE("Compare configuration with json file and string", "[compare_json_types]") {
  GenomicsDBImportConfig loader_config;
  loader_config.read_from_file(loader_json);
  
  VariantQueryConfig query_config_from_file;
  query_config_from_file.update_from_loader(loader_config);
  query_config_from_file.read_from_file(query_json);
  query_config_from_file.subset_query_column_ranges_based_on_partition(loader_config);

  char *query_json_buffer=NULL, *loader_json_buffer=NULL;
  size_t query_json_buffer_length, loader_json_buffer_length;
  CHECK(TileDBUtils::read_entire_file(query_json, (void **)&query_json_buffer, &query_json_buffer_length) == TILEDB_OK);
  CHECK(query_json_buffer);
  CHECK(query_json_buffer_length > 0);

  VariantQueryConfig query_config_from_str;
  query_config_from_str.update_from_loader(loader_config);
  query_config_from_str.read_from_JSON_string(query_json_buffer);
  query_config_from_str.subset_query_column_ranges_based_on_partition(loader_config);
  
  // Compare configs read as JSON file or JSON string - they should be equivalent
  CHECK(query_config_from_file.get_workspace(0) == query_config_from_str.get_workspace(0));
  CHECK(query_config_from_file.get_array_name(0) == query_config_from_str.get_array_name(0));
  CHECK(query_config_from_file.get_column_partition(0) == query_config_from_str.get_column_partition(0));
  CHECK(query_config_from_file.get_row_partition(0) == query_config_from_str.get_row_partition(0));
  CHECK(query_config_from_file.get_query_column_ranges(0).size() ==  query_config_from_str.get_query_column_ranges(0).size());
  CHECK(query_config_from_file.get_query_column_ranges(0)[0] == query_config_from_str.get_query_column_ranges(0)[0]);
  CHECK(query_config_from_file.get_query_row_ranges(0).size() ==  query_config_from_str.get_query_row_ranges(0).size());
  CHECK(query_config_from_file.get_query_row_ranges(0)[0] == query_config_from_str.get_query_row_ranges(0)[0]);
  CHECK(query_config_from_file.get_query_filter() == query_config_from_str.get_query_filter());
  CHECK(query_config_from_file.get_segment_size() == query_config_from_str.get_segment_size());
  CHECK(query_config_from_file.get_determine_sites_with_max_alleles() == query_config_from_str.get_determine_sites_with_max_alleles());
  CHECK(query_config_from_file.get_max_diploid_alt_alleles_that_can_be_genotyped() == query_config_from_str.get_max_diploid_alt_alleles_that_can_be_genotyped());
  CHECK(query_config_from_file.get_max_genotype_count() == query_config_from_str.get_max_genotype_count());
  CHECK(query_config_from_file.get_combined_vcf_records_buffer_size_limit() ==  query_config_from_str.get_combined_vcf_records_buffer_size_limit());
  CHECK(query_config_from_file.get_vcf_output_filename() == query_config_from_str.get_vcf_output_filename());
  CHECK(query_config_from_file.get_reference_genome() ==  query_config_from_str.get_reference_genome());
  CHECK(query_config_from_file.produce_GT_field() == query_config_from_str.produce_GT_field());
  CHECK(query_config_from_file.produce_FILTER_field() ==  query_config_from_str.produce_FILTER_field());
  CHECK(query_config_from_file.sites_only_query() == query_config_from_str.sites_only_query());
  CHECK(query_config_from_file.index_output_VCF() == query_config_from_str.index_output_VCF());
  CHECK(query_config_from_file.produce_GT_with_min_PL_value_for_spanning_deletions() == query_config_from_str.produce_GT_with_min_PL_value_for_spanning_deletions());
  CHECK(query_config_from_file.get_vid_mapper().is_initialized() == query_config_from_str.get_vid_mapper().is_initialized());
  CHECK(query_config_from_file.get_vid_mapper().is_callset_mapping_initialized() == query_config_from_str.get_vid_mapper().is_callset_mapping_initialized());
}

