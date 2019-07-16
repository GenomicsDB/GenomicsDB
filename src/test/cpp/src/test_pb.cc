/**
 * The MIT License (MIT)
 * Copyright (c) 2019 Omics Data Automation, Inc.
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

#include <catch2/catch.hpp>
#include <google/protobuf/util/json_util.h>
#include <google/protobuf/util/type_resolver.h>
#include <google/protobuf/util/type_resolver_util.h>

#include "tiledb_utils.h"
#include "json_config.h"
#include "variant_query_config.h"
#include "genomicsdb_export_config.pb.h"

extern std::string g_pb_query_json_file;

void check_equal_query_config(const GenomicsDBConfigBase& json_config, const GenomicsDBConfigBase& pb_config) {
  CHECK(pb_config.get_workspace(0) == json_config.get_workspace(0));
  CHECK(pb_config.get_array_name(0) == json_config.get_array_name(0));
  CHECK(pb_config.get_query_column_ranges(0) == json_config.get_query_column_ranges(0));
  CHECK(pb_config.get_query_row_ranges(0) == json_config.get_query_row_ranges(0));
  CHECK(pb_config.get_segment_size() == json_config.get_segment_size());
  CHECK(pb_config.get_determine_sites_with_max_alleles() == json_config.get_determine_sites_with_max_alleles());
  CHECK(pb_config.get_max_diploid_alt_alleles_that_can_be_genotyped()
      == json_config.get_max_diploid_alt_alleles_that_can_be_genotyped());
  CHECK(pb_config.get_combined_vcf_records_buffer_size_limit()
      == json_config.get_combined_vcf_records_buffer_size_limit());
  CHECK(pb_config.get_vcf_header_filename() == json_config.get_vcf_header_filename());
  CHECK(pb_config.get_vcf_output_format() == json_config.get_vcf_output_format());
  CHECK(pb_config.get_vcf_output_filename() == json_config.get_vcf_output_filename());
  CHECK(pb_config.get_reference_genome() == json_config.get_reference_genome());
  CHECK(pb_config.produce_GT_field() == json_config.produce_GT_field());
  CHECK(pb_config.produce_FILTER_field() == json_config.produce_FILTER_field());
  CHECK(pb_config.sites_only_query() == json_config.sites_only_query());
  CHECK(pb_config.index_output_VCF() == json_config.index_output_VCF());
  CHECK(pb_config.produce_GT_with_min_PL_value_for_spanning_deletions()
      == json_config.produce_GT_with_min_PL_value_for_spanning_deletions());
  CHECK(pb_config.disable_file_locking_in_tiledb() == json_config.disable_file_locking_in_tiledb());
  CHECK(pb_config.get_query_filter() == json_config.get_query_filter());
  CHECK(pb_config.get_attributes() == json_config.get_attributes());
}

TEST_CASE("pb_query_config_test", "[protobuf_config]")
{
  if(!g_pb_query_json_file.empty()) {
    GenomicsDBConfigBase json_config;
    json_config.read_from_file(g_pb_query_json_file, 0);
    char *json_buffer = 0;
    size_t json_buffer_length;
    if (TileDBUtils::read_entire_file(g_pb_query_json_file, (void **)&json_buffer, &json_buffer_length) != TILEDB_OK
	|| !json_buffer || json_buffer_length == 0) { 
      free(json_buffer);
      throw GenomicsDBConfigException((std::string("Could not open query JSON file \"")+g_pb_query_json_file+"\"").c_str());
    }
    genomicsdb_pb::ExportConfiguration export_config;
    //The function JsonStringToMessage was made available in Protobuf version 3.0.0. However,
    //to maintain compatibility with GATK-4, we need to use 3.0.0-beta-1. This version doesn't have
    //the JsonStringToMessage method. A workaround is as follows.
    //https://stackoverflow.com/questions/41651271/is-there-an-example-of-protobuf-with-text-output
    {
      std::string json_to_binary_output;
      google::protobuf::util::TypeResolver* resolver= google::protobuf::util::NewTypeResolverForDescriptorPool(
	  "", google::protobuf::DescriptorPool::generated_pool());
      auto status = google::protobuf::util::JsonToBinaryString(resolver,
	  "/"+export_config.GetDescriptor()->full_name(), json_buffer,
	  &json_to_binary_output);
      CHECK(status.ok());
      delete resolver;
      auto success = export_config.ParseFromString(json_to_binary_output);
      CHECK(success);
    }
    //auto status = google::protobuf::util::JsonStringToMessage(json_buffer, &export_config);
    //CHECK(status.ok());
    free(json_buffer);
    json_buffer = 0;
    CHECK(export_config.IsInitialized());
    GenomicsDBConfigBase pb_config;
    pb_config.read_from_PB(&export_config, 0);
    check_equal_query_config(json_config, pb_config);
    std::string binary_pb_string;
    //This could be the method to pass information from Python/R/Java to C++
    auto serialize_success = export_config.SerializeToString(&binary_pb_string);
    CHECK(serialize_success);
    genomicsdb_pb::ExportConfiguration deserialize_export_config;
    auto deserialize_success = deserialize_export_config.ParseFromString(binary_pb_string);
    CHECK(deserialize_success);
    CHECK(deserialize_export_config.IsInitialized());
    VariantQueryConfig pb_config2;
    pb_config2.read_from_PB(&deserialize_export_config, 0);
    check_equal_query_config(json_config, pb_config2);
  }
}
