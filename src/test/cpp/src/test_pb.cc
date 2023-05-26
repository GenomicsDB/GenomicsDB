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

#include "json_config.h"
#include "variant_query_config.h"
#include "genomicsdb_export_config.pb.h"
#include "vid_mapper_pb.h"

#include "test_common.h"

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
  CHECK(pb_config.enable_shared_posixfs_optimizations() == json_config.enable_shared_posixfs_optimizations());
  CHECK(pb_config.get_query_filter() == json_config.get_query_filter());
  CHECK(pb_config.get_attributes() == json_config.get_attributes());
}

static std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);
static std::string vid_mapping_asa(ctests_input_dir+"vid_protobuf.json");

void check_equal_vid_info_fields(const VidMapper& json, const VidMapper& pb) {
  const FieldInfo *pb_info, *json_info;
  const std::string fields[] = {"PASS", "LowQual", "END", "BaseQRankSum", "ClippingRankSum", "MQRankSum",
        "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP",
        "GT", "AS_RAW_MQ", "AS_RAW_MQRankSum"};
  unsigned int num_fields = pb.get_num_fields();
  CHECK(num_fields == json.get_num_fields());
  for (unsigned i=0; i<sizeof(fields)/sizeof(fields[0]); i++) {
    pb_info = pb.get_field_info(fields[i]);
    json_info = json.get_field_info(fields[i]);
    FieldElementTypeDescriptor pb_fetp = pb_info->get_vcf_type();
    FieldElementTypeDescriptor json_fetp = json_info->get_vcf_type();
    size_t n = pb_fetp.get_num_elements_in_tuple();
    CHECK(n == json_fetp.get_num_elements_in_tuple());
    for(unsigned j=0; j<n; j++) {
      CHECK(pb_fetp.get_tuple_element_bcf_ht_type(j) ==
              json_fetp.get_tuple_element_bcf_ht_type(j));
      CHECK(pb_fetp.get_tuple_element_type_index(j) ==
              json_fetp.get_tuple_element_type_index(j));
    }
    FieldLengthDescriptor pb_fld = pb_info->m_length_descriptor;
    FieldLengthDescriptor json_fld = json_info->m_length_descriptor;
    if (pb_fld.is_fixed_length_field()) {
      CHECK(json_fld.is_fixed_length_field());
      CHECK(pb_fld.get_num_elements() == json_fld.get_num_elements());
    }
    else {
      n = pb_fld.get_num_dimensions();
      CHECK(n == json_fld.get_num_dimensions());
      for(unsigned j=0; j<n; j++) {
        CHECK(pb_fld.get_length_descriptor(j) == json_fld.get_length_descriptor(j));
        CHECK(pb_fld.get_vcf_delimiter(j) == json_fld.get_vcf_delimiter(j));
      }
    }
    CHECK(pb_info->m_VCF_field_combine_operation ==
            json_info->m_VCF_field_combine_operation);

    CHECK(pb_info->m_disable_remap_missing_with_non_ref ==
            json_info->m_disable_remap_missing_with_non_ref);
  }
}

void check_equal_vid_contig(const VidMapper& json, const VidMapper& pb) {
  ContigInfo pb_contig, json_contig;
  const std::string contigs[] = {"1", "2", "3", "4", "5", "6", "7", "8",
        "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
        "20", "21", "22", "X", "Y", "MT"};
  unsigned int num_contigs = pb.get_num_contigs();
  CHECK(num_contigs == json.get_num_contigs());
  for(unsigned i=0; i<sizeof(contigs)/sizeof(contigs[0]); i++) {
    CHECK(pb.get_contig_info(contigs[i], pb_contig));
    CHECK(json.get_contig_info(contigs[i], json_contig));
    CHECK(pb_contig.m_length == json_contig.m_length);
    CHECK(pb_contig.m_tiledb_column_offset == json_contig.m_tiledb_column_offset);
  }
}

void check_equal_vid_config(const VidMapper& json, const VidMapper& pb) {
  check_equal_vid_contig(json, pb);
  check_equal_vid_info_fields(json, pb);
}

TEST_CASE("pb_query_config_test", "[protobuf_config]")
{
  if(!g_query_json_file.empty()) {
    GenomicsDBConfigBase json_config;
    json_config.read_from_file(g_query_json_file, 0);
    genomicsdb_pb::ExportConfiguration export_config;
    GenomicsDBConfigBase::get_pb_from_json_file(&export_config, g_query_json_file);
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
    VariantQueryConfig pb_config3;
    pb_config3.read_from_PB_binary_string(binary_pb_string, 0);
    check_equal_query_config(json_config, pb_config3);
  }
}

TEST_CASE("pb_vid_mapping_test", "[vid_protobuf_config]")
{
  VidMapper vid_json, vid_pb;
  vid_json = std::move(FileBasedVidMapper(vid_mapping_asa));
  VidMappingPB vid_config;
  GenomicsDBConfigBase::get_pb_from_json_file(&vid_config, vid_mapping_asa);
  CHECK(vid_config.IsInitialized());
  vid_pb.parse_vidmap_protobuf(&vid_config);
  // check against pb generated from file
  check_equal_vid_config(vid_json, vid_pb);
}

static std::string vid_mapping_overlapping_contigs(ctests_input_dir+"vid_protobuf_overlapping_contigs.json");

TEST_CASE("pb_overlapping_vid_mapping_test", "[vid_protobuf_overlapping_contigs]")
{
  VidMapper vid_pb;
  VidMappingPB vid_config;
  GenomicsDBConfigBase::get_pb_from_json_file(&vid_config, vid_mapping_overlapping_contigs);
  CHECK(vid_config.IsInitialized());
  CHECK_THROWS_AS(vid_pb.parse_vidmap_protobuf(&vid_config), ProtoBufBasedVidMapperException);
}
