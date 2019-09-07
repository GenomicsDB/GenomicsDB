/**
 * The MIT License (MIT)
 * Copyright (c) 2019 Omics Data Automation, Inc
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

//Enable asserts
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <zlib.h>
#include <google/protobuf/util/json_util.h>
#include <google/protobuf/util/type_resolver.h>
#include <google/protobuf/util/type_resolver_util.h>
#include "tiledb_utils.h"
#include "json_config.h"
#include "genomicsdb_config_base.h"
#include "genomicsdb_export_config.pb.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBConfigException(#X);

void GenomicsDBConfigBase::read_from_PB_binary_string(const std::string& str, const int rank) {
  genomicsdb_pb::ExportConfiguration export_config;
  auto status = export_config.ParseFromString(str);
  if(!status || !export_config.IsInitialized())
    throw GenomicsDBConfigException("Could not deserialize PB binary string");
  read_from_PB(&export_config, rank);
}

void GenomicsDBConfigBase::read_from_PB(const genomicsdb_pb::ExportConfiguration* export_config, const int rank) {
  VERIFY_OR_THROW(export_config && export_config->IsInitialized());
  read_and_initialize_vid_and_callset_mapping_if_available(export_config);
  VERIFY_OR_THROW(m_vid_mapper.is_initialized() && m_vid_mapper.is_callset_mapping_initialized());
  m_workspaces.clear();
  m_workspaces.emplace_back(export_config->workspace());
  m_single_workspace_path = true;
  m_array_names.clear();
  if(!export_config->has_array_name())
    throw GenomicsDBConfigException("PB export config object must have array_name set");
  m_array_names.emplace_back(export_config->array_name());
  m_single_array_name = true;
  if(export_config->has_scan_full() && export_config->scan_full())
    scan_whole_array();
  m_attributes.resize(export_config->attributes_size());
  for(auto i=0;i<export_config->attributes_size();++i)
    m_attributes[i] = std::move(std::string(export_config->attributes(i)));
  if(export_config->has_query_filter())
    m_query_filter = export_config->query_filter();
  if(export_config->has_segment_size())
    m_segment_size = export_config->segment_size();
  if(export_config->has_vcf_header_filename())
    m_vcf_header_filename = export_config->vcf_header_filename();
  m_vcf_output_filename = export_config->has_vcf_output_filename()
    ? export_config->vcf_output_filename() : "-"; //else stdout
  m_vcf_output_format = export_config->has_vcf_output_format()
    ? export_config->vcf_output_format() : "";
  if(!m_scan_whole_array
      && export_config->query_column_ranges_size() == 0
      && export_config->query_contig_intervals_size() == 0
      && export_config->query_row_ranges_size() == 0
      && export_config->query_sample_names_size() == 0)
    throw GenomicsDBConfigException("Must have one of \"scan_full\" or \"query_column_ranges\" or \
	\"query_contig_intervals\" or \"query_row_ranges\" \
	or \"query_sample_names\" in the Protobuf query config object");
  if(!m_scan_whole_array) {
    if(export_config->query_column_ranges_size() > 0 && export_config->query_contig_intervals_size() > 0)
      throw GenomicsDBConfigException("Protobuf query object cannot have both query_column_ranges and query_contig_intervals at the same time");
    if(export_config->query_column_ranges_size() > 0) {
      if(export_config->query_column_ranges_size() == 1)
	m_single_query_column_ranges_vector = true;
      m_column_ranges.resize(export_config->query_column_ranges_size());
      for(auto i=0;i<export_config->query_column_ranges_size();++i) {
	const auto& inner_list_wrapper = export_config->query_column_ranges(i);
	m_column_ranges[i].resize(inner_list_wrapper.column_or_interval_list_size());
	for(auto j=0;j<inner_list_wrapper.column_or_interval_list_size();++j) {
	  const auto& gdb_column_or_interval = inner_list_wrapper.column_or_interval_list(j);
	  if(gdb_column_or_interval.has_column()) {
	    const auto& gdb_column = gdb_column_or_interval.column();
	    if(gdb_column.has_tiledb_column()) {
	      m_column_ranges[i][j].first = gdb_column.tiledb_column();
	      m_column_ranges[i][j].second = gdb_column.tiledb_column();
	    }
	    else {
	      if(!gdb_column.has_contig_position())
		throw GenomicsDBConfigException("Protobuf object GenomicsDBColumn must have either tiledb_column \
		    or contig_position defined");
	      const auto& contig_position = gdb_column.contig_position();
	      assert(m_vid_mapper.is_initialized());
	      ContigInfo contig_info;
	      if (!m_vid_mapper.get_contig_info(contig_position.contig(), contig_info))
		throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : "
		    + contig_position.contig());
	      m_column_ranges[i][j] = verify_contig_position_and_get_tiledb_column_interval(contig_info,
		  contig_position.position(), contig_position.position());
	    }
	  }
	  else {
	    if(!gdb_column_or_interval.has_column_interval())
	      throw GenomicsDBConfigException("Protobuf object GenomicsDBColumnOrInterval must have either column \
		  or column_interval set");
	    const auto& gdb_column_interval = gdb_column_or_interval.column_interval();
	    if(gdb_column_interval.has_tiledb_column_interval()) {
	      m_column_ranges[i][j].first = gdb_column_interval.tiledb_column_interval().begin();
	      m_column_ranges[i][j].second = gdb_column_interval.tiledb_column_interval().end();
	    }
	    else {
	      if(!gdb_column_interval.has_contig_interval())
		throw GenomicsDBConfigException("GenomicsDBColumnInterval Protobuf object must have tiledb_column_interval\
		    or contig_interval");
	      const auto& contig_interval = gdb_column_interval.contig_interval();
	      assert(m_vid_mapper.is_initialized());
	      ContigInfo contig_info;
	      if (!m_vid_mapper.get_contig_info(contig_interval.contig(), contig_info))
		throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : "
		    + contig_interval.contig());
	      m_column_ranges[i][j] = verify_contig_position_and_get_tiledb_column_interval(contig_info,
		  contig_interval.begin(), contig_interval.end());
	    }
	  }
	  if(m_column_ranges[i][j].first > m_column_ranges[i][j].second)
	    std::swap(m_column_ranges[i][j].first, m_column_ranges[i][j].second);
	}
      }
    }
    if(export_config->query_contig_intervals_size() > 0) {
      m_single_query_column_ranges_vector = true;
      m_column_ranges.resize(1u);
      m_column_ranges[0].resize(export_config->query_contig_intervals_size());
      for(auto i=0;i<export_config->query_contig_intervals_size();++i) {
	const auto& contig_interval = export_config->query_contig_intervals(i);
	ContigInfo contig_info;
	if (!m_vid_mapper.get_contig_info(contig_interval.contig(), contig_info))
	  throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : "
	      + contig_interval.contig());
	if(!contig_interval.has_begin() && contig_interval.has_end())
	  throw GenomicsDBConfigException("ContigInterval must have begin if end is specified");
	auto begin = contig_interval.has_begin() ? contig_interval.begin() : 1ll;
	auto end = contig_interval.has_end() ? contig_interval.end()
	  : (contig_interval.has_begin() ? contig_interval.begin() : contig_info.m_length);
	m_column_ranges[0][i] = verify_contig_position_and_get_tiledb_column_interval(contig_info, begin, end);
      }
    }

    if(export_config->query_row_ranges_size() > 0 && export_config->query_sample_names_size() > 0)
      throw GenomicsDBConfigException("Protobuf query object cannot have both query_row_ranges and query_sample_names at the same time");
    if(export_config->query_sample_names_size() > 0) {
      m_single_query_row_ranges_vector = true;
      m_row_ranges.resize(1u);
      m_row_ranges[0].resize(export_config->query_sample_names_size());
      for(auto i=0;i<export_config->query_sample_names_size();++i) {
	int64_t row_idx = -1;
	auto status = m_vid_mapper.get_tiledb_row_idx(row_idx, export_config->query_sample_names(i));
	if(!status)
	  throw GenomicsDBConfigException(std::string("Unknown sample name ") + export_config->query_sample_names(i));
	m_row_ranges[0][i].first = row_idx;
	m_row_ranges[0][i].second = row_idx;
      }
    }
    if(export_config->query_row_ranges_size() > 0) {
      m_row_ranges.resize(export_config->query_row_ranges_size());
      m_single_query_row_ranges_vector = (m_row_ranges.size() == 1u);
      for(auto i=0;i<export_config->query_row_ranges_size();++i) {
	auto& row_range_vec = m_row_ranges[i];
	const auto& pb_row_range_list = export_config->query_row_ranges(i);
	row_range_vec.resize(pb_row_range_list.range_list_size());
	for(auto j=0;j<pb_row_range_list.range_list_size();++j) {
	  row_range_vec[j].first = pb_row_range_list.range_list(j).low();
	  row_range_vec[j].second = pb_row_range_list.range_list(j).high();
	  if(row_range_vec[j].second < row_range_vec[j].first)
	    std::swap(row_range_vec[j].first, row_range_vec[j].second);
	}
      }
    }
  }
  if(export_config->has_reference_genome())
    m_reference_genome = export_config->reference_genome();
  //Limit on max #alt alleles so that PL fields get re-computed
  m_max_diploid_alt_alleles_that_can_be_genotyped = export_config->has_max_diploid_alt_alleles_that_can_be_genotyped()
    ? export_config->max_diploid_alt_alleles_that_can_be_genotyped() : MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;
  m_max_genotype_count = export_config->has_max_genotype_count()
    ? export_config->max_genotype_count() : MAX_GENOTYPE_COUNT;
  m_combined_vcf_records_buffer_size_limit = export_config->has_combined_vcf_records_buffer_size_limit()
    ? export_config->combined_vcf_records_buffer_size_limit() : DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE;
  //Cannot be 0
  m_combined_vcf_records_buffer_size_limit = std::max<size_t>(1ull, m_combined_vcf_records_buffer_size_limit);
  //GATK CombineGVCF does not produce GT field by default - option to produce GT
  m_produce_GT_field = export_config->has_produce_gt_field() ? export_config->produce_gt_field() : false; 
  //GATK CombineGVCF does not produce FILTER field by default - option to produce FILTER
  m_produce_FILTER_field = export_config->has_produce_filter_field() ? export_config->produce_filter_field() : false; 
  //index output VCF file
  m_index_output_VCF = export_config->has_index_output_vcf() ? export_config->index_output_vcf() : false;
  //sites-only query - doesn't produce any of the FORMAT fields
  m_sites_only_query = export_config->has_sites_only_query() && export_config->sites_only_query();
  //when producing GT, use the min PL value GT for spanning deletions
  m_produce_GT_with_min_PL_value_for_spanning_deletions =
    export_config->has_produce_gt_with_min_pl_value_for_spanning_deletions()
    ? export_config->produce_gt_with_min_pl_value_for_spanning_deletions() : false;
  //Disable file locking in TileDB
  m_disable_file_locking_in_tiledb = export_config->has_disable_file_locking_in_tiledb()
    ? export_config->disable_file_locking_in_tiledb() : false;
}

void GenomicsDBConfigBase::read_and_initialize_vid_and_callset_mapping_if_available(
    const genomicsdb_pb::ExportConfiguration* export_config) {
  //Callset mapping and vid mapping
  if(export_config->has_vid_mapping_file()) {
    m_vid_mapping_file = export_config->vid_mapping_file();
    m_vid_mapper = std::move(FileBasedVidMapper(m_vid_mapping_file));
  }
  else {
    if(export_config->has_vid_mapping()) {
      m_vid_mapper.parse_vidmap_protobuf(&(export_config->vid_mapping()));
    }
    else if (!m_vid_mapper.is_initialized()) {
      throw GenomicsDBConfigException("Protobuf ExportConfiguration must have \
          either vid_mapping_file or vid_mapping object defined");
    }
  }
  if(export_config->has_callset_mapping_file()) {
    m_callset_mapping_file = export_config->callset_mapping_file();
    auto callsets_json_doc = std::move(parse_json_file(m_callset_mapping_file));
    m_vid_mapper.read_callsets_info(callsets_json_doc, 0);
  }
  else {
    if(export_config->has_callset_mapping()) {
      m_vid_mapper.parse_callset_protobuf(&(export_config->callset_mapping()));
    }
    else if (!m_vid_mapper.is_initialized()) {
      throw GenomicsDBConfigException("Protobuf ExportConfiguration must have either \
	  callset_mapping_file or callset_mapping object defined");
    }
  }
}

void GenomicsDBConfigBase::get_pb_from_json_file(
        google::protobuf::Message *pb_config,
        const std::string& json_file) {
  char *json_buffer = 0;
  size_t json_buffer_length;
  if (TileDBUtils::read_entire_file(json_file, (void **)&json_buffer, &json_buffer_length) != TILEDB_OK
        || !json_buffer || json_buffer_length == 0) { 
    free(json_buffer);
    throw GenomicsDBConfigException(std::string("Could not open query JSON file ")+json_file);
  }
  std::string json_to_binary_output;
  google::protobuf::util::JsonParseOptions parse_opt;
  parse_opt.ignore_unknown_fields = true;
  //The function JsonStringToMessage was made available in Protobuf version 3.0.0. However,
  //to maintain compatibility with GATK-4, we need to use 3.0.0-beta-1. This version doesn't have
  //the JsonStringToMessage method. A workaround is as follows.
  //https://stackoverflow.com/questions/41651271/is-there-an-example-of-protobuf-with-text-output
  google::protobuf::util::TypeResolver* resolver= google::protobuf::util::NewTypeResolverForDescriptorPool(
      "", google::protobuf::DescriptorPool::generated_pool());
  auto status = google::protobuf::util::JsonToBinaryString(resolver,
      "/"+pb_config->GetDescriptor()->full_name(), json_buffer,
      &json_to_binary_output, parse_opt);
  if (!status.ok()) {
    delete resolver;
    free(json_buffer);
    throw GenomicsDBConfigException(std::string("Error converting JSON to binary string from file ")+json_file);
  }
  delete resolver;
  free(json_buffer);
  auto success = pb_config->ParseFromString(json_to_binary_output);
  if(!success)
    throw GenomicsDBConfigException(std::string("Could not parse query JSON file to protobuf ")+json_file);
}
