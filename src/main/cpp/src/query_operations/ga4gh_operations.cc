/**
 * @file ga4gh_operations.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
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
 *
 * @section DESCRIPTION
 *
 * Fill up GA4GH Variant and Call structures from GenomicsDB 
 *
 **/


#include "ga4gh_operations.h"
#include <google/protobuf/util/json_util.h>

GA4GHCallCreator::GA4GHCallCreator(const VidMapper& vid_mapper, std::ostream& stream, const size_t buffer_size)
  : m_variant(), m_call() {
  common_initialization(vid_mapper, stream, buffer_size, false); 
}

GA4GHCallCreator::GA4GHCallCreator(const VidMapper& vid_mapper, std::ostream& stream, const bool unlimited)
  : m_variant(), m_call() {
   common_initialization(vid_mapper, stream, 4096, unlimited); 
}

void GA4GHCallCreator::common_initialization(const VidMapper& vid_mapper, std::ostream& stream,
    const size_t buffer_size, const bool unlimited) {
  m_vid_mapper = &vid_mapper;
  m_stream = &stream;
  if(unlimited)
    m_limited_stream_buffer = false;
  else {
    m_limited_stream_buffer = true;
    m_stream_buffer_size = buffer_size;
  }
  m_contig_begin_tiledb_column_offset = -1;
  m_contig_end_tiledb_column_offset = -1;
  auto GT_field_info = m_vid_mapper->get_field_info("GT");
  if(GT_field_info == 0)
    throw VariantOperationException("No info found for GT field in vid mapping");
  m_ploidy_step_value = GT_field_info->m_length_descriptor.get_ploidy_step_value();
  m_first_call = true;
}

void GA4GHCallCreator::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell,
    const VariantQueryConfig& query_config, const VariantArraySchema& schema) {
  if(m_first_call) {
    (*m_stream) << "{\n\"variants\": [\n";
    m_first_call = false;
  }
  else
    (*m_stream) << ",";
  auto coords_ptr = cell.get_coordinates();
  //This if statement is NOT the common case (hopefully)
  if(coords_ptr[1] < m_contig_begin_tiledb_column_offset
      || coords_ptr[1] > m_contig_end_tiledb_column_offset) { //different contig from previous cell
    ContigInfo contig_info;
    auto status = m_vid_mapper->get_contig_info(coords_ptr[1], contig_info);
    if(!status)
      throw VariantOperationException(std::string("Could not find a contig corresponding to column ")
	  + std::to_string(coords_ptr[1])+" in vid mapping");
    m_contig_begin_tiledb_column_offset = contig_info.m_tiledb_column_offset;
    m_contig_end_tiledb_column_offset = contig_info.m_tiledb_column_offset + contig_info.m_length;
    m_contig = contig_info.m_name;
  }
  m_variant.Clear();
  //contig, start and end
  //GA4GH defines contig position as 0 based (unlike VCFs)
  m_variant.set_reference_name(m_contig);
  m_variant.set_start(coords_ptr[1]-m_contig_begin_tiledb_column_offset);
  auto end_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX);
  auto end_position = *(cell.get_field_ptr_for_query_idx<int64_t>(end_query_idx));
  m_variant.set_end(end_position-m_contig_begin_tiledb_column_offset);
  //RSID etc (ID field of VCFs)
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_ID_IDX)) {
    auto ID_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_ID_IDX);
    auto length = cell.get_field_length(ID_query_idx);
    if(length > 0u) {
      auto ptr = cell.get_field_ptr_for_query_idx<char>(ID_query_idx);
      auto curr_ID_ptr = ptr;
      for(auto i=0u;i<length;++i) {
	if(ptr[i] == ';') {
	  m_variant.add_names(curr_ID_ptr, ptr+i-curr_ID_ptr);
	  curr_ID_ptr = ptr+i+1u; //go past the ';'
	}
      }
      if(ptr+length-curr_ID_ptr > 0u)
	m_variant.add_names(curr_ID_ptr, ptr+length-curr_ID_ptr);
    }
  }
  //REF
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_REF_IDX)) {
    auto REF_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_REF_IDX);
    auto length = cell.get_field_length(REF_query_idx);
    if(length > 0u)
      m_variant.set_reference_bases(cell.get_field_ptr_for_query_idx<char>(REF_query_idx),
	  length);
  }
  //ALT
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_ALT_IDX)) {
    auto ALT_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_ALT_IDX);
    auto length = cell.get_field_length(ALT_query_idx);
    if(length > 0u) {
      auto ptr = cell.get_field_ptr_for_query_idx<char>(ALT_query_idx);
      auto curr_ALT_ptr = ptr;
      for(auto i=0u;i<length;++i) {
	if(ptr[i] == TILEDB_ALT_ALLELE_SEPARATOR[0]) {
	  auto allele_length = ptr + i - curr_ALT_ptr;
	  if(allele_length == 1u && IS_NON_REF_ALLELE(curr_ALT_ptr[0]))
	    m_variant.add_alternate_bases("<NON_REF>");
	  else
	    m_variant.add_alternate_bases(curr_ALT_ptr, ptr+i-curr_ALT_ptr);
	  curr_ALT_ptr = ptr+i+1u; //go past the ';'
	}
      }
      auto allele_length = ptr+length-curr_ALT_ptr;
      if(allele_length > 0u) {
	if(allele_length == 1u && IS_NON_REF_ALLELE(curr_ALT_ptr[0]))
	  m_variant.add_alternate_bases("<NON_REF>");
	else
	  m_variant.add_alternate_bases(curr_ALT_ptr, ptr+length-curr_ALT_ptr);
      }
    }
  }
  //FILTER
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_FILTER_IDX)) {
    auto FILTER_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_FILTER_IDX);
    auto length = cell.get_field_length(FILTER_query_idx);
    if(length > 0u) {
      auto ptr = cell.get_field_ptr_for_query_idx<int>(FILTER_query_idx);
      auto pass_filter_exists = false;
      for(auto i=0u;i<length;++i) {
	auto& field_info = m_vid_mapper->get_field_info(ptr[i]);
	if(field_info.m_name == "PASS")
	  pass_filter_exists = true;
	else
	  m_variant.add_filters_failed(field_info.m_name);
      }
      m_variant.set_filters_applied(true);
      m_variant.set_filters_passed(length == 1u && pass_filter_exists);
    }
    else {
      m_variant.set_filters_applied(false);
      m_variant.set_filters_passed(false);
    }
  }
  fill_ga4gh_call(query_config, cell);
  *(m_variant.add_calls()) = m_call;
  assert(m_stream);
  std::string s;
  google::protobuf::util::MessageToJsonString(m_variant, &s);
  (*m_stream) << s << "\n";
}

void GA4GHCallCreator::fill_ga4gh_call(const VariantQueryConfig& query_config,
    const GenomicsDBColumnarCell& cell) {
  m_call.Clear();
  auto coords_ptr = cell.get_coordinates();
  std::string callset_name;
  auto status = m_vid_mapper->get_callset_name(coords_ptr[0], callset_name);
  if(!status)
    throw VariantOperationException(std::string("No callset name found for row idx ")
	+ std::to_string(coords_ptr[0]));
  m_call.set_call_set_name(callset_name);
  m_call.set_call_set_id(callset_name);
  //Genotypes
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX)) {
    auto GT_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_GT_IDX);
    auto ptr = cell.get_field_ptr_for_query_idx<int>(GT_query_idx);
    auto length = cell.get_field_length(GT_query_idx);
    auto pb_GT_list_ptr = m_call.mutable_genotype();
    for(auto i=0u;i<length;i+=m_ploidy_step_value) {
      google::protobuf::Value x;
      //TODO: maybe reduce heap allocations
      x.set_string_value((ptr[i] == get_bcf_gt_no_call_allele_index<int>()) ? std::string(".")
	  : std::to_string(ptr[i]));
      *(pb_GT_list_ptr->add_values()) = std::move(x); //ListValue is weird
    }
  }
  //Genotype likelihoods
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_GL_IDX)) {
    auto GL_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_GL_IDX);
    auto ptr = cell.get_field_ptr_for_query_idx<float>(GL_query_idx);
    auto length = cell.get_field_length(GL_query_idx);
    for(auto i=0u;i<length;++i)
      m_call.add_genotype_likelihood(ptr[i]);
  }
  //Other attributes
  //Skip END, REF, ALT
  for (auto i=query_config.get_first_normal_field_query_idx();
      i<query_config.get_num_queried_attributes(); ++i) {
    auto known_field_enum = query_config.get_known_field_enum_for_query_idx(i);
    if(cell.is_valid(i) //valid
	&&(!query_config.is_defined_known_field_enum_for_query_idx(i) //Not known field
	|| (known_field_enum != GVCF_ID_IDX    //these fields are already handled in the GA4GH schema
	  && known_field_enum != GVCF_FILTER_IDX
	  && known_field_enum != GVCF_GT_IDX
	  && known_field_enum != GVCF_GL_IDX
	  ))) {
      cell.fill_ga4gh_attribute(&m_attr_list, i);
      m_call.mutable_attributes()->mutable_attr()->insert(google::protobuf::MapPair
	  <std::string, ga4gh::AttributeValueList>(query_config.get_query_attribute_name(i),
	  m_attr_list));
    }
  }
}

void GA4GHCallCreator::finalize() {
  (*m_stream) << "]\n}\n";
}
