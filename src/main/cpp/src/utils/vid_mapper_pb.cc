/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#include "vid_mapper_pb.h"
#include "known_field_info.h"
#include <ctype.h>
#include <cstdlib>
#include <cstdio>
#include "genomicsdb_logger.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw ProtoBufBasedVidMapperException(#X);

GenomicsDBProtoBufInitAndCleanup g_genomicsdb_protobuf_init_and_cleanup;

int VidMapper::parse_callset_protobuf(
  const CallsetMappingPB* callset_map_protobuf) {

  assert (callset_map_protobuf->IsInitialized() &&
          callset_map_protobuf->callsets_size()!= 0);

  auto num_callsets = callset_map_protobuf->callsets_size();
  m_row_idx_to_info.resize(num_callsets);
  m_max_callset_row_idx = -1;

  google::protobuf::RepeatedPtrField<SampleIDToTileDBIDMap>::const_iterator it =
    callset_map_protobuf->callsets().cbegin();

  for (; it != callset_map_protobuf->callsets().cend(); ++it) {
    SampleIDToTileDBIDMap sample_info = *it;
    std::string callset_name = it->sample_name();
    int64_t row_idx = sample_info.row_idx();
    int64_t idx_in_file = sample_info.idx_in_file();
    std::string stream_name = sample_info.has_stream_name() ? sample_info.stream_name() : "";

    if (m_callset_name_to_row_idx.find(callset_name) !=
        m_callset_name_to_row_idx.end()) {

      if (m_callset_name_to_row_idx[callset_name] != row_idx) {
        throw ProtoBufBasedVidMapperException(
          std::string("ERROR: Callset/sample ")
          + callset_name
          + " has two entries in the vid mapper Protobuf object with two TileDB row indexes: "
          + std::to_string(m_callset_name_to_row_idx[callset_name])
          + ", "
          + std::to_string(row_idx));
      }
    }

    m_max_callset_row_idx = std::max(m_max_callset_row_idx, row_idx);

    auto file_idx = get_or_append_global_file_idx(stream_name);
    if (static_cast<size_t>(row_idx) >= m_row_idx_to_info.size())
      m_row_idx_to_info.resize(static_cast<size_t>(row_idx)+1u);
    const auto& curr_row_info = m_row_idx_to_info[row_idx];

    // If row information for the same callset is found,
    // check whether it conflicts with current information.
    // If so, throw error
    if (curr_row_info.m_is_initialized) {
      if (curr_row_info.m_file_idx >= 0 && file_idx >= 0
          && curr_row_info.m_file_idx != file_idx)
        throw ProtoBufBasedVidMapperException(
          std::string("Callset/sample ")
          + callset_name
          + " has multiple stream names: "
          + m_file_idx_to_info[curr_row_info.m_file_idx].m_name
          + ", "
          + stream_name);
      if (curr_row_info.m_idx_in_file != idx_in_file)
        throw ProtoBufBasedVidMapperException(
          std::string("Conflicting values of \"idx_in_file\" ")
          + std::string("specified for Callset/sample ")
          + callset_name
          + ": "
          + std::to_string(curr_row_info.m_idx_in_file)
          + ", "
          + std::to_string(idx_in_file));
    } else {
      if(file_idx >= 0) {
	assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
	int64_t other_row_idx = 0ll;
	auto added_successfully = m_file_idx_to_info[file_idx].add_local_tiledb_row_idx_pair(
	    idx_in_file,
	    row_idx,
	    other_row_idx);
	if (!added_successfully)
	  throw ProtoBufBasedVidMapperException(std::string("Attempting to import a sample from file/stream ")+stream_name
	      +" multiple times under aliases '"+m_row_idx_to_info[other_row_idx].m_name
	      +"' and '"+callset_name+"' with row indexes "+std::to_string(other_row_idx)
	      +" and "+std::to_string(row_idx)+" respectively");
      }
    }

    m_callset_name_to_row_idx[callset_name] = row_idx;
    VERIFY_OR_THROW(static_cast<size_t>(row_idx) < m_row_idx_to_info.size());
    if (m_row_idx_to_info[row_idx].m_is_initialized &&
        m_row_idx_to_info[row_idx].m_name != callset_name)
      throw ProtoBufBasedVidMapperException(
        std::string("Callset/sample ")
        + callset_name
        + " has the same TileDB row index as "
        + m_row_idx_to_info[row_idx].m_name);

    m_row_idx_to_info[row_idx].set_info(
      row_idx,
      callset_name,
      file_idx,
      idx_in_file);
  }  // end of for all callsets

  check_for_missing_row_indexes();

  m_is_callset_mapping_initialized = true;
  m_is_initialized = m_is_contig_and_field_info_initialized;

  return GENOMICSDB_VID_MAPPER_SUCCESS;
} // end of parse_callset_protobuf

int VidMapper::parse_vidmap_protobuf(
  const VidMappingPB* vid_map_protobuf) {

  assert (vid_map_protobuf->IsInitialized() &&
          vid_map_protobuf->contigs_size()!=0 &&
          vid_map_protobuf->fields_size()!=0);

  int ret = 0;
  ret = parse_contigs_from_vidmap(vid_map_protobuf);
  assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);
  ret = parse_infofields_from_vidmap(vid_map_protobuf);
  assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);

  m_is_contig_and_field_info_initialized = true;
  m_is_initialized = m_is_callset_mapping_initialized;

  return GENOMICSDB_VID_MAPPER_SUCCESS;
}

int VidMapper::parse_contigs_from_vidmap(
  const VidMappingPB* vid_map_protobuf) {

  auto num_contigs = vid_map_protobuf->contigs_size();
  m_contig_idx_to_info.resize(num_contigs);
  m_contig_begin_2_idx.resize(num_contigs);
  m_contig_end_2_idx.resize(num_contigs);
  auto duplicate_contigs_exist = false;
  std::string contig_name;

  for (auto contig_idx = 0L; contig_idx < num_contigs; ++contig_idx) {
    contig_name = vid_map_protobuf->contigs(contig_idx).name();

    if (m_contig_name_to_idx.find(contig_name) != m_contig_name_to_idx.end()) {
      logger.warn("Contig/chromosome name {} appears more than once in vid map", contig_name);
      duplicate_contigs_exist = true;
      continue;
    }

    auto tiledb_column_offset =
      vid_map_protobuf->contigs(contig_idx).tiledb_column_offset();

    VERIFY_OR_THROW(tiledb_column_offset >= 0LL);
    auto length = vid_map_protobuf->contigs(contig_idx).length();

    VERIFY_OR_THROW(length >= 0LL);

    m_contig_name_to_idx[contig_name] = contig_idx;
    m_contig_idx_to_info[contig_idx].set_info(
      contig_idx,
      contig_name,
      length,
      tiledb_column_offset);
    m_contig_begin_2_idx[contig_idx].first = tiledb_column_offset;
    m_contig_begin_2_idx[contig_idx].second = contig_idx;
    m_contig_end_2_idx[contig_idx].first =
      tiledb_column_offset + length - 1; //inclusive
    m_contig_end_2_idx[contig_idx].second = contig_idx;

    if (duplicate_contigs_exist) {
      throw ProtoBufBasedVidMapperException(
        std::string("Duplicate contigs found: ")
        + contig_name);
    }
  }
  std::sort(
    m_contig_begin_2_idx.begin(),
    m_contig_begin_2_idx.end(),
    contig_offset_idx_pair_cmp);
  std::sort(
    m_contig_end_2_idx.begin(),
    m_contig_end_2_idx.end(),
    contig_offset_idx_pair_cmp);

  // Check that there are no spurious overlaps.
  // If found, throw an exception
  auto last_contig_idx = -1;
  auto last_contig_end_column = -1ll;
  auto overlapping_contigs_exist = false;
  for (auto contig_idx = 0UL; contig_idx < m_contig_begin_2_idx.size();
       ++contig_idx)
//  for (const auto& offset_idx_pair : m_contig_begin_2_idx)
  {
//    auto contig_idx = offset_idx_pair.second;
    const auto& contig_info = m_contig_idx_to_info[contig_idx];
    if (last_contig_idx >= 0) {
      const auto& last_contig_info = m_contig_idx_to_info[last_contig_idx];
      if (contig_info.m_tiledb_column_offset <= last_contig_end_column) {
        logger.info("Contig/chromosome {} begins at TileDB column {} and intersects with contig/chromosome {} that spans columns [{}, {}]",
                     contig_info.m_name,
                     contig_info.m_tiledb_column_offset,
                     last_contig_info.m_name,
                     last_contig_info.m_tiledb_column_offset,
                     last_contig_info.m_tiledb_column_offset + last_contig_info.m_length-1);
        overlapping_contigs_exist = true;
      }
    }
    last_contig_idx = contig_idx;
    last_contig_end_column =
      contig_info.m_tiledb_column_offset + contig_info.m_length - 1;
  }

  if (overlapping_contigs_exist) {
    throw ProtoBufBasedVidMapperException(
      std::string("Overlapping contigs found"));
  }

  return GENOMICSDB_VID_MAPPER_SUCCESS;
} // end of parse_contigs_from_vidmap

int VidMapper::parse_infofields_from_vidmap(
  const VidMappingPB* vid_map_protobuf) {

  auto num_fields = vid_map_protobuf->fields_size();
  m_field_idx_to_info.resize(num_fields);
  auto duplicate_fields_exist = false;

  for (auto pb_field_idx = 0, field_idx = 0; pb_field_idx < num_fields;
       ++pb_field_idx) {
    const auto& field_name = vid_map_protobuf->fields(pb_field_idx).name();
    if (m_field_name_to_idx.find(field_name) != m_field_name_to_idx.end()) {
      logger.warn("Duplicate field name {} found in vid map", field_name);
      duplicate_fields_exist = true;
      continue;
    }

    // Known fields
    auto known_field_enum = 0u;
    auto is_known_field =
      KnownFieldInfo::get_known_field_enum_for_name(
        field_name,
        known_field_enum);
    // Map
    m_field_name_to_idx[field_name] = field_idx;
    m_field_idx_to_info[field_idx].set_info(field_name, field_idx);
    auto& ref = m_field_idx_to_info[field_idx];

    if(vid_map_protobuf->fields(pb_field_idx).has_vcf_name())
      ref.m_vcf_name = vid_map_protobuf->fields(pb_field_idx).vcf_name();

    // VCF class type can be an array of values: INFO, FORMAT and FILTER
    auto class_type_size =
      vid_map_protobuf->fields(pb_field_idx).vcf_field_class_size();

    if (class_type_size > 0L) {
      for (int i = 0; i < class_type_size; ++i) {
        std::string class_name =
          vid_map_protobuf->fields(pb_field_idx).vcf_field_class(i);
        if (class_name == "INFO")
          ref.m_is_vcf_INFO_field = true;
        else if (class_name == "FORMAT") {
          ref.m_is_vcf_FORMAT_field = true;
        } else if (class_name == "FILTER")
          ref.m_is_vcf_FILTER_field = true;
      }
    }

    //Length descriptor
    if (vid_map_protobuf->fields(pb_field_idx).length_size() > 0) {
      ref.m_length_descriptor.resize(vid_map_protobuf->fields(pb_field_idx).length_size());
      for (auto i=0; i<vid_map_protobuf->fields(pb_field_idx).length_size(); ++i) {
        auto& pb_length_descriptor_component = vid_map_protobuf->fields(pb_field_idx).length(i);
        if (pb_length_descriptor_component.has_fixed_length())
          ref.m_length_descriptor.set_num_elements(i, pb_length_descriptor_component.fixed_length());
        else {
          assert(pb_length_descriptor_component.has_variable_length_descriptor());
          parse_string_length_descriptor(field_name.c_str(),
                                         pb_length_descriptor_component.variable_length_descriptor().c_str(),
                                         pb_length_descriptor_component.variable_length_descriptor().size(),
                                         ref.m_length_descriptor, i);
        }
      }
    } else {
      if (is_known_field) {
        auto length_descriptor_code =
          KnownFieldInfo::get_length_descriptor_for_known_field_enum(
            known_field_enum);
        ref.m_length_descriptor.set_length_descriptor(0u, length_descriptor_code);
        if (length_descriptor_code == BCF_VL_FIXED)
          ref.m_length_descriptor.set_num_elements(0u,
              KnownFieldInfo::get_num_elements_for_known_field_enum(
                known_field_enum,
                0u,
                0u)
                                                  );  //don't care about ploidy
      }
    }

    //Field type - int, char, or <int, float> (allele specific annotations in GATK)
    if (vid_map_protobuf->fields(pb_field_idx).type_size() == 0u)
      throw VidMapperException(std::string("Attribute 'type' is mandatory for all fields in GenomicsDB ")
                               +" field "+field_name+" missing 'type' in Protobuf structure");
    FieldElementTypeDescriptor type_descriptor(vid_map_protobuf->fields(pb_field_idx).type_size());
    for (auto i=0u; i<static_cast<unsigned>(vid_map_protobuf->fields(pb_field_idx).type_size()); ++i) {
      const auto& field_type = vid_map_protobuf->fields(pb_field_idx).type(i);
      auto type_index_ht_type_pair = get_type_index_and_bcf_ht_type(field_type.c_str());
      type_descriptor.set_tuple_element_type(i, type_index_ht_type_pair.first, type_index_ht_type_pair.second);
    }
    ref.set_type(type_descriptor);

    //Sometimes the VCF type can be different from the real datatype of the field
    //For example, for multi-D vectors, the VCF type is string: @$@#@#!#$%$%
    if (vid_map_protobuf->fields(pb_field_idx).has_vcf_type()) {
      auto type_index_ht_type_pair = get_type_index_and_bcf_ht_type(vid_map_protobuf->fields(pb_field_idx).vcf_type().c_str());
      type_descriptor.resize_num_elements_in_tuple(1u);
      type_descriptor.set_tuple_element_type(0u, type_index_ht_type_pair.first, type_index_ht_type_pair.second);
      ref.set_vcf_type(type_descriptor);
    }
    ref.modify_field_type_if_multi_dim_field();

    if (vid_map_protobuf->fields(pb_field_idx).has_vcf_field_combine_operation())
      set_VCF_field_combine_operation(ref,
                                      vid_map_protobuf->fields(pb_field_idx).vcf_field_combine_operation().c_str());
    else
      if(is_known_field)
	ref.m_VCF_field_combine_operation =
	  KnownFieldInfo::get_VCF_field_combine_operation_for_known_field_enum(known_field_enum);

    //Generally used when multi-D vectors are represented as delimited strings in the VCF
    //Used mostly in conjunction with the vcf_type attribute
    if (vid_map_protobuf->fields(pb_field_idx).vcf_delimiter_size() > 0) {
      for (auto i=0; i<vid_map_protobuf->fields(pb_field_idx).vcf_delimiter_size(); ++i)
        ref.m_length_descriptor.set_vcf_delimiter(i,
            vid_map_protobuf->fields(pb_field_idx).vcf_delimiter(i).c_str());
    }
    ref.m_disable_remap_missing_with_non_ref = 
      vid_map_protobuf->fields(pb_field_idx).has_disable_remap_missing_with_non_ref() ?
      vid_map_protobuf->fields(pb_field_idx).disable_remap_missing_with_non_ref() : false;

    ++field_idx;
    flatten_field(field_idx, field_idx-1);
  } // for (auto field_idx = 0; field_idx < num_fields; ++field_idx)

  add_mandatory_fields();

  if (duplicate_fields_exist) {
    throw ProtoBufBasedVidMapperException(
      std::string("Duplicate fields exist in vid map"));
  }

  return GENOMICSDB_VID_MAPPER_SUCCESS;
} // end of parse_infofields_from_vidmap

ProtoBufBasedVidMapperException::~ProtoBufBasedVidMapperException() { ; }

ProtoBufBasedVidMapperException::ProtoBufBasedVidMapperException(
  const std::string m) : msg_("ProtoBufBasedVidMapperException : "+m) { ; }
