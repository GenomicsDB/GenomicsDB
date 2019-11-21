/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#ifdef HTSDIR

#include "broad_combined_gvcf.h"
#include "genomicsdb_multid_vector_field.h"

#ifdef DO_MEMORY_PROFILING
#include "memory_measure.h"
#define ONE_MB 1024ull*1024ull
#endif

//INFO fields
#define MAKE_BCF_INFO_TUPLE(enum_idx, query_idx, field_info) \
  INFO_tuple_type(enum_idx, query_idx, field_info)
#define BCF_INFO_GET_KNOWN_FIELD_ENUM(X) (std::get<0>(X))
#define BCF_INFO_GET_QUERY_FIELD_IDX(X) (std::get<1>(X))
#define BCF_INFO_GET_FIELD_INFO_PTR(X) (std::get<2>(X))
#define BCF_INFO_GET_BCF_HT_TYPE(X) (std::get<2>(X)->get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u))
#define BCF_INFO_GET_VCF_FIELD_NAME(X) (std::get<2>(X)->m_vcf_name)
#define BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(X) (std::get<2>(X)->m_VCF_field_combine_operation)
//FORMAT fields
#define MAKE_BCF_FORMAT_TUPLE(enum_idx, query_idx, field_info) \
  FORMAT_tuple_type(enum_idx, query_idx, field_info)
#define BCF_FORMAT_GET_KNOWN_FIELD_ENUM(X) (std::get<0>(X))
#define BCF_FORMAT_GET_QUERY_FIELD_IDX(X) (std::get<1>(X))
#define BCF_FORMAT_GET_BCF_HT_TYPE(X) (std::get<2>(X)->get_vcf_type().get_tuple_element_bcf_ht_type(0u))
#define BCF_FORMAT_GET_VCF_FIELD_NAME(X) (std::get<2>(X)->m_vcf_name)

#define GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(X) (std::get<2>(X))

//Static member
const std::unordered_set<char> BroadCombinedGVCFOperator::m_legal_bases({'A', 'T', 'G', 'C'});

//Utility functions for encoding GT field
template<bool phase_information_in_TileDB, bool produce_GT_field>
int encode_GT_element(const int value, const bool is_phased);

//Specialization - phase information exists in TileDB and produce_GT_field
template<>
inline int encode_GT_element<true, true>(const int value, const bool is_phased) {
  return is_bcf_valid_value(value)
         ? (is_phased ? bcf_gt_phased(value) : bcf_gt_unphased(value))
         : value;
}

//Specialization - phase information exists in TileDB and !produce_GT_field
template<>
inline int encode_GT_element<true, false>(const int value, const bool is_phased) {
  return is_bcf_valid_value(value)
         ? (is_phased ? (bcf_gt_missing | 1) : bcf_gt_missing)
         : value;
}
//Specialization - phase information doesn't exist in TileDB and produce_GT_field
template<>
inline int encode_GT_element<false, true>(const int value, const bool is_phased) {
  return is_bcf_valid_value(value) ? bcf_gt_unphased(value) : value;
}

//Specialization - phase information doesn't exist in TileDB and !produce_GT_field
template<>
inline int encode_GT_element<false, false>(const int value, const bool is_phased) {
  return is_bcf_valid_value(value) ? bcf_gt_missing : value;
}

template<bool phase_information_in_TileDB, bool produce_GT_field>
void encode_GT_vector(int* inout_vec, const uint64_t input_offset,
                      const unsigned num_elements_per_sample, uint64_t& output_idx);

//Specialization - phase information exists in TileDB and produce_GT_field
template<>
inline void encode_GT_vector<true, true>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx) {
  if (num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, true>(inout_vec[input_offset], false);
  //skip over phasing elements
  for (auto k=2u; k<num_elements_per_sample; k+=2u)
    inout_vec[output_idx++] = encode_GT_element<true, true>(inout_vec[input_offset+k], (inout_vec[input_offset+k-1u] > 0));
}

//Specialization - phase information exists in TileDB and !produce_GT_field
template<>
inline void encode_GT_vector<true, false>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx) {
  if (num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, false>(inout_vec[input_offset], false);
  //skip over phasing elements
  for (auto k=2u; k<num_elements_per_sample; k+=2u)
    inout_vec[output_idx++] = encode_GT_element<true, false>(inout_vec[input_offset+k], (inout_vec[input_offset+k-1u] > 0));
}

//Specialization - phase information doesn't exist in TileDB and produce_GT_field
template<>
inline void encode_GT_vector<false, true>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx) {
  if (num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, true>(inout_vec[input_offset], false);
  //skip over phasing elements
  for (auto k=1u; k<num_elements_per_sample; ++k)
    inout_vec[output_idx++] = encode_GT_element<false, true>(inout_vec[input_offset+k], false);
}

//Specialization - phase information doesn't exist in TileDB and !produce_GT_field
template<>
inline void encode_GT_vector<false, false>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx) {
  if (num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, false>(inout_vec[input_offset], false);
  //skip over phasing elements
  for (auto k=1u; k<num_elements_per_sample; ++k)
    inout_vec[output_idx++] = encode_GT_element<false, false>(inout_vec[input_offset+k], false);
}

BroadCombinedGVCFOperator::BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VidMapper& id_mapper,
    const VariantQueryConfig& query_config,
    const bool use_missing_values_only_not_vector_end)
  : GA4GHOperator(query_config, id_mapper, true) {
  clear();
  if (!id_mapper.is_initialized())
    throw BroadCombinedGVCFException("Id mapper is not initialized");
  if (!id_mapper.is_callset_mapping_initialized())
    throw BroadCombinedGVCFException("Callset mapping in id mapper is not initialized");
  m_query_config = &query_config;
  //Initialize VCF structs
  m_vcf_adapter = &vcf_adapter;
  m_vid_mapper = &id_mapper;
  m_use_missing_values_not_vector_end = use_missing_values_only_not_vector_end;
  m_vcf_hdr = vcf_adapter.get_vcf_header();
  m_bcf_out = bcf_init();
  //vector of char*, to avoid frequent reallocs()
  m_alleles_pointer_buffer.resize(100u);
  //DP INFO field - handle after all FORMAT fields have been processed
  FORMAT_tuple_type DP_INFO_as_FORMAT_tuple;
  auto is_DP_INFO_queried = false;
  //Determine queried INFO and FORMAT fields
  for (auto i=0u; i<query_config.get_num_queried_attributes(); ++i) {
    auto* field_info = m_vid_mapper->get_field_info(query_config.get_query_attribute_name(i));
    if (field_info) {
      auto known_field_enum = query_config.is_defined_known_field_enum_for_query_idx(i) ? query_config.get_known_field_enum_for_query_idx(i)
                              : UNDEFINED_ATTRIBUTE_IDX_VALUE;
      auto VCF_field_combine_operation = query_config.get_VCF_field_combine_operation_for_query_attribute_idx(i);
      auto sites_only_query = query_config.sites_only_query();
      auto add_to_INFO_vector = (field_info->m_is_vcf_INFO_field && known_field_enum != GVCF_END_IDX
                                 && (known_field_enum != GVCF_DP_IDX || VCF_field_combine_operation != VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_DP) //not DP or not combined as GATK Combine GVCF DP
                                 && VCF_field_combine_operation != VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT  //not moved to FORMAT
                                );
      auto add_to_FORMAT_vector = (
                                    (field_info->m_is_vcf_FORMAT_field && (!sites_only_query
                                        || known_field_enum == GVCF_DP_FORMAT_IDX
                                        || known_field_enum == GVCF_MIN_DP_IDX)) //FORMAT field when !sites_only_query or FORMAT field specifically GVCF_DP_FORMAT_IDX or GVCF_MIN_DP_IDX
                                    ||
                                    (field_info->m_is_vcf_INFO_field
                                     && ((known_field_enum == GVCF_DP_IDX && VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_DP)
                                         || (VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT && !sites_only_query)
                                        )
                                    )
                                  );
      if (add_to_INFO_vector) {
        switch (VCF_field_combine_operation) {
        case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_NONE:
          break;
        case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION:
          std::cerr << "WARNING: No valid combination operation found for INFO field "<<field_info->m_vcf_name<<" - the field will NOT be part of INFO fields in the generated VCF records\n";
          break;
        case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_HISTOGRAM_SUM: {
          m_INFO_fields_vec.emplace_back(MAKE_BCF_INFO_TUPLE(known_field_enum, i, field_info));
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info->m_vcf_name, BCF_HL_INFO);
          assert(field_info->is_flattened_field());
          auto parent_field_vid_idx = field_info->get_parent_composite_field_idx();
          auto& parent_field_vid_info = m_vid_mapper->get_field_info(parent_field_vid_idx);
          if (parent_field_vid_info.get_genomicsdb_type().get_num_elements_in_tuple() != 2u)
            throw BroadCombinedGVCFException(std::string(
                                               "Operation histogram_sum is only supported for fields whose elements are tuple with 2 constituent elements; field ")
                                             +parent_field_vid_info.m_name+" does not satisfy this requirement");
          if (field_info->get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u) != BCF_HT_INT
              && field_info->get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u) != BCF_HT_REAL)
            throw BroadCombinedGVCFException(std::string(
                                               "Operation histogram_sum is only supported for fields whose elements are tuple with 2 constituent elements")
                                             +" that are int or float; field "
                                             +parent_field_vid_info.m_name+" does not satisfy this requirement");
          auto iter_flag_pair = m_INFO_histogram_field_map.insert(std::pair<unsigned, INFO_histogram_field_tuple_type>(
                                  parent_field_vid_idx, INFO_histogram_field_tuple_type(0u, 0u, 0)));
          auto& curr_INFO_histogram_tuple = (*(iter_flag_pair.first)).second;
          if (field_info->get_element_index_in_tuple() == 0u)
            std::get<0>(curr_INFO_histogram_tuple) = i;
          else
            std::get<1>(curr_INFO_histogram_tuple) = i;
	  if(GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(curr_INFO_histogram_tuple) == 0) { //null, must allocate
	    HistogramFieldHandlerBase* ptr = 0;
	    auto is_bin_real =
	      (parent_field_vid_info.get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u) == BCF_HT_REAL)
	      ? 1 : 0;
	    auto is_count_real =
	      (parent_field_vid_info.get_genomicsdb_type().get_tuple_element_bcf_ht_type(1u) == BCF_HT_REAL)
	      ? 1 : 0;
	    switch((is_bin_real << 1) | is_count_real) {
	      case 0: //both int
		ptr = new HistogramFieldHandler<int, int, int64_t>();
		break;
	      case 1: //int, float
		ptr = new HistogramFieldHandler<int, float>();
		break;
	      case 2: //float, int
		ptr = new HistogramFieldHandler<float, int, int64_t>();
		break;
	      case 3: //both float
		ptr = new HistogramFieldHandler<float, float>();
		break;
	    }
	    assert(ptr);
	    GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(curr_INFO_histogram_tuple) = ptr;
	  }
          break;
        }
        default:
          m_INFO_fields_vec.emplace_back(MAKE_BCF_INFO_TUPLE(known_field_enum, i, field_info));
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info->m_vcf_name, BCF_HL_INFO);
          break;
        }
      }
      if (add_to_FORMAT_vector) {
        auto format_tuple = MAKE_BCF_FORMAT_TUPLE(known_field_enum, i, field_info);
        if ((field_info->m_is_vcf_FORMAT_field)
            || (VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT)
           ) {
          m_FORMAT_fields_vec.emplace_back(format_tuple);
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info->m_vcf_name, BCF_HL_FMT);
        } else { //DP INFO
          DP_INFO_as_FORMAT_tuple = std::move(format_tuple);
          is_DP_INFO_queried = true;
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, "DP", BCF_HL_INFO);
        }
      }
    }
  }
  //is FILTER field queried?
  if (query_config.is_defined_query_idx_for_known_field_enum(GVCF_FILTER_IDX)) {
    for (auto i=0u; i<m_vid_mapper->get_num_fields(); ++i) {
      auto& field_info = m_vid_mapper->get_field_info(i);
      if (field_info.m_is_vcf_FILTER_field)
        VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, m_vid_mapper, field_info.m_vcf_name, BCF_HL_FLT);
    }
  }
  //If DP is queried, it should always be the last field in m_FORMAT_fields_vec
  if (is_DP_INFO_queried) {
    //Move to the last element
    m_FORMAT_fields_vec.emplace_back(DP_INFO_as_FORMAT_tuple);
    if (BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u]) != GVCF_DP_IDX)
      throw BroadCombinedGVCFException("Last queried FORMAT field should be DP, instead it is "
                                       +g_known_variant_field_names[BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u])]);
  }
  //Qual combine operation
  auto* field_info = m_vid_mapper->get_field_info("QUAL");
  m_vcf_qual_tuple = MAKE_BCF_INFO_TUPLE(KnownVariantFieldsEnum::GVCF_QUAL_IDX,
                                         query_config.get_query_idx_for_known_field_enum(KnownVariantFieldsEnum::GVCF_QUAL_IDX),
                                         field_info);
  //Add missing contig names to template header
  for (auto i=0u; i<m_vid_mapper->get_num_contigs(); ++i) {
    auto& contig_info = m_vid_mapper->get_contig_info(i);
    auto& contig_name = contig_info.m_name;
    if (bcf_hdr_name2id(m_vcf_hdr, contig_name.c_str()) < 0) {
      std::string contig_vcf_line = std::string("##contig=<ID=")+contig_name+",length="
                                    +std::to_string(contig_info.m_length)+">";
      int line_length = 0;
      auto hrec = bcf_hdr_parse_line(m_vcf_hdr, contig_vcf_line.c_str(), &line_length);
      bcf_hdr_add_hrec(m_vcf_hdr, hrec);
      bcf_hdr_sync(m_vcf_hdr);
    }
  }
  //Get contig info for position 0, store curr contig in next_contig and call switch_contig function to do all the setup
  auto curr_contig_flag = m_vid_mapper->get_next_contig_location(-1ll, m_next_contig_name, m_next_contig_begin_position);
  assert(curr_contig_flag);
  switch_contig();
  //Add samples to template header if this is not a sites only query
  if (!(m_query_config->sites_only_query())) {
    std::string callset_name;
    for (auto i=0ull; i<query_config.get_num_rows_to_query(); ++i) {
      auto row_idx = query_config.get_array_row_idx_for_query_row_idx(i);
      auto status = m_vid_mapper->get_callset_name(row_idx, callset_name);
      if (!status || callset_name.empty())
        throw BroadCombinedGVCFException(std::string("No sample/CallSet name specified in JSON file/Protobuf object for TileDB row ")
                                         + std::to_string(row_idx));
      auto add_sample_status = bcf_hdr_add_sample(m_vcf_hdr, callset_name.c_str());
      if (add_sample_status < 0)
        throw BroadCombinedGVCFException(std::string("Could not add sample ")
                                         +callset_name+" to the combined VCF/gVCF header");
    }
  }
  bcf_hdr_sync(m_vcf_hdr);
  //Map from vid mapper field idx to hdr field idx
  m_global_field_idx_to_hdr_idx.resize(m_vid_mapper->get_num_fields(), -1);
  for (auto i=0u; i<m_vid_mapper->get_num_fields(); ++i) {
    auto& field_info = m_vid_mapper->get_field_info(i);
    //Could be -1
    m_global_field_idx_to_hdr_idx[i] = bcf_hdr_id2int(m_vcf_hdr, BCF_DT_ID, field_info.m_vcf_name.c_str());
  }
  //Since the header might have been modified and might be streamed into Java and BCF2Codec doesn't
  //deal with BCF IDX attributes, clean up the header completely by serializing and deserializing
  std::vector<uint8_t> tmp_buffer(10000u);
  while (bcf_hdr_serialize(m_vcf_hdr, &(tmp_buffer[0u]), 0u, tmp_buffer.size()-1u, 0, 0) == 0u)
    tmp_buffer.resize(2*tmp_buffer.size());
  bcf_hdr_destroy(m_vcf_hdr);
  m_vcf_hdr = bcf_hdr_init("r");
  auto hdr_offset = bcf_hdr_deserialize(m_vcf_hdr, &(tmp_buffer[0u]), 0u, tmp_buffer.size(), 0);
  assert(hdr_offset > 0u);
  m_vcf_adapter->set_vcf_header(m_vcf_hdr);
  m_vcf_adapter->print_header();
  //vector of field pointers used for handling remapped fields when dealing with spanning deletions
  //Individual pointers will be allocated later
  m_spanning_deletions_remapped_fields.resize(m_remapped_fields_query_idxs.size());
  //Initialize GT encoding function pointer
  if (query_config.is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX)) {
    const auto& GT_length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(m_GT_query_idx);
    if (GT_length_descriptor.contains_phase_information()) {
      //GATK CombineGVCF does not produce GT field by default - option to produce GT
      if (query_config.produce_GT_field())
        m_encode_GT_vector_function_ptr = encode_GT_vector<true, true>;
      else
        m_encode_GT_vector_function_ptr = encode_GT_vector<true, false>;
    } else {
      //GATK CombineGVCF does not produce GT field by default - option to produce GT
      if (query_config.produce_GT_field())
        m_encode_GT_vector_function_ptr = encode_GT_vector<false, true>;
      else
        m_encode_GT_vector_function_ptr = encode_GT_vector<false, false>;
    }
  } else
    m_encode_GT_vector_function_ptr = encode_GT_vector<true, false>; //encode phase, but don't produce GT
#ifdef DO_MEMORY_PROFILING
  m_next_memory_limit = 100*ONE_MB;
#endif
  m_bcf_record_size = 0ull;
  m_should_add_GQ_field = true; //always added in new version of CombineGVCFs
}

BroadCombinedGVCFOperator::~BroadCombinedGVCFOperator() {
  for(auto& key_tuple_pair : m_INFO_histogram_field_map) {
    auto& info_tuple = key_tuple_pair.second;
    if(GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(info_tuple) != 0) {
      delete GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(info_tuple);
      GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(info_tuple) = 0;
    }
  }
  bcf_destroy(m_bcf_out);
  clear();
#ifdef DO_PROFILING
  m_bcf_t_creation_timer.print("bcf_t creation time", std::cerr);
#endif
}

void BroadCombinedGVCFOperator::clear() {
  m_curr_contig_name.clear();
  m_next_contig_name.clear();
  m_alleles_pointer_buffer.clear();
  m_INFO_fields_vec.clear();
  m_FORMAT_fields_vec.clear();
  m_MIN_DP_vector.clear();
  m_DP_FORMAT_vector.clear();
  m_spanning_deletions_remapped_fields.clear();
  m_spanning_deletion_remapped_GT.clear();
  m_spanning_deletion_current_genotype.clear();
  m_global_field_idx_to_hdr_idx.clear();
  m_FILTER_idx_vec.clear();
}

bool BroadCombinedGVCFOperator::remap_if_needed_and_combine(const Variant& variant,
    const unsigned query_field_idxs[],
    const VCFFieldCombineOperationEnum combine_op,
    const void** output_ptr,
    unsigned& num_result_elements) {
  const auto is_histogram = (combine_op == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_HISTOGRAM_SUM);
  auto num_iter = is_histogram ? 2u: 1u;
  const FieldLengthDescriptor* length_descriptors[2];
  const FieldInfo* field_info_ptrs[2];
  for(auto i=0u;i<num_iter;++i) {
   length_descriptors[i] =  &(m_query_config->get_length_descriptor_for_query_attribute_idx(query_field_idxs[i]));
   field_info_ptrs[i] = m_query_config->get_field_info_for_query_attribute_idx(query_field_idxs[i]);
  }
  const auto remapping_needed =  (m_remapping_needed && length_descriptors[0]->is_length_allele_dependent());
  //Fields such as PL should be skipped, if the #alleles is above a threshold
  if (remapping_needed
      && length_descriptors[0u]->is_length_allele_dependent()
      && check_if_too_many_alleles_and_print_message(variant, *(length_descriptors[0u])))
    return false;
  //Remapper for m_remapped_variant
  RemappedVariant remapper_variants[2] = { RemappedVariant(m_remapped_variant, query_field_idxs[0u], true),
    RemappedVariant(m_remapped_variant, query_field_idxs[1u], true) };
  //Get HistogramFieldHandlerBase ptr for histogram combine operation
  HistogramFieldHandlerBase* handler_ptr = 0;
  if(is_histogram) {
    const auto parent_field_vid_idx = field_info_ptrs[0]->get_parent_composite_field_idx();
    assert(m_INFO_histogram_field_map.find(parent_field_vid_idx) != m_INFO_histogram_field_map.end());
    handler_ptr = GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(m_INFO_histogram_field_map[parent_field_vid_idx]);
  }
  bool is_first_iter[] = { true, true };
  bool valid_result_found[] = { false, false };
  //Iterate over valid calls - m_remapped_variant and variant have same list of valid calls
  for (auto iter=m_remapped_variant.begin(); iter!=m_remapped_variant.end(); ++iter) {
    auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
    const auto& original_call = variant.get_call(curr_call_idx_in_variant);
    //Only use the buffer from the first call to store remapped data - saves memory
    auto& remapped_call = m_remapped_variant.get_call(0u);
    for (auto i=0u;i<num_iter;++i) {
      auto query_field_idx = query_field_idxs[i];
      const auto& length_descriptor = *(length_descriptors[i]);;
      const auto& original_field = original_call.get_field(query_field_idx);
      assert(field_info_ptrs[i]->get_genomicsdb_type().get_num_elements_in_tuple() == 1u);
      const auto bcf_ht_type = field_info_ptrs[i]->get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u);
      auto& remapped_field = remapped_call.get_field(query_field_idx);
      if(remapping_needed)
	remap_if_needed(variant,
	    *(m_query_config),
	    curr_call_idx_in_variant,
	    query_field_idx,
	    remapped_field,
	    remapper_variants[i],
	    length_descriptor);
      switch(combine_op) {
	case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_SUM:
	  valid_result_found[i] = m_field_handlers[bcf_ht_type]->get_valid_sum(
	      remapping_needed ? remapped_field : original_field, is_first_iter[i])
	    || valid_result_found[i];
	  break;
	case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM:
	  if (length_descriptor.get_num_dimensions() == 1u)
	    valid_result_found[i] = m_field_handlers[bcf_ht_type]->compute_valid_element_wise_sum(
		remapping_needed ? remapped_field : original_field, is_first_iter[i])
	      || valid_result_found[i];
	  else
	    valid_result_found[i] = m_field_handlers[bcf_ht_type]->compute_valid_element_wise_sum_2D_vector(
		remapping_needed ? remapped_field : original_field, *(field_info_ptrs[i]), is_first_iter[i])
	      || valid_result_found[i];
	  break;
	case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_HISTOGRAM_SUM:
	  //Update histogram only after both components of the histogram - bin and count are remapped
	  break;
	default:
	  throw BroadCombinedGVCFException(std::string("Unknown combine operation for incremental computation ")
		+ std::to_string(combine_op));
      }
    }
    //Update histogram only after both components of the histogram - bin and count are remapped
    if(is_histogram) {
      assert(handler_ptr);
      valid_result_found[0] = handler_ptr->compute_valid_histogram_sum_2D_vector(
	  remapping_needed ? remapped_call.get_field(query_field_idxs[0])
	  : original_call.get_field(query_field_idxs[0]),
	  remapping_needed ? remapped_call.get_field(query_field_idxs[1])
	  : original_call.get_field(query_field_idxs[1]),
	  field_info_ptrs[0],
	  field_info_ptrs[1],
	  is_first_iter[0]
	  ) || valid_result_found[0];
    }
    for(auto i=0u;i<num_iter;++i)
      is_first_iter[i] = false;
  }
  const auto bcf_ht_type = field_info_ptrs[0u]->get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u);
  //Get pointers to results
  switch(combine_op) {
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_SUM:
      *output_ptr = m_field_handlers[bcf_ht_type]->get_pointer_to_sum();
      num_result_elements = 1u;
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM:
      if (length_descriptors[0u]->get_num_dimensions() == 1u)
	m_field_handlers[bcf_ht_type]->get_pointer_to_element_wise_sum(output_ptr, num_result_elements);
      else {
	//caller deals with stringification
	num_result_elements = 0;
      }
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_HISTOGRAM_SUM:
      //caller deals with stringification
      num_result_elements = 0;
      break;
    default:
      throw BroadCombinedGVCFException(std::string("Unknown combine operation for incremental computation ")
	  + std::to_string(combine_op));
  }
  return (valid_result_found[0] || valid_result_found[1]);
}

bool BroadCombinedGVCFOperator::handle_VCF_field_combine_operation(const Variant& variant,
    const INFO_tuple_type& curr_tuple, void*& result_ptr, unsigned& num_result_elements) {
  auto valid_result_found = false;
  auto query_field_idx = BCF_INFO_GET_QUERY_FIELD_IDX(curr_tuple);
  auto length_descriptor = m_query_config->get_length_descriptor_for_query_attribute_idx(query_field_idx);
  //Fields such as PL are skipped if the #alleles is above a certain threshold
  if (length_descriptor.is_length_genotype_dependent()
      && too_many_alt_alleles_for_genotype_length_fields(m_merged_alt_alleles.size()))
    return false;
  //Check if this is a field that was remapped - for remapped fields, we must use field objects from m_remapped_variant
  //else we should use field objects from the original variant
  auto& src_variant = (m_remapping_needed && (length_descriptor.is_length_allele_dependent()
                       || query_field_idx == m_GT_query_idx))
                      ? m_remapped_variant : variant;
  auto bcf_ht_type = BCF_INFO_GET_BCF_HT_TYPE(curr_tuple);
  //valid field handler
  assert(static_cast<size_t>(bcf_ht_type) < m_field_handlers.size()
         && m_field_handlers[bcf_ht_type].get());
  auto num_valid_input_elements = 0u;
  switch (BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(curr_tuple)) {
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_SUM:
      valid_result_found = m_field_handlers[bcf_ht_type]->get_valid_sum(src_variant, *m_query_config,
	  query_field_idx, result_ptr, num_valid_input_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MEAN:
      valid_result_found = m_field_handlers[bcf_ht_type]->get_valid_mean(src_variant, *m_query_config,
	  query_field_idx, result_ptr, num_valid_input_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MEDIAN:
      valid_result_found = m_field_handlers[bcf_ht_type]->get_valid_median(src_variant, *m_query_config,
	  query_field_idx, result_ptr, num_valid_input_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM:
      valid_result_found = remap_if_needed_and_combine(variant, &query_field_idx,
	  VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM,
	  const_cast<const void**>(&result_ptr), num_result_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_CONCATENATE:
      valid_result_found = m_field_handlers[bcf_ht_type]->concatenate_field(src_variant, *m_query_config,
	  query_field_idx, const_cast<const void**>(&result_ptr), num_result_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_HISTOGRAM_SUM:
      //do nothing, will be handled separately
      break;
    default:
      throw BroadCombinedGVCFException(std::string("Unknown VCF field combine operation ")
	  +std::to_string(BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(curr_tuple)));
      break;
  }
  return valid_result_found;
}

void BroadCombinedGVCFOperator::handle_INFO_fields(const Variant& variant) {
  //interval variant, add END tag
  if (m_remapped_variant.get_column_end() > m_remapped_variant.get_column_begin()) {
    int vcf_end_pos = m_remapped_variant.get_column_end() - m_curr_contig_begin_position + 1; //vcf END is 1 based
    bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "END", &vcf_end_pos, 1);
    m_bcf_record_size += sizeof(int);
  }
  for (auto i=0u; i<m_INFO_fields_vec.size(); ++i) {
    auto& curr_tuple = m_INFO_fields_vec[i];
    //Just need a 8-byte value, the contents could be a float or int (determined by the templated median function)
    int64_t result = 0;
    void* result_ptr = reinterpret_cast<void*>(&result);
    //For element wise operations
    auto num_result_elements = 1u;
    auto valid_result_found = handle_VCF_field_combine_operation(variant, curr_tuple, result_ptr, num_result_elements);
    if (valid_result_found) {
      auto bcf_ht_type = BCF_INFO_GET_BCF_HT_TYPE(curr_tuple);
      auto combined_result_bcf_ht_type = (bcf_ht_type == BCF_HT_INT
	  && BCF_INFO_GET_FIELD_INFO_PTR(curr_tuple)->is_VCF_field_combine_operation_sum()) ? BCF_HT_INT64
	: bcf_ht_type;
      if (BCF_INFO_GET_FIELD_INFO_PTR(curr_tuple)->m_length_descriptor.get_num_dimensions() == 1u) {
        bcf_update_info(m_vcf_hdr, m_bcf_out, BCF_INFO_GET_VCF_FIELD_NAME(curr_tuple).c_str(),
	    result_ptr, num_result_elements,
	    combined_result_bcf_ht_type);
        m_bcf_record_size += num_result_elements*VariantFieldTypeUtil::size(combined_result_bcf_ht_type);
      } else {
        auto stringified = std::move(m_field_handlers[bcf_ht_type]->stringify_2D_vector(*(BCF_INFO_GET_FIELD_INFO_PTR(curr_tuple))));
        bcf_update_info(m_vcf_hdr, m_bcf_out, BCF_INFO_GET_VCF_FIELD_NAME(curr_tuple).c_str(), stringified.c_str(), stringified.length(),
                        BCF_HT_STR);
        m_bcf_record_size += stringified.length();
      }
    }
  }
  for (auto& composite_vid_field_idx_INFO_histogram_tuple_pair : m_INFO_histogram_field_map) {
    const auto parent_field_vid_idx = composite_vid_field_idx_INFO_histogram_tuple_pair.first;
    auto& curr_INFO_histogram_tuple = composite_vid_field_idx_INFO_histogram_tuple_pair.second;
    unsigned query_field_idxs[2];
    query_field_idxs[0] = std::get<0>(curr_INFO_histogram_tuple);
    query_field_idxs[1] = std::get<1>(curr_INFO_histogram_tuple);
    unsigned num_result_elements = 0u;
    auto valid_result_found = remap_if_needed_and_combine(variant, query_field_idxs,
	  VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_HISTOGRAM_SUM,
	  0, num_result_elements);
    if (valid_result_found) {
      const auto& parent_vid_field_info = m_vid_mapper->get_field_info(parent_field_vid_idx);
      const auto& length_descriptor = parent_vid_field_info.m_length_descriptor;
      assert(length_descriptor.get_num_dimensions() == 2u);
      auto result_str = std::move(
	  GET_HISTOGRAM_FIELD_HANDLER_PTR_FROM_TUPLE(curr_INFO_histogram_tuple)->stringify_histogram(
	    length_descriptor.get_vcf_delimiter(0u),
	    length_descriptor.get_vcf_delimiter(1u)
	    ));
      bcf_update_info(m_vcf_hdr, m_bcf_out,
	  parent_vid_field_info.m_vcf_name.c_str(),
	  result_str.c_str(), result_str.length(), BCF_HT_STR);
      m_bcf_record_size += result_str.length();
    }
  }
}

void BroadCombinedGVCFOperator::handle_FORMAT_fields(const Variant& variant) {
  //For weird DP field handling
  auto valid_DP_found = false;
  auto valid_MIN_DP_found = false;
  m_MIN_DP_vector.resize(m_remapped_variant.get_num_calls()); //will do nothing after the first resize
  auto valid_DP_FORMAT_found = false;
  m_DP_FORMAT_vector.resize(m_remapped_variant.get_num_calls());
  //Pointer to extended vector inside field handler object
  const void* ptr = 0;
  uint64_t num_elements = 0ull;
  int* int_vec = 0;     //will point to same address as ptr, but as int*
  //Handle all fields - simply need to extend to the largest size
  for (auto i=0u; i<m_FORMAT_fields_vec.size(); ++i) {
    auto& curr_tuple = m_FORMAT_fields_vec[i];
    auto query_field_idx = BCF_FORMAT_GET_QUERY_FIELD_IDX(curr_tuple);
    auto length_descriptor = m_query_config->get_length_descriptor_for_query_attribute_idx(query_field_idx);
    //Fields such as PL are skipped if the #alleles is above a certain threshold
    if (length_descriptor.is_length_genotype_dependent()
        && too_many_alt_alleles_for_genotype_length_fields(m_merged_alt_alleles.size()))
      continue;
    auto known_field_enum = BCF_FORMAT_GET_KNOWN_FIELD_ENUM(curr_tuple);
    auto bcf_ht_type = BCF_FORMAT_GET_BCF_HT_TYPE(curr_tuple);
    auto is_char_type = (bcf_ht_type == BCF_HT_CHAR || bcf_ht_type == BCF_HT_STR);
    //valid field handler
    assert(static_cast<size_t>(bcf_ht_type) < m_field_handlers.size()
           && m_field_handlers[bcf_ht_type].get());
    //Check if this is a field that was remapped - for remapped fields, we must use field objects from m_remapped_variant
    //else we should use field objects from the original variant
    auto& src_variant = (m_remapping_needed && (length_descriptor.is_length_allele_dependent()
                         || query_field_idx == m_GT_query_idx))
                        ? m_remapped_variant : variant;
    auto valid_field_found = m_field_handlers[bcf_ht_type]->collect_and_extend_fields(src_variant, *m_query_config,
                             query_field_idx, &ptr, num_elements,
                             m_use_missing_values_not_vector_end && !is_char_type, m_use_missing_values_not_vector_end && is_char_type,
                             known_field_enum == GVCF_GT_IDX);
    if (valid_field_found) {
      auto j=0u;
      auto do_insert = !(m_query_config->sites_only_query());    //by default, insert into VCF record if !sites_only
      switch (known_field_enum) {
      case GVCF_GT_IDX: { //GT field is a pita
        int_vec = const_cast<int*>(reinterpret_cast<const int*>(ptr));
        auto num_elements_per_sample = num_elements/variant.get_num_calls();
        //GT of type 0/2 is stored as [0,0,2] in GenomicsDB if phased ploidy is used, else [0, 2]
        //GT of type 0|2 is stored as [0,1,2] in GenomicsDB if phased ploidy is used, else [0, 2]
        auto max_ploidy = length_descriptor.get_ploidy(num_elements_per_sample);
        uint64_t output_idx = 0ull;
        for (j=0ull; j<num_elements; j+=num_elements_per_sample)
          (*m_encode_GT_vector_function_ptr)(int_vec, j, num_elements_per_sample, output_idx);
        num_elements = max_ploidy*variant.get_num_calls();
        break;
      }
      case GVCF_GQ_IDX:
        do_insert = m_should_add_GQ_field;
        break;
      case GVCF_MIN_DP_IDX: //simply copy over min-dp values
        memcpy_s(&(m_MIN_DP_vector[0]), m_remapped_variant.get_num_calls()*sizeof(int),
                 ptr, m_remapped_variant.get_num_calls()*sizeof(int));
        valid_MIN_DP_found = true;
        break;
      case GVCF_DP_FORMAT_IDX:
        memcpy_s(&(m_DP_FORMAT_vector[0]), m_remapped_variant.get_num_calls()*sizeof(int),
                 ptr, m_remapped_variant.get_num_calls()*sizeof(int));
        valid_DP_FORMAT_found = true;
        do_insert = false; //Do not insert DP_FORMAT, wait till DP is read
        break;
      case GVCF_DP_IDX:
        int_vec = const_cast<int*>(reinterpret_cast<const int*>(ptr));
        valid_DP_found = true;
        do_insert = false; //Do not insert DP, handle DP related garbage at the end
        break;
      default:
        break;
      }
      if (do_insert) {
        bcf_update_format(m_vcf_hdr, m_bcf_out, BCF_FORMAT_GET_VCF_FIELD_NAME(curr_tuple).c_str(), ptr, num_elements,
                          bcf_ht_type);
        m_bcf_record_size += num_elements*VariantFieldTypeUtil::size(bcf_ht_type);
      }
    }
  }
  //Update DP fields
  if (valid_DP_found || valid_DP_FORMAT_found) {
    int64_t sum_INFO_DP = 0;
    auto found_one_valid_DP_FORMAT = false;
    for (auto j=0ull; j<m_remapped_variant.get_num_calls(); ++j) {
      int dp_info_val = valid_DP_found ? int_vec[j] : get_bcf_missing_value<int>();
      int dp_format_val = valid_DP_FORMAT_found ? m_DP_FORMAT_vector[j] : get_bcf_missing_value<int>();
      assert(is_bcf_valid_value<int>(dp_format_val) || dp_format_val == get_bcf_missing_value<int>());
      assert(is_bcf_valid_value<int>(dp_info_val) || dp_info_val == get_bcf_missing_value<int>());
      //If DP_FORMAT is invalid, use dp_info value for DP_FORMAT
      //dp_format_val = is_bcf_valid_value<int>(dp_format_val) ? dp_format_val : dp_info_val;
      if (!is_bcf_valid_value<int>(dp_info_val)) { //no valid DP info value found
        //MIN_DP gets higher priority
        if (valid_MIN_DP_found && is_bcf_valid_value<int>(m_MIN_DP_vector[j]))
          dp_info_val = m_MIN_DP_vector[j];
        else //else DP_FORMAT
          dp_info_val = dp_format_val;
      }
      m_DP_FORMAT_vector[j] = dp_format_val;
      found_one_valid_DP_FORMAT = is_bcf_valid_value<int>(dp_format_val) || found_one_valid_DP_FORMAT;
      sum_INFO_DP += (is_bcf_valid_value<int>(dp_info_val) ? dp_info_val : 0);
    }
    //Add DP to FORMAT if one valid DP FORMAT value found and this is not a sites only query
    if (found_one_valid_DP_FORMAT && !(m_query_config->sites_only_query())) {
      bcf_update_format_int32(m_vcf_hdr, m_bcf_out, "DP", &(m_DP_FORMAT_vector[0]), m_DP_FORMAT_vector.size()); //add DP FORMAT field
      m_bcf_record_size += m_DP_FORMAT_vector.size()*sizeof(int);
    }
    //If at least one valid DP value found from (DP or DP_FORMAT or MIN_DP), add DP to INFO
    if (sum_INFO_DP > 0 && !m_is_reference_block_only) {
      bcf_update_info(m_vcf_hdr, m_bcf_out, "DP", &sum_INFO_DP, 1, BCF_HT_INT64);
      m_bcf_record_size += sizeof(int64_t);
    }
  }
}

//FIXME: totally naive implementation - too many reallocations etc
void BroadCombinedGVCFOperator::merge_ID_field(const Variant& variant, const unsigned query_idx) {
#ifdef DEBUG
  //Produces ID fields in sorted order - consistency useful in CI testing
  //I'm sorry :(
  std::set<std::string> id_set;
#else
  std::unordered_set<std::string> id_set;
#endif
  for (const auto& curr_call : variant) {
    auto& field_ptr = curr_call.get_field(query_idx);
    if (field_ptr.get() && field_ptr->is_valid()) {
      auto* ptr = dynamic_cast<VariantFieldString*>(field_ptr.get());
      assert(ptr);
      const auto& curr_ID_value = ptr->get();
      auto last_begin_value = 0u;
      for (auto i=0u; i<curr_ID_value.length(); ++i)
        if (curr_ID_value[i] == ';') {
          id_set.insert(curr_ID_value.substr(last_begin_value, i-last_begin_value));
          last_begin_value = i+1u;
        }
      if (curr_ID_value.length() > last_begin_value)
        id_set.insert(curr_ID_value.substr(last_begin_value, curr_ID_value.length()-last_begin_value));
    }
  }
  m_ID_value.clear();
  for (const auto& str : id_set)
    m_ID_value += (str + ';');
  if (!m_ID_value.empty())
    m_ID_value.pop_back(); //delete last ';'
}

void BroadCombinedGVCFOperator::operate(Variant& variant, const VariantQueryConfig& query_config) {
#ifdef DO_PROFILING
  m_bcf_t_creation_timer.start();
#endif
  //Handle spanning deletions - change ALT alleles in calls with deletions to *, <NON_REF>
  handle_deletions(variant, query_config);
  GA4GHOperator::operate(variant, query_config);
  //Moved to new contig
  if (static_cast<int64_t>(m_remapped_variant.get_column_begin()) >= m_next_contig_begin_position) {
    std::string contig_name;
    int64_t contig_position;
    auto status = m_vid_mapper->get_contig_location(m_remapped_variant.get_column_begin(), contig_name, contig_position);
    if (status) {
      int64_t contig_begin_position = m_remapped_variant.get_column_begin() - contig_position;
      if (contig_begin_position != m_next_contig_begin_position) {
        m_next_contig_name = std::move(contig_name);
        m_next_contig_begin_position = contig_begin_position;
      }
    } else
      throw BroadCombinedGVCFException("Unknown contig for position "+std::to_string(m_remapped_variant.get_column_begin()));
    switch_contig();
  }
  //clear out
  bcf_clear(m_bcf_out);
  m_bcf_record_size = 0ull;
  //If no valid FORMAT fields exist, this value is never set
  m_bcf_out->n_sample = bcf_hdr_nsamples(m_vcf_hdr);
  //position
  m_bcf_out->rid = m_curr_contig_hdr_idx;
  m_bcf_out->pos = m_remapped_variant.get_column_begin() - m_curr_contig_begin_position;
  //ID field
  if (m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_ID_IDX)) {
    auto ID_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_ID_IDX);
    merge_ID_field(variant, ID_query_idx);
    if (!m_ID_value.empty())
      bcf_update_id(m_vcf_hdr, m_bcf_out, m_ID_value.c_str());
  }
  m_bcf_record_size += m_ID_value.length();
  //GATK combined GVCF does not care about QUAL value
  m_bcf_out->qual = get_bcf_missing_value<float>();
  if (BCF_INFO_GET_FIELD_INFO_PTR(m_vcf_qual_tuple)
      && m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_QUAL_IDX)
      && (BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(m_vcf_qual_tuple)
          != VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION)) {
    unsigned num_result_elements = 1u;
    auto qual_result = 1.0f;
    void* result_ptr = reinterpret_cast<void*>(&qual_result);
    auto valid_result_found = handle_VCF_field_combine_operation(variant, m_vcf_qual_tuple, result_ptr, num_result_elements);
    if (valid_result_found)
      m_bcf_out->qual = qual_result;
  }
  m_bcf_record_size += 3*sizeof(int);
  //Update alleles
  auto& ref_allele = dynamic_cast<VariantFieldString*>(m_remapped_variant.get_common_field(0u).get())->get();
  if (ref_allele.length() == 1u && ref_allele[0] == 'N') {
    ref_allele[0] = m_vcf_adapter->get_reference_base_at_position(m_curr_contig_name.c_str(), m_bcf_out->pos);
    if (BroadCombinedGVCFOperator::m_legal_bases.find(ref_allele[0]) == BroadCombinedGVCFOperator::m_legal_bases.end())
      ref_allele[0] = 'N';
  }
  const auto& alt_alleles = dynamic_cast<VariantFieldALTData*>(m_remapped_variant.get_common_field(1u).get())->get();
  auto total_num_merged_alleles = alt_alleles.size() + 1u;      //+1 for REF
  if (total_num_merged_alleles > m_alleles_pointer_buffer.size())
    m_alleles_pointer_buffer.resize(total_num_merged_alleles);
  //REF
  m_alleles_pointer_buffer[0] = ref_allele.c_str();
  m_bcf_record_size += ref_allele.length()*sizeof(char);
  //ALT
  for (auto i=1u; i<total_num_merged_alleles; ++i) {
    m_alleles_pointer_buffer[i] = alt_alleles[i-1u].c_str();
    m_bcf_record_size += alt_alleles[i-1u].length()*sizeof(char);
  }
  bcf_update_alleles(m_vcf_hdr, m_bcf_out, &(m_alleles_pointer_buffer[0]), total_num_merged_alleles);
  //FILTER fields
  if (m_query_config->produce_FILTER_field() &&
      m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_FILTER_IDX)) {
    auto FILTER_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_FILTER_IDX);
    //Remove duplicates across samples
    std::unordered_set<int> filter_idx_set;
    for (const auto& call : variant) {
      auto curr_FILTER_field = call.get_field<VariantFieldPrimitiveVectorData<int>>(FILTER_query_idx);
      if (curr_FILTER_field && curr_FILTER_field->is_valid()) {
        auto& curr_FILTER_idx_vec = curr_FILTER_field->get();
        filter_idx_set.insert(curr_FILTER_idx_vec.begin(), curr_FILTER_idx_vec.end());
      }
    }
    if (filter_idx_set.size()) {
      m_FILTER_idx_vec.resize(filter_idx_set.size());
      auto idx = 0u;
      for (auto global_field_idx : filter_idx_set) {
        assert(static_cast<size_t>(global_field_idx) < m_global_field_idx_to_hdr_idx.size());
        auto hdr_field_idx = m_global_field_idx_to_hdr_idx[global_field_idx];
        if (hdr_field_idx >= 0)
          m_FILTER_idx_vec[idx++] = hdr_field_idx;
      }
      bcf_update_filter(m_vcf_hdr, m_bcf_out, &(m_FILTER_idx_vec[0]), m_FILTER_idx_vec.size());
    }
  }
  //Flag that determines when to add GQ field - only when <NON_REF> is the only alternate allele
  //m_should_add_GQ_field = (m_NON_REF_exists && alt_alleles.size() == 1u);
  m_should_add_GQ_field = true; //always added in new version of CombineGVCFs
  //INFO fields
  handle_INFO_fields(variant);
  //FORMAT fields
  handle_FORMAT_fields(variant);
#ifdef DO_PROFILING
  m_bcf_t_creation_timer.stop();
#endif
#ifdef DO_MEMORY_PROFILING
  statm_t mem_result;
  read_off_memory_status(mem_result);
  if (mem_result.resident >= m_next_memory_limit) {
    std::cerr << "Crossed "<<m_next_memory_limit
              << " memory used (bytes) " << mem_result.resident
              <<" at contig "
              << bcf_hdr_id2name(m_vcf_hdr, m_curr_contig_hdr_idx)
              << " position " << m_bcf_out->pos+1
              << "\n";
    m_next_memory_limit = std::max<uint64_t>(m_next_memory_limit, mem_result.resident)
                          + 100*ONE_MB;
  }
#endif
  m_vcf_adapter->handoff_output_bcf_line(m_bcf_out, m_bcf_record_size);
}

void BroadCombinedGVCFOperator::switch_contig() {
  m_curr_contig_name = std::move(m_next_contig_name);
  m_curr_contig_begin_position = m_next_contig_begin_position;
  m_curr_contig_hdr_idx = bcf_hdr_id2int(m_vcf_hdr, BCF_DT_CTG, m_curr_contig_name.c_str());
  m_vid_mapper->get_next_contig_location(m_next_contig_begin_position, m_next_contig_name, m_next_contig_begin_position);
}

//Modifies original Variant object
void BroadCombinedGVCFOperator::handle_deletions(Variant& variant, const VariantQueryConfig& query_config) {
  //#merged alleles in spanning deletion can be at most 3 - REF,*,<NON_REF>; however, the #input alleles
  //can be much larger. LUT gets resized later
  m_reduced_alleles_LUT.resize_luts_if_needed(variant.get_num_calls(), 3u);
  m_reduced_alleles_LUT.reset_luts();
  auto GT_length_descriptor = m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX)
                              ? m_query_config->get_length_descriptor_for_query_attribute_idx(m_GT_query_idx)
                              : FieldLengthDescriptor();
  m_remapped_variant.resize(variant.get_num_calls(),  variant.get_call(0).get_num_fields());
  for (auto iter=variant.begin(), e=variant.end(); iter != e; ++iter) {
    auto& curr_call = *iter;
    auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
    //Deletion and not handled as spanning deletion
    //So replace the ALT with *,<NON_REF> and REF with "N"
    //Remap PL, AD fields
    if (curr_call.contains_deletion() && variant.get_column_begin() > curr_call.get_column_begin()) {
      auto& ref_allele = get_known_field<VariantFieldString, true>(curr_call, query_config, GVCF_REF_IDX)->get();
      auto& alt_alleles = get_known_field<VariantFieldALTData, true>(curr_call, query_config, GVCF_ALT_IDX)->get();
      assert(alt_alleles.size() > 0u);
      //Invalidate INFO fields
      for (const auto& tuple : m_INFO_fields_vec) {
        auto query_idx = BCF_INFO_GET_QUERY_FIELD_IDX(tuple);
        auto& field = curr_call.get_field(query_idx);
        if (field.get())
          field->set_valid(false);
      }
      //Already handled as a spanning deletion, nothing to do
      if (alt_alleles[0u] == g_vcf_SPANNING_DELETION &&
          (alt_alleles.size() == 1u || (alt_alleles.size() == 2u && IS_NON_REF_ALLELE(alt_alleles[1u]))))
        continue;
      m_reduced_alleles_LUT.resize_luts_if_needed(variant.get_num_calls(), alt_alleles.size()+1u); //+1 for REF
      //Reduced allele list will be REF="N", ALT="*, <NON_REF>"
      m_reduced_alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, 0, 0);  //REF-REF mapping
      //For each deletion allele, find the PL value corresponding to the genotype in which all the elements
      //of the genotype are equal to the deletion allele. The deletion allele with the lowest PL value is
      //mapped to "*" allele
      //GT field - for ploidy
      auto ploidy = 0u;
      auto* original_GT_field_ptr = (m_GT_query_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE)
                                    ? curr_call.get_field<VariantFieldPrimitiveVectorData<int>>(m_GT_query_idx) : 0;
      if (original_GT_field_ptr && original_GT_field_ptr->is_valid())
        ploidy = GT_length_descriptor.get_ploidy(original_GT_field_ptr->get().size());
      else
        original_GT_field_ptr = 0;
      auto lowest_deletion_allele_idx = -1;
      int lowest_PL_value = INT_MAX;
      auto PL_field_ptr = get_known_field_if_queried<VariantFieldPrimitiveVectorData<int>, true>(curr_call, query_config, GVCF_PL_IDX);
      auto has_NON_REF = false;
      //PL field exists flag
      auto PL_field_exists = (PL_field_ptr && PL_field_ptr->is_valid());
      std::vector<int> empty_int_vec;
      auto& PL_vector = PL_field_exists ? PL_field_ptr->get() : empty_int_vec;
      m_spanning_deletion_current_genotype.resize(ploidy);
      for (auto i=0u; i<alt_alleles.size(); ++i) {
        auto allele_idx = i+1;  //+1 for REF
        if (VariantUtils::is_deletion(ref_allele, alt_alleles[i])) {
          if (lowest_deletion_allele_idx < 0) //uninitialized, assign to first deletion found
            lowest_deletion_allele_idx = allele_idx;
          if (PL_field_exists) {
            //Genotype with all elements set to the deletion allele
            m_spanning_deletion_current_genotype.assign(ploidy, allele_idx);
            auto gt_idx = VariantOperations::get_genotype_index(m_spanning_deletion_current_genotype, true);
            //Truncations - dropped values etc
            if (gt_idx < PL_vector.size() && PL_vector[gt_idx] < lowest_PL_value) {
              lowest_PL_value = PL_vector[gt_idx];
              lowest_deletion_allele_idx = allele_idx;
            }
          }
        } else if (IS_NON_REF_ALLELE(alt_alleles[i])) {
          m_reduced_alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, allele_idx, 2);
          has_NON_REF = true;
        }
      }
      assert(lowest_deletion_allele_idx >= 1);    //should be an ALT allele
      //first ALT allele in reduced list is *
      m_reduced_alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, lowest_deletion_allele_idx, 1);
      if (has_NON_REF) {
        alt_alleles.resize(2u);
        alt_alleles[1u] = TILEDB_NON_REF_VARIANT_REPRESENTATION;
      } else
        alt_alleles.resize(1u); //only spanning deletion
      ref_allele = "N"; //set to unknown REF for now
      alt_alleles[0u] = g_vcf_SPANNING_DELETION;
      unsigned num_reduced_alleles = alt_alleles.size() + 1u;   //+1 for REF
      //Remap fields that need to be remapped
      for (auto i=0u; i<m_remapped_fields_query_idxs.size(); ++i) {
        auto query_field_idx = m_remapped_fields_query_idxs[i];
        auto length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(query_field_idx);
        //field whose length is dependent on #alleles
        assert(length_descriptor.is_length_allele_dependent());
        auto& curr_field = curr_call.get_field(query_field_idx);
        if (curr_field.get() && curr_field->is_valid()) {   //Not null
          if (length_descriptor.get_num_dimensions() > 1u) {
            auto& remapped_field =
              m_remapped_variant.get_call(curr_call_idx_in_variant).get_field(query_field_idx);
	    copy_field(remapped_field, curr_call.get_field(query_field_idx));
            remap_allele_specific_annotations(curr_field,
                                              remapped_field,
                                              curr_call_idx_in_variant,
                                              m_reduced_alleles_LUT, num_reduced_alleles, has_NON_REF, ploidy,
                                              query_config, query_field_idx);
            std::swap(curr_field, remapped_field);
          } else {
            unsigned num_reduced_elements =
              length_descriptor.get_num_elements(num_reduced_alleles-1u, ploidy, 0u);     //#alt alleles
            //Remapper for variant
            RemappedVariant remapper_variant(variant, query_field_idx);
            //Copy field to pass to remap function
            assert(i < m_spanning_deletions_remapped_fields.size());
            copy_field(m_spanning_deletions_remapped_fields[i], curr_field);
            curr_field->resize(num_reduced_elements);
            //Get handler for current type
            auto& handler = get_handler_for_type(query_config.get_element_type(query_field_idx));
            assert(handler.get());
            //Call remap function
            handler->remap_vector_data(
              m_spanning_deletions_remapped_fields[i], curr_call_idx_in_variant,
              m_reduced_alleles_LUT, num_reduced_alleles, has_NON_REF, ploidy,
              length_descriptor, num_reduced_elements, remapper_variant);
          }
        }
      }
      //GT field
      if (original_GT_field_ptr) {
        auto& input_GT =
          original_GT_field_ptr->get();
        m_spanning_deletion_remapped_GT.resize(input_GT.size());
        auto null_PL_ptr = std::move(std::unique_ptr<VariantFieldBase>(nullptr));
        auto remap_GT_based_on_input_GT =
          update_GT_to_correspond_to_min_PL_value(
            query_config,
            query_config.is_defined_query_idx_for_known_field_enum(GVCF_PL_IDX)
            ? curr_call.get_field(query_config.get_query_idx_for_known_field_enum(GVCF_PL_IDX))
            : null_PL_ptr,
            input_GT,
            GT_length_descriptor,
            num_reduced_alleles,
            has_NON_REF
          );
        if (remap_GT_based_on_input_GT) { //update_GT_to_correspond_to_min_PL_value didn't work, remap based on input_GT only
          VariantOperations::remap_GT_field(input_GT, m_spanning_deletion_remapped_GT, m_reduced_alleles_LUT, curr_call_idx_in_variant,
                                            num_reduced_alleles, has_NON_REF, GT_length_descriptor);
          //Copy back
          memcpy_s(&(input_GT[0]), input_GT.size()*sizeof(int), &(m_spanning_deletion_remapped_GT[0]), input_GT.size()*sizeof(int));
        }
      }
    }
  }
}

bool BroadCombinedGVCFOperator::update_GT_to_correspond_to_min_PL_value(
  const VariantQueryConfig& query_config,
  std::unique_ptr<VariantFieldBase>& PL_field,
  std::vector<int>& input_GT,
  const FieldLengthDescriptor& GT_length_descriptor,
  const unsigned num_alleles,
  const bool has_NON_REF) {
  auto remap_GT_based_on_input_GT = true;
  auto PL_field_ptr = PL_field.get();
  if (PL_field_ptr && PL_field_ptr->is_valid()
      && m_query_config->produce_GT_with_min_PL_value_for_spanning_deletions()) {
    remap_GT_based_on_input_GT = false;
    //Get handler for current type
    auto& handler = get_handler_for_type(query_config.get_element_type(m_GT_query_idx));
    assert(handler.get());
    //Get the tuple containing data for min value of PL
    auto min_value_tuple = handler->determine_allele_combination_and_genotype_index_for_min_value(
                             PL_field, num_alleles, has_NON_REF,
                             GT_length_descriptor.get_ploidy(input_GT.size()));
    //Found one valid PL value - use allele idx vec
    if (GenotypeForMinValueTracker<int>::found_at_least_one_valid_value(min_value_tuple)) {
      //if phased ploidy, allele idx and phase information alternate in input_GT vector
      //So, step value will be 2, if no phase, then step == 1
      auto step_value = GT_length_descriptor.get_ploidy_step_value();
      auto& allele_idx_vec = GenotypeForMinValueTracker<int>::get_allele_idx_vec(min_value_tuple);
      for (auto i=0u,j=0u; i<input_GT.size(); i+=step_value,++j) {
        assert(j < allele_idx_vec.size());
        input_GT[i] =  allele_idx_vec[j];
      }
    } else
      remap_GT_based_on_input_GT = true; //no valid PL values found, remap based on input GT
  }
  return remap_GT_based_on_input_GT;
}

#endif //ifdef HTSDIR
