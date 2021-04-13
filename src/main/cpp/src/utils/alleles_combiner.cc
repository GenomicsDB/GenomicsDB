/**
 * The MIT License (MIT)
 * Copyright (c) 2020-21 Omics Data Automation, Inc.
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

#include "genomicsdb_iterators.h"
#include "alleles_combiner.h"
#include "known_field_info.h"

//Convenience class grouping variables that're useful while parsing ALT alleles for samples
class AlleleConfig {
  public:
    AlleleConfig() {
      begin_allele(0, 0);
    }
    void begin_allele(const unsigned arg_allele_idx, const char* arg_ptr) {
      is_symbolic_allele = false;
      is_deletion = false;
      is_MNV = false;
      is_NON_REF_allele = false;
      num_mismatched_bases = 0u;
      length = 0u;
      ptr = arg_ptr;
      allele_idx = arg_allele_idx;
    }
    bool is_symbolic_allele;
    bool is_deletion;
    bool is_MNV;
    bool is_NON_REF_allele;
    unsigned num_mismatched_bases;
    unsigned allele_idx;
    size_t length;
    const char* ptr;
};

AllelesCombiner::AllelesCombiner(const GenomicsDBGVCFIterator* iterator, const size_t num_queried_rows) {
  m_iterator = iterator;
  m_contains_deletion.resize(num_queried_rows, false);
  m_contains_MNV.resize(num_queried_rows, false);
  m_is_REF_block.resize(num_queried_rows, false);
  m_contains_deletion_or_MNV_spanning_current_location.resize(num_queried_rows, false);
  m_allele_idx_for_deletion_or_MNV_spanning_current_location.resize(num_queried_rows, UNDEFINED_ATTRIBUTE_IDX_VALUE);
  m_NON_REF_idx_vec.resize(num_queried_rows, UNDEFINED_ATTRIBUTE_IDX_VALUE);
#ifdef DEBUG //in DEBUG mode, catch out of bounds errors (if bugs)
  m_alleles_LUT.resize_luts_if_needed(num_queried_rows, 0u);
#else
  m_row_query_idx_with_deletion_MNV_at_current_location_vec.reserve(num_queried_rows);
  m_index_in_deletions_MNVs_vec.reserve(num_queried_rows);
  m_deletion_MNV_allele_idx_vec.reserve(num_queried_rows); //no particular reason for using this size
  m_alleles_LUT.resize_luts_if_needed(num_queried_rows, 100u);
  m_merged_alleles_vec.reserve(100u); //no particular reason for using 100, most locations should have <100 alleles
#endif
  reset_for_next_query_interval();
}

void AllelesCombiner::reset_before_adding_new_sample_info_at_current_position() {
  m_merged_alleles_vec.resize(1u, MergedAllelesVecEntry(0u, STRING_VIEW(0, 0), false)); //1 for REF
  m_merged_alleles_vec[0u].REF_length = 0u;
  m_merged_alleles_vec[0u].allele = STRING_VIEW(0, 0u);
  m_spanning_deletion_allele_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  //This should cover all single-base alleles, '&' is GenomicsDB's representation for NON_REF
  for(const auto b : "atgcnATGCN*&")
    m_single_base_ALT_allele_to_merged_idx[static_cast<uint8_t>(b)] = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  m_single_base_REF_ALT_allele_to_index.clear();
  m_deletion_REF_length_to_index.clear();
  m_MNV_ALT_to_index.clear();
  //Not resetting full m_alleles_LUT here since we can do a reset only for
  //rows newly added in this iteration of the GVCF iterator. Saves some memory accesses
  //Mark rows that had deletions/MNVs starting at previous locations as spanning
  for(auto i=0ull;i<m_row_query_idx_with_deletion_MNV_at_current_location_vec.size();++i) {
    const auto row_query_idx = m_row_query_idx_with_deletion_MNV_at_current_location_vec[i];
    if(m_iterator && m_iterator->is_valid_row_query_idx(row_query_idx)) {
      m_contains_deletion_or_MNV_spanning_current_location[row_query_idx] = true;
      assert(m_index_in_deletions_MNVs_vec[i+1u] > m_index_in_deletions_MNVs_vec[i]); //at least 1 deletion/MNV
      const auto index = m_index_in_deletions_MNVs_vec[i];
      //FIXME: this should be computed based on PL values. Functionality will be added later
      m_allele_idx_for_deletion_or_MNV_spanning_current_location[row_query_idx]
        = m_deletion_MNV_allele_idx_vec[index];
    }
  }
  m_row_query_idx_with_deletion_MNV_at_current_location_vec.clear();
  m_index_in_deletions_MNVs_vec.resize(1u);
  m_index_in_deletions_MNVs_vec[0u] = 0u;
  m_deletion_MNV_allele_idx_vec.clear();
}

void AllelesCombiner::reset_for_next_query_interval() {
  m_num_calls_with_deletions_or_MNVs = 0u;
  m_num_calls_with_NON_REF_allele = 0u;
  reset_before_adding_new_sample_info_at_current_position();
}

std::pair<unsigned, bool> AllelesCombiner::handle_single_base_ALT(const char base) {
  const auto merged_allele_idx = m_single_base_ALT_allele_to_merged_idx[static_cast<uint8_t>(base)];
  if(merged_allele_idx == UNDEFINED_ATTRIBUTE_IDX_VALUE) {
    const auto num_merged_alleles = m_merged_alleles_vec.size();
    m_single_base_ALT_allele_to_merged_idx[static_cast<uint8_t>(base)] = num_merged_alleles;
    return std::pair<unsigned, bool>(num_merged_alleles, true);
  }
  return std::pair<unsigned, bool>(merged_allele_idx, false);
}

std::pair<unsigned, bool> AllelesCombiner::handle_insertion_or_symbolic_allele(const STRING_VIEW& ALT) {
  //insert if doesn't exist
  auto result_pair = m_single_base_REF_ALT_allele_to_index.insert(
      std::pair<STRING_VIEW, unsigned>(ALT, m_merged_alleles_vec.size()));
  return std::pair<unsigned, bool>((*(result_pair.first)).second, result_pair.second);
}

void AllelesCombiner::handle_allele(const size_t row_query_idx,  const STRING_VIEW& orig_REF, AlleleConfig& allele_config) {
  if(allele_config.length == 0u || allele_config.is_NON_REF_allele) //NON REF should be the last allele - don't add it to map or LUT
    return;
  STRING_VIEW REF(orig_REF);
  STRING_VIEW ALT(allele_config.ptr, allele_config.length);
  //result_pair.first will contain index in m_merged_alleles_vec where the current ALT allele is stored
  //result_pair.second is true if this is a new ALT allele, false otherwise
  std::pair<unsigned, bool> result_pair; 
  const auto is_REF_single_base = (REF.length() == 1u);
  if(is_REF_single_base || allele_config.is_symbolic_allele) { //treat the REF of symbolic alleles as single base
    REF = STRING_VIEW(REF.data(), 1u);
    //REF and alleles are already 'normalized' since single base REF or symbolic
    result_pair = (ALT.length() == 1u) ? handle_single_base_ALT(ALT[0u])
      : handle_insertion_or_symbolic_allele(ALT);
  }
  else { //multi-base REF and non-symbolic ALT allele
    //Example   REF=ATGC  ALT=A,AGC,TTGC,TAGC
    //The first 2 are deletions, the 3rd a SNV and the 4th an MNV
    //We store 'normalized' versions of these alleles by removing the longest common suffix
    //from the REF and ALT. So the normalized versions would be:
    //ATGC->A, AT->A, A->T, AT->TA.
    //Normalized versions are uniquely defined and independent of the other ALT alleles
    const auto same_length = (REF.length() == ALT.length());
    //SNV
    if(same_length && allele_config.num_mismatched_bases == 1u) {
      assert(ALT[0u] != REF[0u]);
      REF = STRING_VIEW(REF.data(), 1u);
      ALT = STRING_VIEW(ALT.data(), 1u);
      result_pair = handle_single_base_ALT(ALT[0u]);
    }
    else {
      if(ALT.length() > REF.length()) { //insertion
        //Normalized form of insertion has 1 base REF and 1+insertion length ALT
        //So longest common suffix would be REF.length-1
        ALT = STRING_VIEW(ALT.data(), ALT.length() - (REF.length()-1u));
        REF = STRING_VIEW(REF.data(), 1u);
        result_pair = handle_insertion_or_symbolic_allele(ALT);
      }
      else {
        allele_config.is_deletion = (ALT.length() < REF.length()); //'*' handled above
        allele_config.is_MNV = (same_length && allele_config.num_mismatched_bases > 1u);
        assert(allele_config.is_deletion || allele_config.is_MNV);
        if(allele_config.is_deletion) {
          //In a normalized deletion, ALT has a single base. So longest common suffix
          //length is ALT.length-1
          REF = STRING_VIEW(REF.data(), REF.length()-(ALT.length()-1u));
          ALT = STRING_VIEW(ALT.data(), 1u);
          //Since ALT.length is 1, we can distinguish between different deletions
          //by just noting the length of the REF allele in the normalized representation
          //Example TAGC,->T TA->T, TAG->T
          //So m_deletion_REF_length_to_index is simply a map<REF_length, index>
          auto insert_result_pair = m_deletion_REF_length_to_index.insert(
              std::pair<unsigned, unsigned>(REF.length(), m_merged_alleles_vec.size()));
          result_pair = std::pair<unsigned, bool>((*(insert_result_pair.first)).second,
              insert_result_pair.second);
        }
        else { //MNV
          const auto suffix_length = VariantUtils::find_longest_common_suffix_length(REF, ALT);
          REF = STRING_VIEW(REF.data(), REF.length()-suffix_length);
          ALT = STRING_VIEW(ALT.data(), ALT.length()-suffix_length);
          auto insert_result_pair = m_MNV_ALT_to_index.insert(
              std::pair<STRING_VIEW,unsigned>(ALT, m_merged_alleles_vec.size()));
          result_pair = std::pair<unsigned, bool>((*(insert_result_pair.first)).second,
              insert_result_pair.second);
        }
      }
    }
  }
  if(result_pair.second) {//new alt allele
    //index of spanning deletion in merged vec
    if(allele_config.is_symbolic_allele && allele_config.is_deletion)
      m_spanning_deletion_allele_idx = result_pair.first;
    m_merged_alleles_vec.push_back(MergedAllelesVecEntry(REF.length(), ALT, allele_config.is_symbolic_allele));
    m_alleles_LUT.resize_luts_if_needed(m_merged_alleles_vec.size());
  }
  m_alleles_LUT.add_input_merged_idx_pair(row_query_idx, allele_config.allele_idx, result_pair.first);
  //Track indexes of deletions and MNVs
  if(allele_config.is_deletion || allele_config.is_MNV)
    m_deletion_MNV_allele_idx_vec.push_back(allele_config.allele_idx);
}

void AllelesCombiner::insert_allele_info(const size_t row_query_idx, const STRING_VIEW& REF,
    const STRING_VIEW& delimited_ALT_str, const bool begins_before_curr_start_position) {
  assert(row_query_idx < m_contains_MNV.size());
  assert(m_merged_alleles_vec.size() > 0u); //element 0 for REF allele
  //Updated merged REF allele - longest REF
  if(!begins_before_curr_start_position && m_merged_alleles_vec[0].REF_length < REF.length()) {
    m_merged_alleles_vec[0u].allele = REF;
    m_merged_alleles_vec[0u].REF_length = REF.length();
  }
  //New set of alleles - reset LUT for this row
  m_alleles_LUT.reset_luts_for_sample(row_query_idx);
  m_alleles_LUT.add_input_merged_idx_pair(row_query_idx, 0u, 0u); //for REF
  //Properties for this sample (over all alleles) 
  auto contains_deletion = false;
  auto contains_MNV = false;
  auto num_deletions_or_MNVs = 0u;
  unsigned NON_REF_allele_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  //split up ALT - delimited by |
  const char* ptr = delimited_ALT_str.data();
  //Per allele properties
  AlleleConfig allele_config;
  auto num_alleles = 1u; //1 since 0 is REF
  allele_config.begin_allele(num_alleles, ptr);
  for(auto i=0u;i<delimited_ALT_str.length();++i) {
    const auto curr_char = ptr[i];
    switch(curr_char) {
      case TILEDB_NON_REF_VARIANT_REPRESENTATION[0u]:
        allele_config.is_symbolic_allele = true;
        allele_config.is_NON_REF_allele = true;
        NON_REF_allele_idx = allele_config.allele_idx;
        break;
      case TILEDB_ALT_ALLELE_SEPARATOR[0u]:
        {
          handle_allele(row_query_idx, REF, allele_config);
          contains_deletion = contains_deletion || allele_config.is_deletion;
          contains_MNV = contains_MNV || allele_config.is_MNV;
          num_deletions_or_MNVs += (allele_config.is_deletion || allele_config.is_MNV);
          ++num_alleles;
          //Reset per allele values
          allele_config.begin_allele(num_alleles, ptr+i+1u);
          allele_config.length = -1; //increment at the end of the loop body fixes this
          break;
        }
      case '<':
      case '>':
      case '[':
      case ']':
        allele_config.is_symbolic_allele = true;
        break;
      case '*':
        allele_config.is_symbolic_allele = true;
        allele_config.is_deletion = true;
        break;
      default:
        allele_config.num_mismatched_bases += (allele_config.length < REF.length()
            && curr_char != REF[allele_config.length]);
        break;
    }
    ++(allele_config.length);
  }
  //For last ALT allele which doesn't have a separator
  if(allele_config.length) {
    handle_allele(row_query_idx, REF, allele_config);
    contains_deletion = contains_deletion || allele_config.is_deletion;
    contains_MNV = contains_MNV || allele_config.is_MNV;
    num_deletions_or_MNVs += (allele_config.is_deletion || allele_config.is_MNV);
    ++num_alleles;
  }
  //Update flags and indexes
  m_contains_deletion[row_query_idx] = contains_deletion;
  m_contains_MNV[row_query_idx] = contains_MNV;
  m_is_REF_block[row_query_idx] = (num_alleles == 2u && NON_REF_allele_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE);
  //since this row has an entry at current location, it cannot be a spanning deletion that begins before current location
  m_contains_deletion_or_MNV_spanning_current_location[row_query_idx] = false;
  m_NON_REF_idx_vec[row_query_idx] = NON_REF_allele_idx;
  m_num_calls_with_deletions_or_MNVs += (contains_deletion || contains_MNV) ;
  m_num_calls_with_NON_REF_allele += (NON_REF_allele_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE);
  if(num_deletions_or_MNVs) {
    m_row_query_idx_with_deletion_MNV_at_current_location_vec.push_back(row_query_idx);
    m_index_in_deletions_MNVs_vec.push_back(m_index_in_deletions_MNVs_vec.back()+num_deletions_or_MNVs);
  }
}

void AllelesCombiner::finished_updating_allele_info_for_current_position() {
  //If there are deletions/MNVs that begin before current position and no samples beginning at the
  //current position have "*" in ALT list, add "*" to merged list
  assert(m_num_calls_with_deletions_or_MNVs >=
      m_row_query_idx_with_deletion_MNV_at_current_location_vec.size());
  const auto num_calls_with_deletions_or_MNVs_that_begin_before_current_position
    = m_num_calls_with_deletions_or_MNVs
    - m_row_query_idx_with_deletion_MNV_at_current_location_vec.size();
  if(num_calls_with_deletions_or_MNVs_that_begin_before_current_position > 0u) {
    if(m_spanning_deletion_allele_idx == UNDEFINED_ATTRIBUTE_IDX_VALUE) {
      m_spanning_deletion_allele_idx = m_merged_alleles_vec.size();
      m_merged_alleles_vec.emplace_back(1u, STRING_VIEW("*", 1u), true);
    }
    if(m_iterator) {
      //Update lut for spanning deletions/MNVs
      for(auto iter=m_iterator->begin_valid_row_query_idx();iter!=m_iterator->end_valid_row_query_idx();++iter) {
        const auto row_query_idx = *iter;
        if(m_contains_deletion_or_MNV_spanning_current_location[row_query_idx]) {
          m_alleles_LUT.reset_luts_for_sample(row_query_idx);
          m_alleles_LUT.add_input_merged_idx_pair(row_query_idx, 0u, 0u); //for REF
          m_alleles_LUT.add_input_merged_idx_pair(row_query_idx,
              m_allele_idx_for_deletion_or_MNV_spanning_current_location[row_query_idx],
              m_spanning_deletion_allele_idx);
          //not adding NON_REF - will be handled by contains_NON_REF_allele()
        }
      }
    }
  }
  //Add NON_REF as the last allele in the merged vec
  if(m_num_calls_with_NON_REF_allele > 0u)
    m_merged_alleles_vec.emplace_back(1u, STRING_VIEW("<NON_REF>", 9u), true);
}

void AllelesCombiner::remove_allele_info(const size_t row_query_idx) {
  m_num_calls_with_deletions_or_MNVs -= contains_deletion_or_MNV(row_query_idx);
  m_num_calls_with_NON_REF_allele -= contains_NON_REF_allele(row_query_idx); 
}

void AllelesCombiner::get_merged_VCF_spec_alleles_vec(std::string& buffer, std::vector<STRING_VIEW>& alleles_vec) const {
  const auto& REF = m_merged_alleles_vec[0u].allele;
  alleles_vec.push_back(REF);
  //First reserve capacity in buffer before setting pointers in alleles_vec - else dangling references
  auto num_chars_to_add = 0ull;
  for(auto i=1u;i<m_merged_alleles_vec.size();++i) {
    const auto& entry = m_merged_alleles_vec[i];
    assert(entry.is_symbolic_allele || REF.length() >= entry.REF_length);
    //add suffix
    num_chars_to_add += entry.allele.length() + (entry.is_symbolic_allele ? 0u : REF.length()-entry.REF_length);
  }
  buffer.reserve(buffer.size()+num_chars_to_add);
  for(auto i=1u;i<m_merged_alleles_vec.size();++i) {
    const auto& entry = m_merged_alleles_vec[i];
    const auto begin_idx = buffer.size();
    buffer.append(entry.allele.data(), entry.allele.length());
    //denormalize ALT by adding suffix
    if(!entry.is_symbolic_allele) {
      assert(REF.length() >= entry.REF_length);
      //suffix - everything in REF after entry.REF
      buffer.append(REF.data()+entry.REF_length, REF.length()-entry.REF_length);
    }
    alleles_vec.emplace_back(&(buffer[begin_idx]), buffer.size()-begin_idx);
  }
}
