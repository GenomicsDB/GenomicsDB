/**
 * The MIT License (MIT)
 * Copyright (c) 2022 Omics Data Automation, Inc.
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
#ifndef PL_REMAPPER_TEMPLATE_DEFINITION_HPP
#define PL_REMAPPER_TEMPLATE_DEFINITION_HPP 1

#include "genomicsdb_logger.h"
#include "known_field_info.h"
#include "pl_remapper.h"
#include "vcf.h"

template <typename ValidRowAndGTDataProviderTy>
PLRemapper<ValidRowAndGTDataProviderTy>::PLRemapper(
    const unsigned GT_query_idx, const ValidRowAndGTDataProviderTy& iter,
    const AllelesCombiner<ValidRowAndGTDataProviderTy>& alleles_combiner)
    : m_GT_query_idx(GT_query_idx), m_provider(&iter), m_alleles_combiner(&alleles_combiner) {}

template <typename ValidRowAndGTDataProviderTy>
template <typename OperatorTy, typename element_ty, PloidyEnum ploidy_enum, bool do_remap, bool is_REF_block,
          bool contains_NON_REF_allele>
bool PLRemapper<ValidRowAndGTDataProviderTy>::remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx,
                                                                      unsigned ploidy) {
  const auto PL_ptr_length_pair = m_provider->get_raw_pointer_and_length_for_query_idx(row_query_idx, m_PL_query_idx);
  auto length = PL_ptr_length_pair.second;
  auto PL_field_ptr = reinterpret_cast<const element_ty*>(PL_ptr_length_pair.first);
  auto no_overflow = true;
  if (length == 0u || ploidy == 0u) {
    // PL field missing, upto the operator to decide what to do. VCFWriter will put a '.' or just blank string
    no_overflow = no_overflow && op.template write_empty<element_ty>(m_PL_query_idx, row_query_idx);
    return no_overflow;
  }
  if (!do_remap) {
    no_overflow =
        no_overflow && op.template write_element<element_ty, true>(m_PL_query_idx, row_query_idx, PL_field_ptr[0u]);
    for (auto i = 1u; i < length; ++i)
      no_overflow =
          no_overflow && op.template write_element<element_ty, false>(m_PL_query_idx, row_query_idx, PL_field_ptr[i]);
    return no_overflow;
  }
  auto num_merged_alleles = m_alleles_combiner->get_num_merged_alleles();
  m_merged_to_input_allele_idx.resize(num_merged_alleles);
  for (auto i = 0u; i < num_merged_alleles; ++i)
    m_merged_to_input_allele_idx[i] =
        m_alleles_combiner->template get_input_allele_idx<do_remap, is_REF_block, contains_NON_REF_allele>(
            row_query_idx, i);
  auto num_gts = KnownFieldInfo::get_number_of_genotypes(num_merged_alleles - 1u, ploidy);  //-1 for #ALT alleles
  // First element - all elements of ploidy == REF
  auto input_gt_idx = 0;
  if (static_cast<unsigned>(input_gt_idx) >= length)
    no_overflow = no_overflow && op.template write_missing<element_ty, true>(m_PL_query_idx, row_query_idx);
  else
    no_overflow = no_overflow && op.template write_element<element_ty, true>(m_PL_query_idx, row_query_idx,
                                                                             PL_field_ptr[input_gt_idx]);
  // FIXME: genotype count limits
  switch (ploidy_enum) {
    case PloidyEnum::HAPLOID:  // haploid, number of genotypes == number of alleles
    {
      for (auto i = 1u; i < num_gts; ++i) {
        auto input_gt_idx = m_merged_to_input_allele_idx[i];
        if (input_gt_idx < 0 || static_cast<unsigned>(input_gt_idx) >= length)
          no_overflow = no_overflow && op.template write_missing<element_ty, false>(m_PL_query_idx, row_query_idx);
        else
          no_overflow = no_overflow && op.template write_element<element_ty, false>(m_PL_query_idx, row_query_idx,
                                                                                    PL_field_ptr[input_gt_idx]);
      }
      break;
    }
    case PloidyEnum::DIPLOID: {
      for (auto allele_j = 0u; allele_j < num_merged_alleles; ++allele_j) {
        auto input_allele_j = m_merged_to_input_allele_idx[allele_j];
        if (input_allele_j < 0) {  // merged allele has no counterpart in input (including NON_REF)
          for (auto allele_k = (allele_j == 0u) ? 1u : allele_j; allele_k < num_merged_alleles; ++allele_k) {
            no_overflow = no_overflow && op.template write_missing<element_ty, false>(m_PL_query_idx, row_query_idx);
          }
          continue;
        }
        // First element of PL vector has been dealt outside the switch statement
        for (auto allele_k = (allele_j == 0u) ? 1u : allele_j; allele_k < num_merged_alleles; ++allele_k) {
          auto input_allele_k = m_merged_to_input_allele_idx[allele_k];
          auto input_gt_idx = bcf_alleles2gt(input_allele_j, input_allele_k);
          if (input_gt_idx < 0 || static_cast<unsigned>(input_gt_idx) >= length)
            no_overflow = no_overflow && op.template write_missing<element_ty, false>(m_PL_query_idx, row_query_idx);
          else
            no_overflow = no_overflow && op.template write_element<element_ty, false>(m_PL_query_idx, row_query_idx,
                                                                                      PL_field_ptr[input_gt_idx]);
        }
      }
      break;
    }
    case PloidyEnum::GENERAL: {
#define SET_PLOIDY_INDEX_IN_STACK_ELEMENT(X, Y) ((X).first = (Y))
#define SET_ALLELE_INDEX_IN_STACK_ELEMENT(X, Y) ((X).second = (Y))
#define GET_PLOIDY_INDEX_IN_STACK_ELEMENT(X) ((X).first)
#define GET_ALLELE_INDEX_IN_STACK_ELEMENT(X) ((X).second)
      m_merged_allele_idx_vec_for_current_gt_combination.resize(
          ploidy + 1u);  //+1 to avoid unnecessary if statements in the while loop
      m_input_call_allele_idx_vec_for_current_gt_combination.resize(ploidy);
      // Enumerate genotypes based on method described in
      // http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
      // Use custom "stack" instead of a recursive function call
      //"top" of the stack is the last element of the vector
      // resize to max #genotypes - avoids frequent dynamic memory allocations/frees
      m_ploidy_index_allele_index_stack.resize(num_gts);
      // In each iteration, generate all genotypes where the last ploidy corresponds to the allele
      // corresponding to top_level_allele_idx
      auto allele_idx = 0;
      auto ploidy_idx = 0;
      // In each iteration of the loop, let the top of the stack contain (P=x, A=y)
      // The subsequent iterations will generate all genotype combinations corresponding to
      // gt[x] == y first before moving on to elements in the stack below
      // Initializer element in the stack set (P=ploidy, A=#alleles-1)
      // Thus, the while loop will generate all genotype combinations for ploidies 0..ploidy-1
      // with alleles 0..alleles-1
      SET_PLOIDY_INDEX_IN_STACK_ELEMENT(m_ploidy_index_allele_index_stack[0u], ploidy);
      SET_ALLELE_INDEX_IN_STACK_ELEMENT(m_ploidy_index_allele_index_stack[0u], num_merged_alleles - 1);
      auto num_elements_in_stack = 1u;
      auto remapped_gt_idx = 0ull;
      while (num_elements_in_stack > 0u) {
        auto& top_stack_element = m_ploidy_index_allele_index_stack[num_elements_in_stack - 1u];
        allele_idx = GET_ALLELE_INDEX_IN_STACK_ELEMENT(top_stack_element);
        ploidy_idx = GET_PLOIDY_INDEX_IN_STACK_ELEMENT(top_stack_element);
        --num_elements_in_stack;  // popped stack
        assert(ploidy_idx >= 0 &&
               static_cast<size_t>(ploidy_idx) < m_merged_allele_idx_vec_for_current_gt_combination.size());
        m_merged_allele_idx_vec_for_current_gt_combination[ploidy_idx] = allele_idx;
        // Assigned one allele idx for all ploidys
        if (ploidy_idx == 0) {
          if (remapped_gt_idx > 0u)  // first element (all REF) has been handled outside the switch
          {
            auto curr_genotype_combination_contains_missing_allele_for_input = false;
            for (auto i = 0u; i < ploidy; ++i) {
              auto input_allele_idx =
                  m_merged_to_input_allele_idx[m_merged_allele_idx_vec_for_current_gt_combination[i]];
              m_input_call_allele_idx_vec_for_current_gt_combination[i] = input_allele_idx;
              if (input_allele_idx < 0)  // no mapping found for current allele in input gvcf
                curr_genotype_combination_contains_missing_allele_for_input = true;
            }
            auto input_gt_idx =
                KnownFieldInfo::get_genotype_index(m_input_call_allele_idx_vec_for_current_gt_combination, false);
            if (curr_genotype_combination_contains_missing_allele_for_input || input_gt_idx >= length)
              no_overflow = no_overflow && op.template write_missing<element_ty, false>(m_PL_query_idx, row_query_idx);
            else
              no_overflow = no_overflow && op.template write_element<element_ty, false>(m_PL_query_idx, row_query_idx,
                                                                                        PL_field_ptr[input_gt_idx]);
          }
          ++remapped_gt_idx;
        } else {
          --ploidy_idx;  // current ploidy_idx
          // Reverse order so that alleles with lower idx are closer to the top of the stack
          for (auto i = allele_idx; i >= 0; --i) {
            assert(num_elements_in_stack < m_ploidy_index_allele_index_stack.size());
            auto& curr_stack_element = m_ploidy_index_allele_index_stack[num_elements_in_stack];
            SET_PLOIDY_INDEX_IN_STACK_ELEMENT(curr_stack_element, ploidy_idx);
            SET_ALLELE_INDEX_IN_STACK_ELEMENT(curr_stack_element, i);
            ++num_elements_in_stack;
          }
        }
      }
      break;
    }
    default:
      throw PLRemapperException(std::string("Unhandled ploidy type enum ") + std::to_string(ploidy_enum));
  }
  return no_overflow;
}

template <typename ValidRowAndGTDataProviderTy>
template <typename OperatorTy, typename element_ty, bool contains_phase, bool do_remap, bool is_REF_block,
          bool contains_NON_REF_allele>
bool PLRemapper<ValidRowAndGTDataProviderTy>::remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) {
  assert(m_provider->is_valid_row_query_idx(row_query_idx));
  const auto GT_ptr_length_pair = m_provider->get_raw_pointer_and_length_for_query_idx(row_query_idx, m_GT_query_idx);
  auto ploidy = KnownFieldInfo::get_ploidy<contains_phase>(GT_ptr_length_pair.second);
  switch (ploidy) {
    case 0u:
      return remap_for_row_query_idx<OperatorTy, element_ty, PloidyEnum::NONE, do_remap, is_REF_block,
                                     contains_NON_REF_allele>(op, row_query_idx, ploidy);
    case 1u:
      return remap_for_row_query_idx<OperatorTy, element_ty, PloidyEnum::HAPLOID, do_remap, is_REF_block,
                                     contains_NON_REF_allele>(op, row_query_idx, ploidy);
    case 2u:
      return remap_for_row_query_idx<OperatorTy, element_ty, PloidyEnum::DIPLOID, do_remap, is_REF_block,
                                     contains_NON_REF_allele>(op, row_query_idx, ploidy);
    default:
      return remap_for_row_query_idx<OperatorTy, element_ty, PloidyEnum::GENERAL, do_remap, is_REF_block,
                                     contains_NON_REF_allele>(op, row_query_idx, ploidy);
  }
}

template <typename ValidRowAndGTDataProviderTy>
template <typename OperatorTy, typename element_ty, bool contains_phase, bool do_remap>
bool PLRemapper<ValidRowAndGTDataProviderTy>::remap_for_row_query_idx(OperatorTy& op, unsigned PL_query_idx,
                                                                      const size_t row_query_idx) {
  m_PL_query_idx = PL_query_idx;
  assert(m_provider->is_valid_row_query_idx(row_query_idx));
  // Push conditions as high as possible to void runtime if conditions in low level functions
  switch ((m_alleles_combiner->is_REF_block(row_query_idx) << 1u) |
          (m_alleles_combiner->contains_NON_REF_allele(row_query_idx))) {
    case 0u:  // neither REF block nor contains NON_REF
      return remap_for_row_query_idx<OperatorTy, element_ty, contains_phase, do_remap, false, false>(op, row_query_idx);
    case 1u:  // no REF block, but contains NON_REF
      return remap_for_row_query_idx<OperatorTy, element_ty, contains_phase, do_remap, false, true>(op, row_query_idx);
    case 3u:  // REF block and contains NON_REF
      return remap_for_row_query_idx<OperatorTy, element_ty, contains_phase, do_remap, true, true>(op, row_query_idx);
    default:  // illegal
    {
      auto msg =
          std::string("Is REF block but doesn't contain valid NON_REF allele index ") + std::to_string(row_query_idx);
      logger.fatal(PLRemapperException(msg), msg.c_str());
      break;
    }
  }
  return false;
}

template <typename ValidRowAndGTDataProviderTy>
template <typename OperatorTy, typename element_ty, bool contains_phase, bool do_remap>
bool PLRemapper<ValidRowAndGTDataProviderTy>::remap_all_queried_valid_rows(OperatorTy& op, unsigned PL_query_idx) {
  auto no_overflow = true;
  // Iterate over valid rows
  for (auto iter = m_provider->begin_valid_row_query_idx(); iter != m_provider->end_valid_row_query_idx(); ++iter) {
    no_overflow = no_overflow &&
                  remap_for_row_query_idx<OperatorTy, element_ty, contains_phase, do_remap>(op, PL_query_idx, *iter);
    if (!no_overflow) break;
  }
  return no_overflow;
}

#endif
