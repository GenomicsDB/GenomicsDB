/**
 * The MIT License (MIT)
 * Copyright (c) 2021 Omics Data Automation, Inc.
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
#ifndef GT_REMAPPER_TEMPLATE_DEFINITION_HPP
#define GT_REMAPPER_TEMPLATE_DEFINITION_HPP 1

#include "gt_remapper.h"
#include "vcf.h"

template<typename ValidRowAndGTDataProviderTy>
GTRemapper<ValidRowAndGTDataProviderTy>::GTRemapper(const unsigned GT_query_idx,
    const ValidRowAndGTDataProviderTy& iter,
    const AllelesCombiner<ValidRowAndGTDataProviderTy>& alleles_combiner)
  : m_GT_query_idx(GT_query_idx),
    m_provider(&iter),
    m_alleles_combiner(&alleles_combiner) {
}

template<typename ValidRowAndGTDataProviderTy> template<bool produce_GT_field, bool do_remap, bool is_REF_block, bool contains_NON_REF_allele>
int GTRemapper<ValidRowAndGTDataProviderTy>::get_merged_allele_idx(const size_t row_query_idx, const int allele_idx) const {
  assert(m_provider->is_valid_row_query_idx(row_query_idx));
  if(!produce_GT_field) //static condition
    return get_bcf_gt_no_call_allele_index<int>();
  else
    return (do_remap ?
        m_alleles_combiner->template get_merged_allele_idx<do_remap, is_REF_block, contains_NON_REF_allele>(row_query_idx, allele_idx)
        : allele_idx);
}

template<typename ValidRowAndGTDataProviderTy> template< typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap,
  bool is_REF_block, bool contains_NON_REF_allele>
bool GTRemapper<ValidRowAndGTDataProviderTy>::remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const {
  assert(m_provider->is_valid_row_query_idx(row_query_idx));
  const auto ptr_length_pair = m_provider->get_raw_pointer_and_length_for_query_idx(row_query_idx, m_GT_query_idx);
  const auto GT_ptr = reinterpret_cast<const int*>(ptr_length_pair.first);
  const auto length = ptr_length_pair.second;
  auto no_overflow = true;
  if(length) {
    no_overflow = no_overflow && op.write_GT_allele_index(row_query_idx,
        get_merged_allele_idx<produce_GT_field, do_remap, is_REF_block, contains_NON_REF_allele>(
          row_query_idx, GT_ptr[0u]));
    const auto step_value = contains_phase ? 2u : 1u; //decided statically
    for(auto i=1u;i<length;i+=step_value) {
      no_overflow = no_overflow && op.write_GT_phase(row_query_idx, contains_phase ? GT_ptr[i] : 0); //decided statically
      no_overflow = no_overflow && op.write_GT_allele_index(row_query_idx,
          get_merged_allele_idx<produce_GT_field, do_remap, is_REF_block, contains_NON_REF_allele>(
            row_query_idx, GT_ptr[i + contains_phase]));
    }
  }
  else {
    //GT field missing, upto the operator to decide what to do. VCFWriter will put a '.'
    no_overflow = no_overflow && op.write_GT_empty(row_query_idx);
  }
  return no_overflow;
}

template<typename ValidRowAndGTDataProviderTy> template< typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
bool GTRemapper<ValidRowAndGTDataProviderTy>::remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const {
  assert(m_provider->is_valid_row_query_idx(row_query_idx));
  //Push conditions as high as possible to void runtime if conditions in low level functions
  switch((m_alleles_combiner->is_REF_block(row_query_idx) << 1u)
      | (m_alleles_combiner->contains_NON_REF_allele(row_query_idx))) {
    case 0u: //neither REF block nor contains NON_REF
      return remap_for_row_query_idx<OperatorTy, contains_phase, produce_GT_field, do_remap, false, false>(op, row_query_idx);
      break;
    case 1u: //no REF block, but contains NON_REF
      return remap_for_row_query_idx<OperatorTy, contains_phase, produce_GT_field, do_remap, false, true>(op, row_query_idx);
      break;
    case 3u: //REF block and contains NON_REF
      return remap_for_row_query_idx<OperatorTy, contains_phase, produce_GT_field, do_remap, true, true>(op, row_query_idx);
      break;
    default: //illegal
      throw GTRemapperException(std::string("Is REF block but doesn't contain valid NON_REF allele index ")+std::to_string(row_query_idx));
      //throw GTRemapperException(std::string("Both REF block and spanning deletion enabled for row query idx ")
      //+ std::to_string(row_query_idx));
      break;
  }
  return false;
}

template<typename ValidRowAndGTDataProviderTy>
template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
bool GTRemapper<ValidRowAndGTDataProviderTy>::remap_all_queried_valid_rows(OperatorTy& op) const {
  auto no_overflow = true;
  //Iterate over valid rows
  for(auto iter=m_provider->begin_valid_row_query_idx();iter!=m_provider->end_valid_row_query_idx();++iter)
    no_overflow = no_overflow && remap_for_row_query_idx<OperatorTy, contains_phase, produce_GT_field, do_remap>(op, *iter);
  return no_overflow;
}

#endif
