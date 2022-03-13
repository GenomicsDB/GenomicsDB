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

#ifndef PL_REMAPPER_H
#define PL_REMAPPER_H 1

#include "headers.h"

template <typename T>
class AllelesCombiner;

// Exceptions thrown
class PLRemapperException : public std::exception {
 public:
  PLRemapperException(const std::string m = "") : msg_("PLRemapperException exception : " + m) { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept { return msg_.c_str(); }

 private:
  std::string msg_;
};

/*
 * PL remapping module for columnar gvcf iterator - remaps PL like fields based on order of alleles using
 * AllelesCombiner ie fields whose length is equal to the number of genotypes
 */

/*
 The template parameter ValidRowAndGTDataProviderTy should define three functions:
 bool is_valid_row_query_idx(uint64_t row_query_idx)
 begin_valid_row_query_idx(), end_valid_row_query_idx() - this should return a forward iterator that iterates
 over valid row query idxs. Dereferencing the iterator should give you row_query_idx (uint64_t)
 std::pair<uint8_t*, size_t> get_raw_pointer_and_length_for_query_idx(uint64_t row_query_idx, const int field_query_idx)
 which returns the pointer containing GT data
 GenomicsDBGVCFIterator is an example of a class that satisfies the requirements of ValidRowAndGTDataProviderTy
 Why template, why not just use GenomicsDBGVCFIterator? Makes unit testing difficult since I have to deal
 with a complex GenomicsDBGVCFIterator object. I can define simple classes in the unit tests for
 ValidRowAndGTDataProviderTy
*/

enum PloidyEnum : unsigned {
  NONE = 0u,
  HAPLOID = 1u,
  DIPLOID = 2u,
  GENERAL  // every other ploidy
};

template <typename ValidRowAndGTDataProviderTy>
class PLRemapper {
 public:
  PLRemapper(unsigned GT_query_idx, const ValidRowAndGTDataProviderTy& iter,
             const AllelesCombiner<ValidRowAndGTDataProviderTy>& alleles_combiner);

 private:
  template <typename OperatorTy, typename element_ty, PloidyEnum ploidy_enum, bool do_remap, bool is_REF_block,
            bool contains_NON_REF_allele>
  bool remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx, unsigned ploidy);

  template <typename OperatorTy, typename element_ty, bool contains_phase, bool do_remap, bool is_REF_block,
            bool contains_NON_REF_allele>
  bool remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx);

 public:
  /*
   * OperatorTy should define a function:
   * template<typename element_ty>
   * bool write_empty(unsigned PL_query_idx, const uint64_t row_query_idx) - called when the whole field is empty
   * template<typename element_ty, bool is_first_element>
   * bool write_element<element_ty, is_first_element>(unsigned PL_query_idx, const uint64_t row_query_idx,
   * const element_ty v);
   * template<typename element_ty, bool is_first_element>
   * bool write_missing<is_first_element>(unsigned PL_query_idx, const uint64_t row_query_idx) - missing element
   * remap_for_row_query_idx will call these functions for each element of the remapped vector
   * The function should return false if there is an overflow or the operator couldn't accept the data
   */

  // PL_query_idx is not necessarily PL, just the query idx of a field whose length is equal to the number
  // of genotypes which we want to remap

  template <typename OperatorTy, typename element_ty, bool contains_phase, bool do_remap>
  bool remap_for_row_query_idx(OperatorTy& op, unsigned PL_query_idx, const size_t row_query_idx);

  template <typename OperatorTy, typename element_ty, bool contains_phase, bool do_remap>
  bool remap_all_queried_valid_rows(OperatorTy& op, unsigned PL_query_idx);

 private:
  // PL_query_idx is not necessarily PL, just the field idx which we want to remap
  unsigned m_PL_query_idx;
  unsigned m_GT_query_idx;
  const ValidRowAndGTDataProviderTy* m_provider;
  const AllelesCombiner<ValidRowAndGTDataProviderTy>* m_alleles_combiner;
  std::vector<int> m_merged_to_input_allele_idx{128u};  // for a single row, tmp variable for opt
  // Vector to store allele indexes for a given genotype index - avoid dynamic memory ops
  // One for merged, other for input
  // Sized for ploidy 32
  std::vector<int> m_merged_allele_idx_vec_for_current_gt_combination{32u};
  std::vector<int> m_input_call_allele_idx_vec_for_current_gt_combination{32u};
  // Ploidy-allele index stack to mimic recursive function behavior - save mallocs
  std::vector<std::pair<int, int>> m_ploidy_index_allele_index_stack{128u};
};

#endif
