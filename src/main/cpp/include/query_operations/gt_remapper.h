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

#ifndef GT_REMAPPER_H
#define GT_REMAPPER_H 1

#include "headers.h"

template<typename T>
class AllelesCombiner;

//Exceptions thrown
class GTRemapperException : public std::exception {
 public:
  GTRemapperException(const std::string m="") : msg_("GTRemapperException exception : "+m) { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};


/*
 * GT remapping module for columnar gvcf iterator - provide the allele index values as
 * per the remapping computed in AllelesCombiner
 */

// The template parameter ValidRowAndGTDataProviderTy should define three functions:
// bool is_valid_row_query_idx(uint64_t row_query_idx)
// begin_valid_row_query_idx(), end_valid_row_query_idx() - this should return a forward iterator that iterates
// over valid row query idxs. Dereferencing the iterator should give you row_query_idx (uint64_t)
// std::pair<uint8_t*, size_t> get_raw_pointer_and_length_for_query_idx(uint64_t row_query_idx, const int field_query_idx)
// which returns the pointer containing GT data
// GenomicsDBGVCFIterator is an example of a class that satisfies the requirements of ValidRowAndGTDataProviderTy
// Why template, why not just use GenomicsDBGVCFIterator? Makes unit testing difficult since I have to deal
// with a complex GenomicsDBGVCFIterator object. I can define simple classes in the unit tests for ValidRowAndGTDataProviderTy

template<typename ValidRowAndGTDataProviderTy>
class GTRemapper {
  public:
    GTRemapper(const unsigned GT_query_idx, const ValidRowAndGTDataProviderTy& iter,
        const AllelesCombiner<ValidRowAndGTDataProviderTy>& alleles_combiner);
  private:
    /*
     * Wrapper around AllelesCombiner get_merged_allele_idx with an extra template parameter
     * specific to the GT field. Same principle of avoiding runtime if conditions in inner
     * loops
     * produce_GT_field - should match value of the field produce_GT_field in m_query_config
     */
    template<bool produce_GT_field, bool do_remap, bool is_REF_block, bool contains_NON_REF_allele>
    int get_merged_allele_idx(const size_t row_query_idx, const int allele_idx) const;

    template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap,
      bool is_REF_block, bool contains_NON_REF_allele>
    bool remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const;
  public:
    /*
     * OperatorTy should define three functions:
     * bool write_GT_allele_index(const uint64_t row_query_idx, const int gt);
     * bool write_GT_phase(const uint64_t row_query_idx, const int v); // v is 0 for unphased, 1 for phased
     * bool write_GT_empty(const uint64_t row_query_idx) //when GT field is missing for row_query_idx
     * remap_for_row_query_idx will call these functions for each element of the ploidy/phase
     * in order
     * The functions should return false if there is an overflow or the operator couldn't accept the data
     */
    template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
    bool remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const;

    template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
    bool remap_all_queried_valid_rows(OperatorTy& op) const;

  private:
    unsigned m_GT_query_idx;
    const ValidRowAndGTDataProviderTy* m_provider;
    const AllelesCombiner<ValidRowAndGTDataProviderTy>* m_alleles_combiner;
};

#endif
