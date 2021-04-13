#ifndef GT_REMAPPER_H
#define GT_REMAPPER_H 1

#include "headers.h"

class GenomicsDBGVCFIterator;
class AllelesCombiner;

//Exceptions thrown
class GTRemapperException : public std::exception {
 public:
  GTRemapperException(const std::string m="") : msg_("GTRemapperException exception : "+m) { ; }
  ~GTRemapperException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};


/*
 * GT remapping module for columnar gvcf iterator
 */
class GTRemapper {
  public:
    GTRemapper(const unsigned GT_query_idx, const GenomicsDBGVCFIterator& iter);
    /*
     * OperatorTy should define two functions
     * write_GT_allele_index(const int gt);
     * write_GT_phase(const int v); // v is 0 for unphased, 1 for phased
     * remap_for_row_query_idx will call these functions for each element of the ploidy/phase
     * in order
     */
    template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
    bool remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const;
    template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
    bool remap_all_queried_valid_rows(OperatorTy& op) const;
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
  private:
    unsigned m_GT_query_idx;
    const GenomicsDBGVCFIterator* m_iterator;
    const AllelesCombiner* m_alleles_combiner;
};

#endif
