#include "gt_remapper.h"
#include "vcf_fmt_writer.h"
#include "genomicsdb_iterators.h"

GTRemapper::GTRemapper(const unsigned GT_query_idx,
    const GenomicsDBGVCFIterator& iter)
  : m_GT_query_idx(GT_query_idx),
    m_iterator(&iter),
    m_alleles_combiner(&(iter.get_alleles_combiner())) {
}

template<bool produce_GT_field, bool do_remap, bool is_REF_block, bool contains_NON_REF_allele>
inline int GTRemapper::get_merged_allele_idx(const size_t row_query_idx, const int allele_idx) const {
  assert(m_iterator->is_valid_row_query_idx(row_query_idx));
  if(!produce_GT_field) //static condition
    return get_bcf_gt_no_call_allele_index<int>();
  else
    return m_alleles_combiner->get_merged_allele_idx<do_remap, is_REF_block, contains_NON_REF_allele>(row_query_idx, allele_idx);
}

template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap,
  bool is_REF_block, bool contains_NON_REF_allele>
bool GTRemapper::remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const {
  assert(m_iterator->is_valid_row_query_idx(row_query_idx));
  const auto ptr_length_pair = m_iterator->get_raw_pointer_and_length_for_query_idx(row_query_idx, m_GT_query_idx);
  const auto GT_ptr = reinterpret_cast<const int*>(ptr_length_pair.first);
  const auto length = ptr_length_pair.second;
  auto no_overflow = true;
  if(length)
    no_overflow = no_overflow && op.write_GT_allele_index(
        get_merged_allele_idx<produce_GT_field, do_remap, is_REF_block, contains_NON_REF_allele>(
          row_query_idx, GT_ptr[0u]));
  const auto step_value = contains_phase ? 2u : 1u; //decided statically
  for(auto i=1u;i<length;i+=step_value) {
    no_overflow = no_overflow && op.write_GT_phase(contains_phase ? GT_ptr[i] : 0); //decided statically
    no_overflow = no_overflow && op.write_GT_allele_index(
        get_merged_allele_idx<produce_GT_field, do_remap, is_REF_block, contains_NON_REF_allele>(
        row_query_idx, GT_ptr[i + contains_phase]));
  }
  return no_overflow;
}

template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
bool GTRemapper::remap_for_row_query_idx(OperatorTy& op, const size_t row_query_idx) const {
  assert(m_iterator->is_valid_row_query_idx(row_query_idx));
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

template<typename OperatorTy, bool contains_phase, bool produce_GT_field, bool do_remap>
bool GTRemapper::remap_all_queried_valid_rows(OperatorTy& op) const {
  auto no_overflow = true;
  //Iterate over valid rows
  for(auto iter=m_iterator->begin_valid_row_query_idx();iter!=m_iterator->end_valid_row_query_idx();++iter)
    no_overflow = no_overflow && remap_for_row_query_idx<OperatorTy, contains_phase, produce_GT_field, do_remap>(op, *iter);
  return no_overflow;
}

//Explicit template instantiations
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, true, true, true>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, true, true, false>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, true, false, true>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, true, false, false>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, false, true, true>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, false, true, false>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, false, false, true>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterFSB, false, false, false>(VCFWriterFSB& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, true, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, true, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, false, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, false, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, true, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, true, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, false, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, false, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, true, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, true, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, false, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, false, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, true, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, true, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, false, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, false, false>(VCFWriterNoOverflow<std::ostream>& op) const;
