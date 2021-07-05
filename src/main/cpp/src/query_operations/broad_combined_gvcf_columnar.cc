/**
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation Inc 
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

#include "variant_operations.h"
#include "vcf.h"
#include "broad_combined_gvcf.h"
#include "gt_remapper.h"

template<class WriterTy, bool contains_phase, bool produce_GT_field, bool do_remap>
bool BroadCombinedGVCFOperator::write_vcf_line(WriterTy& writer, const GenomicsDBGVCFCell& variant) {
  //The moment no_overflow becomes false, stop writing
  //By writing code in the following way, we avoid tedious if-else statements in the code
  //and let the compiler handle it
  auto no_overflow = true;
  //contig and position
  auto curr_interval = m_iterator->get_current_variant_interval();
  assert(m_contig_info_ptr && m_contig_info_ptr->m_tiledb_column_offset <= curr_interval.first && 
      curr_interval.second < m_contig_info_ptr->m_tiledb_column_offset + m_contig_info_ptr->m_length);
  auto contig_position = curr_interval.first - m_contig_info_ptr->m_tiledb_column_offset;
  auto contig_position_end = curr_interval.second - m_contig_info_ptr->m_tiledb_column_offset;
  no_overflow = no_overflow
    && writer. template write<char, true>(m_contig_info_ptr->m_name.c_str(), m_contig_info_ptr->m_name.length());
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  no_overflow = no_overflow && writer. template write<int64_t, true>(contig_position+1); //VCF is based
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  //ID
  no_overflow = no_overflow && writer. template write<char, true>('.');
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  const auto& alleles_combiner = variant.get_iterator()->get_alleles_combiner();
  const auto& norm_REF_ALT_vec = alleles_combiner.get_merged_normalized_REF_ALT_vec();
  assert(norm_REF_ALT_vec.size()); //at least REF
  const auto& REF = norm_REF_ALT_vec[0].allele;
  //REF
  if(REF.length() == 0u) //interval cut, get base from ref genome
    no_overflow = no_overflow && writer. template write<char, true>(
        m_vcf_adapter->get_reference_base_at_position(m_contig_info_ptr->m_name.c_str(), contig_position));
  else
    no_overflow = no_overflow
      && writer. template write<char, true>(REF.data(), REF.length());
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  //ALT
  if(norm_REF_ALT_vec.size() == 1u)
    no_overflow = no_overflow && writer. template write<char, true>('.');
  else {
    const auto& entry = norm_REF_ALT_vec[1u];
    const auto& allele = entry.allele;
    no_overflow = no_overflow && writer. template write<char, true>(allele.data(), allele.length());
    //add suffix
    if(!(entry.is_symbolic_allele)) {
      assert(REF.length() >= entry.REF_length);
      no_overflow = no_overflow && writer. template write<char, true>(REF.data()+entry.REF_length, REF.length()-entry.REF_length);
    }
    for(auto i=2u;i<norm_REF_ALT_vec.size();++i) {
      no_overflow = no_overflow && writer. template write<char, true>(',');
      const auto& entry = norm_REF_ALT_vec[i];
      const auto& allele = entry.allele;
      no_overflow = no_overflow && writer. template write<char, true>(allele.data(), allele.length());
      //add suffix
      if(!(entry.is_symbolic_allele)) {
        assert(REF.length() >= entry.REF_length);
        no_overflow = no_overflow && writer. template write<char, true>(REF.data()+entry.REF_length, REF.length()-entry.REF_length);
      }
    }
  }
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  //QUAL
  no_overflow = no_overflow && writer. template write<char, true>('.');
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  //FILTER
  no_overflow = no_overflow && writer. template write<char, true>('.');
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  //INFO
  if(contig_position_end > contig_position) {
    no_overflow = no_overflow && writer. template write<char, true>(static_cast<const char*>("END="), 4u);
    no_overflow = no_overflow && writer. template write<int, true>(contig_position_end+1);
  }
  else
    no_overflow = no_overflow && writer. template write<char, true>('.');
  if(!no_overflow)
    return false;
  //FORMAT fields
  if(m_query_config->sites_only_query()) {
    no_overflow = no_overflow && writer. template write<char, true>('\n');
    return no_overflow;
  }
  //GT only for now
  if(!(m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX))) {
    no_overflow = no_overflow && writer. template write<char, true>('\n');
    return no_overflow;
  }
  no_overflow = no_overflow && writer. template write<char, true>('\t');
  no_overflow = no_overflow && writer. template write<char, true>(static_cast<const char*>("GT"), 2u);
  const auto num_queried_rows = m_query_config->get_num_rows_to_query();
  const auto& gt_remapper = m_iterator->get_GT_remapper();
  for(auto row_query_idx=0ull;row_query_idx<num_queried_rows;++row_query_idx) {
    no_overflow = no_overflow && writer. template write<char, true>('\t');
    if(!(m_iterator->is_valid_row_query_idx(row_query_idx))) {
      no_overflow = no_overflow && writer. template write<char, true>('.');
      continue;
    }
    assert(m_iterator->is_valid_row_query_idx(row_query_idx));
    no_overflow = no_overflow
      && gt_remapper.remap_for_row_query_idx<WriterTy, contains_phase, produce_GT_field, do_remap>(
          writer, row_query_idx);
  }
  no_overflow = no_overflow && writer. template write<char, true>('\n');
  for(const auto& field_info : m_FORMAT_fields_vec) {
  }
  return no_overflow;
}

void BroadCombinedGVCFOperator::operate_on_columnar_cell(const GenomicsDBGVCFCell& variant) {
  SingleVariantOperatorBase::operate_on_columnar_cell(variant);
  // The objective is to minimize the number of if-else conditions at runtime in the lower levels of the
  // callgraph. For example, the condition produce_GT_field is determined at the time of query and doesn't
  // change through the lifetime of the query. Hence, we would like to check for it at the highest possible
  // level of the callgraph and execute only code that's required for the specific value of produce_GT_field
  // used in this query. Using C++ template parameters, this becomes a lot easier than manually defining N versions
  // of a function differing in their implementation based on the value of produce_GT_field (and other parameters).
  // By using template parameters (a) we let the compiler generate the N versions needed (b) using if statements with
  // template parameters in function bodies allows us to write a single function body and let the compiler optimize
  // each version based on the value of the template parameter (c) pass template parameters to lower level functions
  // We have a big switch near the top of the callgraph to choose the right callgraph version to invoke (this switch
  // statement was produced using a script)
  const bool contains_phase = m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX)
    && m_query_config->get_length_descriptor_for_query_attribute_idx(m_GT_query_idx).contains_phase_information();
  const bool produce_GT_field = m_query_config->produce_GT_field();
  const bool do_remap = m_iterator->get_alleles_combiner().remapping_needed();
  //Get an integer representation of the parameters - if C++ had tuple based switch statements, this would
  //have been easier
  switch((static_cast<uint64_t>(m_writer_type_enum) << 3u)
      | (static_cast<uint64_t>(contains_phase) << 2u)
      | (static_cast<uint64_t>(produce_GT_field) << 1u)
      | (static_cast<uint64_t>(do_remap))
      )
  {
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, true, true, true>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, true, true, false>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, true, false, true>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, true, false, false>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, false, true, true>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, false, true, false>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, false, false, true>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_STRING_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::string>, false, false, false>(m_vcf_writer_to_string, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, true, true, true>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, true, true, false>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, true, false, true>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(true) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, true, false, false>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, false, true, true>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(true) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, false, true, false>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(true)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, false, false, true>(m_vcf_writer_to_ostream, variant);
    break;
    case ( (static_cast<uint64_t>(VCFWRITER_ENUM::STL_OSTREAM_NO_LIMIT) << 3u) 
        | (static_cast<uint64_t>(false) << 2u)
        | (static_cast<uint64_t>(false) << 1u)
        | (static_cast<uint64_t>(false)) ):
      write_vcf_line<VCFWriterNoOverflow<std::ostream>, false, false, false>(m_vcf_writer_to_ostream, variant);
    break;
  }
}
