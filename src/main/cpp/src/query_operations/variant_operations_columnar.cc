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

void SingleVariantOperatorBase::operate_on_columnar_cell(const GenomicsDBGVCFCell& variant) {
  m_iterator = variant.get_iterator();
  assert(m_iterator);
  auto curr_interval = m_iterator->get_column_range();
  //New contig - get contig info from vid
  if(!(m_contig_info_ptr && curr_interval.first >= m_contig_info_ptr->m_tiledb_column_offset
        && curr_interval.second < m_contig_info_ptr->m_tiledb_column_offset+m_contig_info_ptr->m_length)) {
    auto status = m_vid_mapper->get_contig_info_for_location(curr_interval.first, m_contig_info_ptr);
    if(!status)
      throw VidMapperException(std::string("Could not find contig for column ")+std::to_string(curr_interval.first));
  }
}

void BroadCombinedGVCFOperator::operate_on_columnar_cell(const GenomicsDBGVCFCell& variant) {
  SingleVariantOperatorBase::operate_on_columnar_cell(variant); 
}

template<class WriterTy>
bool BroadCombinedGVCFOperator::write_vcf_line(const GenomicsDBGVCFCell& variant) {
  static_assert(std::is_base_of<AbstractVCFWriterBase, WriterTy>::value);
  auto* writer = static_cast<WriterTy*>(m_writer_base_ptr);
  //The moment no_overflow becomes false, stop writing
  //By writing code in the following way, we avoid tedious if-else statements in the code
  //and let the compiler handle it
  auto no_overflow = true;
  //contig and position
  auto curr_interval = m_iterator->get_column_range();
  auto contig_position = curr_interval.first - m_contig_info_ptr->m_tiledb_column_offset;
  auto contig_position_end = curr_interval.second - m_contig_info_ptr->m_tiledb_column_offset;
  no_overflow = no_overflow
    && writer-> template write<char, true>(m_contig_info_ptr->m_name.c_str(), m_contig_info_ptr->m_name.length());
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  no_overflow = no_overflow && writer-> template write<int64_t, true>(contig_position+1); //VCF is based
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //ID
  no_overflow = no_overflow && writer-> template write<char, true>('.');
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //REF
  no_overflow = no_overflow
    && writer-> template write<char, true>(m_merged_reference_allele.c_str(), m_merged_reference_allele.length());
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //ALT
  if(m_merged_alt_alleles.size() == 0u)
    no_overflow = no_overflow && writer-> template write<char, true>('.');
  else {
    no_overflow = no_overflow && writer-> template write<char, true>(m_merged_alt_alleles[0u].data(),
          m_merged_alt_alleles[0u].length());
    for(auto i=1u;i<m_merged_alt_alleles.size();++i) {
      no_overflow = no_overflow && writer-> template write<char, true>(',');
      no_overflow = no_overflow && writer-> template write<char, true>(m_merged_alt_alleles[i].data(),
          m_merged_alt_alleles[i].length());
    }
  }
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //QUAL
  no_overflow = no_overflow && writer-> template write<char, true>('.');
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //FILTER
  no_overflow = no_overflow && writer-> template write<char, true>('.');
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //INFO
  if(contig_position_end > contig_position) {
    no_overflow = no_overflow && writer-> template write<char, true>("END=");
    no_overflow = no_overflow && writer-> template write<int, true>(contig_position_end+1);
  }
  else
    no_overflow = no_overflow && writer-> template write<char, true>('.');
  no_overflow = no_overflow && writer-> template write<char, true>('\t');
  //FORMAT fields
  for(const auto& field_info : m_FORMAT_fields_vec) {
  }
  return no_overflow;
}
