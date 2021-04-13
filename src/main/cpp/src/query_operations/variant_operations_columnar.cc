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
