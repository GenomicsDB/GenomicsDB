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

#include "gt_remapper_template_definition.h"
#include "vcf_fmt_writer.h"
#include "genomicsdb_iterators.h"

//Explicit template instantiations
template class GTRemapper<GenomicsDBGVCFIterator>;

template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, true, true, true>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, true, true, true>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, true, true, false>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, true, true, false>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, true, false, true>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, true, false, true>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, true, false, false>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, true, false, false>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, false, true, true>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, false, true, true>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, false, true, false>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, false, true, false>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, false, false, true>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, false, false, true>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, false, false, false>(VCFWriterFSB& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterFSB, false, false, false>(VCFWriterFSB& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, true, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, true, true, true>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, true, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, true, true, false>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, false, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, true, false, true>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, true, false, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, true, false, false>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, true, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, false, true, true>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, true, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, false, true, false>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, false, true>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, false, false, true>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, false, false, false>(VCFWriterNoOverflow<std::string>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::string>, false, false, false>(VCFWriterNoOverflow<std::string>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, true, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, true, true, true>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, true, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, true, true, false>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, false, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, true, false, true>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, true, false, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, true, false, false>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, true, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, false, true, true>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, true, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, false, true, false>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, false, true>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, false, false, true>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, false, false, false>(VCFWriterNoOverflow<std::ostream>& op) const;
template
bool GTRemapper<GenomicsDBGVCFIterator>::remap_for_row_query_idx<VCFWriterNoOverflow<std::ostream>, false, false, false>(VCFWriterNoOverflow<std::ostream>& op, const size_t row_query_idx) const;
