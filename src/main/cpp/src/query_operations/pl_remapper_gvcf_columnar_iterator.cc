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

#include "pl_remapper_template_definition.h"
#include "vcf_fmt_writer.h"
#include "genomicsdb_iterators.h"

//Explicit instantiation
template class PLRemapper<GenomicsDBGVCFIterator>;

template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, int, true, true>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, int, true, false>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, int, false, true>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, int, false, false>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, float, true, true>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, float, true, false>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, float, false, true>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterFSB, float, false, false>(VCFWriterFSB& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, int, true, true>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, int, true, false>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, int, false, true>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, int, false, false>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, float, true, true>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, float, true, false>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, float, false, true>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::string>, float, false, false>(VCFWriterNoOverflow<std::string>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, int, true, true>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, int, true, false>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, int, false, true>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, int, false, false>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, float, true, true>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, float, true, false>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, float, false, true>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);
template
bool PLRemapper<GenomicsDBGVCFIterator>::remap_all_queried_valid_rows<VCFWriterNoOverflow<std::ostream>, float, false, false>(VCFWriterNoOverflow<std::ostream>& op, unsigned PL_query_idx);

