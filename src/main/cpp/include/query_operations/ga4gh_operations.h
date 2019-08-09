/**
 * @file ga4gh_operations.h
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Omics Data Automation, Inc.
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
 *
 * @section DESCRIPTION
 *
 * Fill up GA4GH Variant and Call structures from GenomicsDB 
 *
 **/

#ifndef GA4GH_OPERATIONS_H
#define GA4GH_OPERATIONS_H

#include "variant_operations.h"
#include "variants.pb.h"

class GA4GHCallCreator : public SingleCellOperatorBase {
  public:
    GA4GHCallCreator(const VidMapper& vid_mapper, std::ostream& stream, const size_t buffer_size);
    GA4GHCallCreator(const VidMapper& vid_mapper, std::ostream& stream, const bool unlimited);
    void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
	const VariantArraySchema& schema);
    void fill_ga4gh_call(const VariantQueryConfig& query_config, const GenomicsDBColumnarCell& cell);
    void finalize();
  private:
    void common_initialization(const VidMapper& vid_mapper, std::ostream& stream,
	const size_t buffer_size, const bool unlimited);
    //Members
    std::ostream* m_stream;
    size_t m_stream_buffer_size;
    bool m_limited_stream_buffer;
    const VidMapper* m_vid_mapper;
    ga4gh::Variant m_variant;
    ga4gh::Call m_call;
    ga4gh::AttributeValueList m_attr_list;
  private:
    //For performance reasons
    std::string m_contig;
    int64_t m_contig_begin_tiledb_column_offset;
    int64_t m_contig_end_tiledb_column_offset;
    unsigned m_ploidy_step_value;
    bool m_first_call;
};

#endif
