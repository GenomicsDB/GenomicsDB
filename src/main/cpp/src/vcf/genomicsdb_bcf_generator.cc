/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#include "genomicsdb_bcf_generator.h"

unsigned GenomicsDBBCFGenerator_NUM_ENTRIES_IN_CIRCULAR_BUFFER=1u;
GenomicsDBBCFGenerator::GenomicsDBBCFGenerator(const std::string& loader_config_file, 
    const genomicsdb_pb::ExportConfiguration *query_config_pb,
    const char* chr, const long start, const long end,
    int my_rank, size_t buffer_capacity, size_t tiledb_segment_size, const char* output_format,
    const bool produce_header_only,
    const bool use_missing_values_only_not_vector_end, const bool keep_idx_fields_in_bcf_header)
  : m_buffer_control(GenomicsDBBCFGenerator_NUM_ENTRIES_IN_CIRCULAR_BUFFER),
    m_vcf_adapter(keep_idx_fields_in_bcf_header),
    m_produce_header_only(produce_header_only)
#ifdef DO_PROFILING
  , m_timer()
#endif
{
  m_done = false;
  //Buffer sizing
  m_buffers.resize(GenomicsDBBCFGenerator_NUM_ENTRIES_IN_CIRCULAR_BUFFER, RWBuffer(buffer_capacity+32768u)); //pad buffer to minimize reallocations
  GenomicsDBImportConfig loader_config;
  //Parse loader JSON file
  //If the loader JSON is not specified, vid_mapping_file and callset_mapping_file must be specified in the query JSON
  if (!(loader_config_file.empty())) {
    loader_config.read_from_file(loader_config_file, my_rank);
    m_query_config.update_from_loader(loader_config, my_rank);
  }
  //Parse query protobuf - ensures that vid, callset may be obtained from loader (if exists)
  m_query_config.read_from_PB(query_config_pb, my_rank);
  m_query_config.validate();
  if (!(loader_config_file.empty()))
    m_query_config.subset_query_column_ranges_based_on_partition(loader_config, my_rank);
  m_query_config.set_vcf_output_format(output_format);
  m_vcf_adapter.initialize(m_query_config);
  auto& vid_mapper = m_query_config.get_vid_mapper();
  //Specified chromosome and start end
  if (chr && strlen(chr) > 0u) {
    ContigInfo contig_info;
    auto found_contig = vid_mapper.get_contig_info(chr, contig_info);
    if (!found_contig)
      throw GenomicsDBJNIException(std::string("Could not find TileDB column interval for contig: ")+chr);
    int64_t column_begin = contig_info.m_tiledb_column_offset + static_cast<int64_t>(start) - 1; //since VCF positions are 1 based
    // check that end position doesn't go beyond contig length
    int64_t end_check = std::min<int64_t>(end, contig_info.m_length);
    int64_t column_end = contig_info.m_tiledb_column_offset + static_cast<int64_t>(end_check) - 1; //since VCF positions are 1 based
    m_query_config.set_column_interval_to_query(column_begin, column_end);
  }
  else if (start != 0 || end != 0) {
    m_query_config.set_column_interval_to_query(start, end);
  }
  m_storage_manager = new VariantStorageManager(m_query_config.get_workspace(my_rank),
                                                m_query_config.get_segment_size(),
                                                m_query_config.enable_shared_posixfs_optimizations());
  m_query_processor = new VariantQueryProcessor(m_storage_manager,
      m_query_config.get_array_name(my_rank),
      vid_mapper);
  m_query_processor->do_query_bookkeeping(m_query_processor->get_array_schema(), m_query_config, vid_mapper, true);
  //Must set buffer before constructing BroadCombinedGVCFOperator
  set_write_buffer();
  m_combined_bcf_operator = new BroadCombinedGVCFOperator(m_vcf_adapter, vid_mapper, m_query_config,
      use_missing_values_only_not_vector_end);
  m_query_column_interval_idx = 0u;
  if (produce_header_only)
    m_scan_state.set_done(true);
  else {
    m_query_processor->scan_and_operate(m_query_processor->get_array_descriptor(), m_query_config, *m_combined_bcf_operator, m_query_column_interval_idx,
                                        true, &m_scan_state);
  }
#ifdef DO_PROFILING
  m_timer.stop();
#endif
}

GenomicsDBBCFGenerator::~GenomicsDBBCFGenerator() {
  //Delete iterator before storage manager is deleted
  if (m_scan_state.get_iterator()) {
    delete m_scan_state.get_iterator();
    m_scan_state.set_iterator(0);
  }
  m_buffers.clear();
  if (m_combined_bcf_operator)
    delete m_combined_bcf_operator;
  m_combined_bcf_operator = 0;
  if (m_query_processor)
    delete m_query_processor;
  m_query_processor = 0;
  if (m_storage_manager)
    delete m_storage_manager;
  m_storage_manager = 0;
#ifdef DO_PROFILING
  m_timer.print("GenomicsDBBCFGenerator", std::cerr);
#endif
}

void GenomicsDBBCFGenerator::produce_next_batch() {
  if (m_done)
    return;
  auto num_bytes_produced = 0ull;
  while (num_bytes_produced == 0u) {
    if (m_scan_state.end()) {
      ++m_query_column_interval_idx;
      if (m_produce_header_only ||
          m_query_column_interval_idx >= m_query_config.get_num_column_intervals()) {
        reset_read_buffer();
        m_done = true;
        return;
      }
      m_scan_state.reset();
    }
    reset_read_buffer();
    set_write_buffer();
    m_query_processor->scan_and_operate(m_query_processor->get_array_descriptor(), m_query_config, *m_combined_bcf_operator, m_query_column_interval_idx,
                                        true, &m_scan_state);
    num_bytes_produced = m_buffers[m_buffer_control.get_read_idx()].m_num_valid_bytes;
  }
}

size_t GenomicsDBBCFGenerator::read_and_advance(uint8_t* dst, size_t offset, size_t n) {
#ifdef DO_PROFILING
  m_timer.start();
#endif
  auto total_bytes_advanced = 0ull;
  if (n == SIZE_MAX)
    produce_next_batch();
  else
    while (total_bytes_advanced < n && m_buffer_control.get_num_entries_with_valid_data() > 0u) {
      auto& curr_buffer = m_buffers[m_buffer_control.get_read_idx()];
      //Minimum of space left in buffer and #bytes to fetch
      auto num_bytes_in_curr_buffer = std::min<size_t>(n-total_bytes_advanced, curr_buffer.get_num_remaining_bytes());
      if (dst)
        memcpy_s(dst+offset+total_bytes_advanced, num_bytes_in_curr_buffer,
                 &(curr_buffer.m_buffer[curr_buffer.m_next_read_idx]), num_bytes_in_curr_buffer);
      //Advance marker
      curr_buffer.m_next_read_idx += num_bytes_in_curr_buffer;
      total_bytes_advanced += num_bytes_in_curr_buffer;
      //Reached end of buffer
      if (curr_buffer.m_next_read_idx >= curr_buffer.m_num_valid_bytes)
        produce_next_batch();
    }
#ifdef DO_PROFILING
  m_timer.stop();
#endif
  return total_bytes_advanced;
}

uint8_t GenomicsDBBCFGenerator::read_next_byte() {
  uint8_t tmp;
  auto num_bytes_read = read_and_advance(&tmp, 0u, 1u);
  return (num_bytes_read > 0u) ? tmp : -1;
}

void GenomicsDBBCFGenerator::set_write_buffer() {
  m_vcf_adapter.set_buffer(m_buffers[m_buffer_control.get_write_idx()]);
  m_buffer_control.advance_write_idx();
}

void GenomicsDBBCFGenerator::reset_read_buffer() {
  auto& curr_buffer = m_buffers[m_buffer_control.get_read_idx()];
  curr_buffer.m_next_read_idx = 0ull;
  curr_buffer.m_num_valid_bytes = 0ull;
  m_buffer_control.advance_read_idx();
}
