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

#ifndef VCF_ADAPTER_H
#define VCF_ADAPTER_H

#include "headers.h"
#include "gt_common.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "timer.h"
#include "genomicsdb_config_base.h"
#include "fmt_writer.h"

enum VCFIndexType {
  VCF_INDEX_CSI=0u,
  VCF_INDEX_TBI,
  VCF_INDEX_NONE
};

//Exceptions thrown
class VCFAdapterException : public std::exception {
 public:
  VCFAdapterException(const std::string m="") : msg_("VCFAdapterException : "+m) { ; }
  ~VCFAdapterException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};

//Struct for accessing reference genome
//reference genome access - required for gVCF merge
class ReferenceGenomeInfo {
 public:
  ReferenceGenomeInfo() {
    clear();
    m_reference_faidx = 0;
    m_reference_last_read_pos = -1;
    m_reference_num_bases_read = 0;
    m_reference_last_seq_read = "";
  }
  void clear() {
    m_reference_last_seq_read.clear();
    m_buffer.clear();
  }
  ~ReferenceGenomeInfo() {
    clear();
    if (m_reference_faidx)
      fai_destroy(m_reference_faidx);
  }
  void initialize(const std::string& reference_genome);
  char get_reference_base_at_position(const char* contig, const int64_t pos);
 private:
  int m_reference_last_read_pos;
  int m_reference_num_bases_read;
  std::string m_reference_last_seq_read;
  std::vector<char> m_buffer;
  faidx_t* m_reference_faidx;
};

class VidMapper;
class VCFAdapter {
 public:
  //Returns true if new field added
  static bool add_field_to_hdr_if_missing(bcf_hdr_t* hdr, const VidMapper* id_mapper, const std::string& field_name, int field_type_idx);
 public:
  VCFAdapter(bool open_output=true);
  virtual ~VCFAdapter();
  void clear();
  virtual void initialize(const GenomicsDBConfigBase& config_base);
  //Allocates header
  bcf_hdr_t* initialize_default_header();
  bcf_hdr_t* get_vcf_header() {
    return m_template_vcf_hdr;
  }
  void set_vcf_header(bcf_hdr_t* hdr) {
    m_template_vcf_hdr = hdr;
  }
  /*
   * The line is ready for output
   * Child classes might actually just swap out the pointer so that the actual output is performed by
   * a thread off the critical path
   **/
  virtual void handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size);
  virtual void print_header();
  /*
   * Return true in child class if some output causes buffer to be full. Default: return false
   */
  virtual bool overflow() const {
    return false;
  }
  char get_reference_base_at_position(const char* contig, const int64_t pos) {
    return m_reference_genome_info.get_reference_base_at_position(contig, pos);
  }
  void close_file();
 protected:
  bool m_open_output;
  bcf_hdr_t* m_template_vcf_hdr;
  //Reference genome info
  ReferenceGenomeInfo m_reference_genome_info;
  //Output fptr
  htsFile* m_output_fptr;
  bool m_is_bcf;
  //Index output VCF
  VCFIndexType m_output_VCF_index_type;
  //GenomicsDBConfigBase object
  const GenomicsDBConfigBase* m_config_base_ptr;
#ifdef DO_PROFILING
  //Timer
  Timer m_vcf_serialization_timer;
#endif
};

class BufferedVCFAdapter : public VCFAdapter, public CircularBufferController {
 public:
  BufferedVCFAdapter(unsigned num_circular_buffers, unsigned max_num_entries);
  virtual ~BufferedVCFAdapter();
  void clear();
  virtual void handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size);
  void advance_write_idx();
  void do_output();
  inline bool overflow() const {
    assert(m_config_base_ptr);
    return (m_combined_vcf_records_buffer_sizes[get_write_idx()] >=
            m_config_base_ptr->get_combined_vcf_records_buffer_size_limit());
  }
 private:
  void resize_line_buffer(std::vector<bcf1_t*>& line_buffer, unsigned new_size);
  std::vector<std::vector<bcf1_t*>> m_line_buffers;   //Outer vector for double-buffering
  std::vector<unsigned> m_num_valid_entries;  //One per double-buffer
  std::vector<size_t> m_combined_vcf_records_buffer_sizes;
};

class VCFSerializedBufferAdapter: public VCFAdapter {
 public:
  VCFSerializedBufferAdapter(bool keep_idx_fields_in_bcf_header=true, bool do_output=false)
    : VCFAdapter(false) {
    m_keep_idx_fields_in_bcf_header = keep_idx_fields_in_bcf_header;
    m_rw_buffer = 0;
    //Temporary hts string
    m_hts_string.l = 0u;
    m_hts_string.m = 4096u;
    m_hts_string.s = (char*)malloc(m_hts_string.m);
    assert(m_hts_string.s);
    m_write_fptr = 0;
    m_do_output = do_output;
  }
  ~VCFSerializedBufferAdapter() {
    if (m_hts_string.s && m_hts_string.m > 0)
      free(m_hts_string.s);
    m_hts_string.s = 0;
    m_hts_string.m = 0;
    if (m_write_fptr && m_write_fptr != stdout && m_write_fptr != stderr)
      fclose(m_write_fptr);
    m_write_fptr = 0;
  }
  //Delete copy and move constructors
  VCFSerializedBufferAdapter(const VCFSerializedBufferAdapter& other) = delete;
  VCFSerializedBufferAdapter& operator=(const VCFSerializedBufferAdapter& other) = delete;
  VCFSerializedBufferAdapter(VCFSerializedBufferAdapter&& other) = delete;
  void initialize(const GenomicsDBConfigBase& config_base);
  void set_buffer(RWBuffer& buffer) {
    m_rw_buffer = &buffer;
  }
  void print_header();
  void handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size);
  inline bool overflow() const {
    assert(m_rw_buffer);
    return (m_rw_buffer->m_num_valid_bytes >=
            m_config_base_ptr->get_combined_vcf_records_buffer_size_limit());
  }
  void do_output() {
    assert(m_write_fptr);
    assert(m_rw_buffer);
    auto write_size = fwrite(&(m_rw_buffer->m_buffer[0]), 1u,  m_rw_buffer->m_num_valid_bytes, m_write_fptr);
    assert(write_size == m_rw_buffer->m_num_valid_bytes);
  }
 private:
  bool m_keep_idx_fields_in_bcf_header;
  RWBuffer* m_rw_buffer;
  kstring_t m_hts_string;
  bool m_do_output;
  FILE* m_write_fptr;
};

#endif
