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

#ifndef VCF2BINARY_H
#define VCF2BINARY_H

#ifdef HTSDIR

#include "headers.h"
#include "vid_mapper.h"
#include "htslib/synced_bcf_reader.h"
#include "gt_common.h"
#include "histogram.h"
#include "tiledb_loader_file_base.h"

//Exceptions thrown
class VCF2BinaryException : public std::exception {
 public:
  VCF2BinaryException(const std::string m="") : msg_("VCF2BinaryException : "+m) { ; }
  ~VCF2BinaryException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};

//Base class for VCFReader and VCFBufferReader
//Contains header, line and a buffer for fetching data
class VCFReaderBase : public virtual GenomicsDBImportReaderBase {
 public:
  VCFReaderBase(const bool is_file_reader)
    : GenomicsDBImportReaderBase(is_file_reader) {
    m_hdr = 0;
    m_line = bcf_init();
  }
  ~VCFReaderBase() {
    if (m_hdr)
      bcf_hdr_destroy(m_hdr);
    m_hdr = 0;
    if (m_line)
      bcf_destroy(m_line);
    m_line = 0;
  }
  //Delete copy and move constructors
  VCFReaderBase(const VCFReaderBase& other) = delete;
  VCFReaderBase& operator=(const VCFReaderBase& other) = delete;
  VCFReaderBase(VCFReaderBase&& other) = delete;
  virtual void initialize(const char* filename,
                          const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, const bool open_file);
  bcf_hdr_t* get_header() {
    return m_hdr;
  }
  bcf1_t* get_line() {
    return m_is_record_valid ? m_line : 0;
  }
 protected:
  bcf_hdr_t* m_hdr;
  bcf1_t* m_line;
};

//Data is read from the buffer that is filled by an external agent
class VCFBufferReader : public BufferReaderBase, public VCFReaderBase {
 public:
  VCFBufferReader(const size_t buffer_capacity, const bool is_bcf,
                  const uint8_t* init_buffer, const size_t init_num_valid_bytes)
    : GenomicsDBImportReaderBase(false), BufferReaderBase(buffer_capacity), VCFReaderBase(false), m_is_bcf(is_bcf) {
    m_buffer.resize(std::max(buffer_capacity, init_num_valid_bytes));
    memcpy_s(&(m_buffer[0]), init_num_valid_bytes, init_buffer, init_num_valid_bytes);
    m_num_valid_bytes_in_buffer = init_num_valid_bytes;
  }
  //Delete copy and move constructors
  VCFBufferReader(const VCFBufferReader& other) = delete;
  VCFBufferReader& operator=(const VCFBufferReader& other) = delete;
  VCFBufferReader(VCFBufferReader&& other) = delete;
  virtual ~VCFBufferReader() = default;
  void initialize(const char* stream_name,
                  const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, const bool open_file);
  void read_and_advance();
 private:
  bool m_is_bcf;
};

//Wrapper around VCF's file I/O functions
//Capability of using index only during seek to minimize memory consumption
class VCFReader : public FileReaderBase, public VCFReaderBase {
 public:
  VCFReader();
  //Delete move and copy constructors
  VCFReader(const VCFReader& other) = delete;
  VCFReader& operator=(const VCFReader& other) = delete;
  VCFReader(VCFReader&& other) = delete;
  //Destructor
  virtual ~VCFReader();
  void initialize(const char* filename,
                  const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, bool open_file);
  //Abstract virtual functions from base class that must be over-ridden
  void add_reader();
  void remove_reader();
  void read_and_advance();
  //Helper functions
  void seek_read_advance(const char* contig, const int pos, bool discard_index);
 private:
  bcf_srs_t* m_indexed_reader;
  htsFile* m_fptr;
  kstring_t m_vcf_file_buffer;
};

class VCFColumnPartition : public File2TileDBBinaryColumnPartitionBase {
  friend class VCF2Binary;
 public:
  class VCFGetBufferWrapper {
    //buffers for vcf get functions - 4 KB
   public:
    VCFGetBufferWrapper(const size_t num_bytes=4096u) {
      m_buffer = (uint8_t*)malloc(num_bytes);
      if (m_buffer == 0)
        throw VCF2BinaryException("Malloc failure");
      m_capacity = num_bytes;
      m_num_values = 0;
    }
    //Delete copy constructor
    VCFGetBufferWrapper(const VCFGetBufferWrapper& other) = delete;
    VCFGetBufferWrapper& operator=(const VCFGetBufferWrapper& other) = delete;
    //Define move constructor explicitly
    VCFGetBufferWrapper(VCFGetBufferWrapper&& other);
    //Destructor
    ~VCFGetBufferWrapper() {
      if (m_buffer && m_capacity)
        free(m_buffer);
      m_buffer = 0;
      m_capacity = 0u;
      m_num_values = 0;
    }
    uint8_t* m_buffer;
    size_t m_capacity;
    int64_t m_num_values;
  };
  /*
   * Primary constructor
   */
  VCFColumnPartition(const bool prefetch_fields,
                     const size_t num_INFO_fields, const size_t num_FORMAT_fields)
    : File2TileDBBinaryColumnPartitionBase() {
    //Initialize as invalid
    m_local_contig_idx = -1;
    m_contig_position = -1;
    m_contig_tiledb_column_offset = -1;
    //If prefetch is enabled, allocate buffer per field, else single buffer
    //Add 1 buffer for END field
    m_vcf_get_buffer_vec.emplace_back(prefetch_fields ? num_INFO_fields+1u : 1u);
    //if prefetch is enabled, allocate buffer per field
    if (prefetch_fields)
      m_vcf_get_buffer_vec.emplace_back(num_FORMAT_fields);
    m_split_output_fptr = 0;
  }
  //Delete copy constructor
  VCFColumnPartition(const VCFColumnPartition& other) = delete;
  VCFColumnPartition& operator=(const VCFColumnPartition& other) = delete;
  //Define move constructor
  VCFColumnPartition(VCFColumnPartition&& other);
  ~VCFColumnPartition();
  int64_t get_column_position_in_record() const {
    auto vcf_reader_ptr = dynamic_cast<VCFReader*>(m_base_reader_ptr);
    assert(vcf_reader_ptr);
    auto line = vcf_reader_ptr->get_line();
    assert(line);
    return (m_contig_tiledb_column_offset + static_cast<int64_t>(line->pos));
  }
  bcf_hdr_t* get_header() {
    auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(m_base_reader_ptr);
    assert(vcf_reader_ptr);
    return vcf_reader_ptr->get_header();
  }
  VCFGetBufferWrapper& get_vcf_get_buffer_wrapper(const bool prefetch_fields,
      const bool is_INFO_field, const bool is_END_field,
      const unsigned idx_in_vcf_fields_vector);
 protected:
  //Position in contig from which to fetch next batch of cells
  int m_local_contig_idx;
  int64_t m_contig_position;  //position in contig (0-based)
  int64_t m_contig_tiledb_column_offset;
  //Buffers for obtaining data from htslib
  //Outer vector of size 2 - INFO, FORMAT - if prefetch enabled, else size 1
  //Inner vector depends on #INFO and FORMAT fields imported if prefetch enabled, else size 1
  std::vector<std::vector<VCFGetBufferWrapper> > m_vcf_get_buffer_vec;
  //File pointer to output partition data - useful when splitting files
  htsFile* m_split_output_fptr;
};

class VCF2Binary : public File2TileDBBinaryBase {
 public:
  VCF2Binary(const std::string& vcf_filename, const std::vector<std::vector<std::string>>& vcf_fields,
             unsigned file_idx, VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
             size_t max_size_per_callset,
             bool treat_deletions_as_intervals,
             bool parallel_partitions=false, bool noupdates=true, bool close_file=false, bool discard_index=false);
  VCF2Binary(const std::string& stream_name, const std::vector<std::vector<std::string>>& vcf_fields,
             unsigned file_idx, const int64_t buffer_stream_idx,
             VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
             const size_t vcf_buffer_reader_buffer_size, const bool vcf_buffer_reader_is_bcf,
             const uint8_t* vcf_buffer_reader_init_buffer, const size_t vcf_buffer_reader_init_num_valid_bytes,
             size_t max_size_per_callset,
             bool treat_deletions_as_intervals);
  //Delete default copy constructor as it is incorrect
  VCF2Binary(const VCF2Binary& other) = delete;
  VCF2Binary& operator=(const VCF2Binary& other) = delete;
  //Define move constructor explicitly
  VCF2Binary(VCF2Binary&& other);
  ~VCF2Binary();
  void clear();
  //Initialization functions
  void initialize(const std::vector<ColumnRange>& partition_bounds);
  //Abstract virtual functions in base class that must be defined
  /*
   * Set order of enabled callsets
   * In VCFs, each line contains data from all callsets
   */
  void set_order_of_enabled_callsets(int64_t& order_value, std::vector<int64_t>& tiledb_row_idx_to_order) const;
  /*
   * List active row idxs
   * In VCFs, each line contains data from all callsets
   */
  void list_active_row_idxs(const ColumnPartitionBatch& partition_batch, int64_t& row_idx_offset, std::vector<int64_t>& row_idx_vec) const;
  void initialize_column_partitions(const std::vector<ColumnRange>& partition_bounds);
  /*
   * Create the subclass of File2TileDBBinaryColumnPartitionBase that must be used
   */
  File2TileDBBinaryColumnPartitionBase* create_new_column_partition_object() const {
    return dynamic_cast<File2TileDBBinaryColumnPartitionBase*>(new VCFColumnPartition(
             m_prefetch_all_VCF_fields_in_record, (*m_vcf_fields)[BCF_HL_INFO].size(),
             (*m_vcf_fields)[BCF_HL_FMT].size()));
  }
  /*
   * Create the subclass of GenomicsDBImportReaderBase that must be used
   */
  GenomicsDBImportReaderBase* create_new_reader_object(const std::string& filename, bool open_file) const;
  bool convert_record_to_binary(std::vector<uint8_t>& buffer, File2TileDBBinaryColumnPartitionBase& partition_info);
  bool seek_and_fetch_position(File2TileDBBinaryColumnPartitionBase& partition_info, bool& is_read_buffer_empty, bool force_seek, bool advance_reader);
  uint64_t get_num_callsets_in_record(const File2TileDBBinaryColumnPartitionBase& partition_info) const {
    return m_enabled_local_callset_idx_vec.size();
  }
  //Helper functions
  void update_local_contig_idx(VCFColumnPartition& vcf_partition, const bcf1_t* line);
  //VCF->TileDB conversion functions
  bool convert_VCF_to_binary_for_callset(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
                                         size_t size_per_callset, uint64_t enabled_callsets_idx);
  /*
   * field_type_idx: BCF_HL_*
   */
  template<class FieldType>
  bool convert_field_to_tiledb(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
                               int64_t& buffer_offset, const int64_t buffer_offset_limit, int local_callset_idx,
                               const std::string& field_name, unsigned field_type_idx, const unsigned idx_in_vcf_fields_vector);
  template<typename FieldType>
  void fetch_field_from_vcf_record(VCFColumnPartition::VCFGetBufferWrapper& vcf_get_buffer_wrapper,
                                   const bcf_hdr_t* hdr, bcf1_t* line,
                                   const std::string& field_name, const int field_type_idx, const int bcf_ht_type);
  //Print partitions of the file - useful when splitting files into partitions
  /*
   * Opens the file for partition - useful when printing data for a specific partition (splitting files)
   */
  bool open_partition_output_file(const std::string& results_directory, std::string& output_filename,
                                  const std::string& output_type, File2TileDBBinaryColumnPartitionBase& partition_info, const unsigned partition_idx);
  /*
   * Prints data of the partition
   */
  void write_partition_data(File2TileDBBinaryColumnPartitionBase& partition_info);
  /*
   * Closes the file for partition - useful when printing data for a specific partition (splitting files)
   */
  void close_partition_output_file(File2TileDBBinaryColumnPartitionBase& partition_info);
  void set_discard_missing_GTs_flag(const bool value) {
    m_discard_missing_GTs = value;
  }
 private:
  bool m_discard_index;
  bool m_import_ID_field;
  bool m_discard_missing_GTs;
  bool m_discard_current_record;
  bool m_prefetch_all_VCF_fields_in_record;
  //Vector of vector of strings, outer vector has 2 elements - 0 for INFO, 1 for FORMAT
  const std::vector<std::vector<std::string>>* m_vcf_fields;
  //Local contig idx to global contig idx
  std::vector<int> m_local_contig_idx_to_global_contig_idx;
  //Local field idx to global field idx
  std::vector<int> m_local_field_idx_to_global_field_idx;
  //For VCFBufferReader
  size_t m_vcf_buffer_reader_buffer_size;
  bool m_vcf_buffer_reader_is_bcf;
  size_t m_vcf_buffer_reader_init_num_valid_bytes;
  //The buffer in with the serialized VCF header used for initialization
  const uint8_t* m_vcf_buffer_reader_init_buffer;
};

#endif //ifdef HTSDIR

#endif
