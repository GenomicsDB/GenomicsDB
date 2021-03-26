/**
 * The MIT License (MIT)
 * Copyright (c) 2018 Intel Corporation
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

#ifndef GENOMICSDB_CONFIG_BASE_H
#define GENOMICSDB_CONFIG_BASE_H

#include "vid_mapper.h"
#include <google/protobuf/message.h>

namespace genomicsdb_pb {
  class ExportConfiguration;
}

// Overridding env variables
#define ENABLE_SHARED_POSIXFS_OPTIMIZATIONS "GENOMICSDB_SHARED_POSIXFS_OPTIMIZATIONS"
#define DISABLE_DELTA_ENCODE_OFFSETS "GENOMICSDB_OFFSETS_DISABLE_DELTA_ENCODE"
#define DISABLE_DELTA_ENCODE_COORDS "GENOMICSDB_COORDS_DISABLE_DELTA_ENCODE"
#define ENABLE_BIT_SHUFFLE_GT "GENOMICSDB_GT_ENABLE_BIT_SHUFFLE"
#define ENABLE_LZ4_GT "GENOMICSDB_GT_ENABLE_LZ4_COMPRESSION"

//Exceptions thrown
class GenomicsDBConfigException : public std::exception {
 public:
  GenomicsDBConfigException(const std::string m="") : msg_("GenomicsDBConfigException : "+m) { ; }
  ~GenomicsDBConfigException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};

class GenomicsDBImportConfig;

class GenomicsDBConfigBase {
 public:
  GenomicsDBConfigBase();
  const std::string& get_workspace(const int rank) const;
  void set_workspace(const std::string& workspace);
  const std::string& get_array_name(const int rank) const;
  void set_array_name(const std::string& array_name);
  ColumnRange get_column_partition(const int rank, const unsigned idx=0u) const;
  TileDBRowRange get_row_partition(const int rank, const unsigned idx=0u) const;
  const std::vector<ColumnRange> get_sorted_column_partitions() const {
    return m_sorted_column_partitions;
  }
  const std::vector<ColumnRange>& get_query_column_ranges(const int rank) const;
  void set_query_column_ranges(const std::vector<ColumnRange>& column_ranges);
  const std::vector<TileDBRowRange>& get_query_row_ranges(const int rank) const;
  void set_query_row_ranges(const std::vector<TileDBRowRange>& row_ranges);
  inline const std::string& get_query_filter() const {
    return m_query_filter;
  }
  inline size_t get_segment_size() const {
    return m_segment_size;
  }
  void set_segment_size(const size_t v) {
    m_segment_size = v;
  }
  inline unsigned get_determine_sites_with_max_alleles() const {
    return m_determine_sites_with_max_alleles;
  }
  inline unsigned get_max_diploid_alt_alleles_that_can_be_genotyped() const {
    return m_max_diploid_alt_alleles_that_can_be_genotyped;
  }
  inline unsigned get_max_genotype_count() const {
    return m_max_genotype_count;
  }
  void set_combined_vcf_records_buffer_size_limit(const size_t val) {
    m_combined_vcf_records_buffer_size_limit = val;
  }
  inline size_t get_combined_vcf_records_buffer_size_limit() const {
    return m_combined_vcf_records_buffer_size_limit;
  }
  void set_vcf_header_filename(const std::string& vcf_header_filename);
  const std::string& get_vcf_header_filename() const {
    return m_vcf_header_filename;
  }
  const std::string& get_vcf_output_filename() const {
    return m_vcf_output_filename;
  }
  void set_vcf_output_filename(const std::string& output_filename) {
     m_vcf_output_filename = output_filename;
  }
  void set_vcf_output_format(const std::string& output_format);
  const std::string& get_vcf_output_format() const {
    return m_vcf_output_format;
  }
  const std::string& get_reference_genome() const {
    return m_reference_genome;
  }
  void set_reference_genome(const std::string& reference_genome) {
    m_reference_genome = reference_genome;
  }
  const bool produce_GT_field() const {
    return m_produce_GT_field;
  }
  const bool produce_FILTER_field() const {
    return m_produce_FILTER_field;
  }
  const bool sites_only_query() const {
    return m_sites_only_query;
  }
  const bool index_output_VCF() const {
    return m_index_output_VCF;
  }
  void set_index_output_VCF(bool index_output_vcf) {
    m_index_output_VCF = index_output_vcf;
  }
  const bool produce_GT_with_min_PL_value_for_spanning_deletions() const {
    return m_produce_GT_with_min_PL_value_for_spanning_deletions;
  }
  const VidMapper& get_vid_mapper() const {
    return m_vid_mapper;
  }
  //Utility functions
  static ColumnRange verify_contig_position_and_get_tiledb_column_interval(const ContigInfo& contig_info,
      const int64_t begin, const int64_t end);
  const std::string& get_callset_mapping_file() const {
    return m_callset_mapping_file;
  }
  void set_callset_mapping_file(const std::string& callset_mapping_file) {
    m_callset_mapping_file = callset_mapping_file;
  }
  const std::string& get_vid_mapping_file() const {
    return m_vid_mapping_file;
  }
  void set_vid_mapping_file(const std::string& vid_mapping_file) {
    m_vid_mapping_file = vid_mapping_file;
  }
  //Sometimes information is present in the loader - copy over
  void update_from_loader(const GenomicsDBImportConfig& loader_config, const int rank=0);
  void subset_query_column_ranges_based_on_partition(const GenomicsDBImportConfig& loader_config, const int rank=0);
  inline TileDBRowRange get_row_bounds() const {
    return TileDBRowRange(m_lb_callset_row_idx, m_ub_callset_row_idx);
  }
  inline uint64_t get_num_rows_within_bounds() const {
    return m_ub_callset_row_idx - m_lb_callset_row_idx + 1ull;
  }
  inline bool enable_shared_posixfs_optimizations() const {
    return is_set_with_env_override(m_enable_shared_posixfs_optimizations, ENABLE_SHARED_POSIXFS_OPTIMIZATIONS);
  }
  inline bool produce_lowercase_alleles_in_soft_masked_regions() const {
    return m_produce_lowercase_alleles_in_soft_masked_regions;
  }
  void scan_whole_array();
  const std::vector<std::string>& get_attributes() const { return m_attributes; }
  //JSON parsing functions
  rapidjson::Document read_from_file(const std::string& filename, const int rank=0);
  rapidjson::Document read_from_JSON_string(const std::string& str, const int rank=0);
  void read_from_JSON(const rapidjson::Document& json_doc, const int rank=0);
  void read_and_initialize_vid_and_callset_mapping_if_available(const rapidjson::Document& json_doc, const int rank);
  //Protobuf parsing functions
  void read_from_PB_binary_string(const std::string& str, const int rank=0);
  void read_from_PB(const genomicsdb_pb::ExportConfiguration* x, const int rank=0);
  void read_and_initialize_vid_and_callset_mapping_if_available(const genomicsdb_pb::ExportConfiguration*);
  static void get_pb_from_json_file(google::protobuf::Message* pb_config, const std::string& json_file);
 protected:
  bool m_single_workspace_path;
  bool m_single_array_name;
  bool m_single_query_column_ranges_vector;
  bool m_column_partitions_specified;
  bool m_single_query_row_ranges_vector;
  bool m_row_partitions_specified;
  bool m_scan_whole_array;
  //GATK CombineGVCF does not produce GT field by default - option to produce GT
  bool m_produce_GT_field;
  //GATK CombineGVCF does not produce FILTER field by default - option to produce FILTER
  bool m_produce_FILTER_field;
  //index output VCF file
  bool m_index_output_VCF;
  //sites-only query - doesn't produce any of the FORMAT fields
  bool m_sites_only_query;
  //when producing GT, use the min PL value GT for spanning deletions
  bool m_produce_GT_with_min_PL_value_for_spanning_deletions;
  std::vector<std::string> m_workspaces;
  std::vector<std::string> m_array_names;
  std::vector<std::vector<ColumnRange>> m_column_ranges;
  std::vector<std::vector<TileDBRowRange>> m_row_ranges;
  std::vector<std::string> m_attributes;
  std::vector<ColumnRange> m_sorted_column_partitions;
  std::vector<TileDBRowRange> m_sorted_row_partitions;
  //Lower and upper bounds of callset row idx to import in this invocation
  int64_t m_lb_callset_row_idx;
  int64_t m_ub_callset_row_idx;
  //Filter Expression for reading arrays
  std::string m_query_filter;
  //TileDB segment size
  size_t m_segment_size;
  //VCF output parameters
  std::string m_vcf_header_filename;
  std::string m_reference_genome;
  std::string m_vcf_output_filename;
  std::string m_vcf_output_format;
  //Count max #alt alleles , don't create combined gVCF
  unsigned m_determine_sites_with_max_alleles;
  //Max diploid alleles for which fields whose length is equal to the number of genotypes can be produced (such as PL)
  unsigned m_max_diploid_alt_alleles_that_can_be_genotyped;
  //Max #genotypes for which fields whose length is equal to the number of genotypes can be produced (such as PL)
  unsigned m_max_genotype_count;
  //Buffer size for combined vcf records
  size_t m_combined_vcf_records_buffer_size_limit;
  //VidMapper
  VidMapper m_vid_mapper;
  //Might be empty strings if using Protobuf
  std::string m_vid_mapping_file;
  std::string m_callset_mapping_file;
  //Enable optimizations (disable file locking and enable keep file handles open until finalization
  bool m_enable_shared_posixfs_optimizations;
  bool m_produce_lowercase_alleles_in_soft_masked_regions;
 public:
  //Static convenience member
  static std::unordered_map<std::string, bool> m_vcf_output_format_to_is_bcf_flag;
  static bool is_set_with_env_override(const bool field, const std::string& env) {
    auto env_var = getenv(env.c_str());
    if (env_var) {
      if (strcasecmp(env_var, "true") == 0 || strcmp(env_var, "1") == 0) {
        return true;
      } else {
        return false;
      }
    } else {
      return field;
    }
  }
};

class GenomicsDBImportConfig : public GenomicsDBConfigBase {
 public:
  GenomicsDBImportConfig();
  void read_from_file(const std::string& filename, int rank=0);
  inline bool is_partitioned_by_row() const {
    return m_row_based_partitioning;
  }
  inline bool is_partitioned_by_column() const {
    return !m_row_based_partitioning;
  }
  inline int64_t get_max_num_rows_in_array() const {
    return m_max_num_rows_in_array;
  }
  inline bool offload_vcf_output_processing() const {
    return m_offload_vcf_output_processing;
  }
  inline bool ignore_cells_not_in_partition() const {
    return m_ignore_cells_not_in_partition;
  }
  inline bool compress_tiledb_array() const {
    return m_compress_tiledb_array;
  }
  inline bool disable_synced_writes() const {
    return m_disable_synced_writes;
  }
  inline bool delete_and_create_tiledb_array() const {
    return m_delete_and_create_tiledb_array;
  }
  inline size_t get_segment_size() const {
    return m_segment_size;
  }
  inline size_t get_num_cells_per_tile() const {
    return m_num_cells_per_tile;
  }
  inline int64_t get_tiledb_compression_level() const {
    return m_tiledb_compression_level;
  }
  inline bool disable_delta_encode_offsets() const {
    return is_set_with_env_override(m_disable_delta_encode_offsets, DISABLE_DELTA_ENCODE_OFFSETS);
  }
  inline bool disable_delta_encode_coords() const {
    return is_set_with_env_override(m_disable_delta_encode_coords, DISABLE_DELTA_ENCODE_COORDS);
  }
  inline bool enable_bit_shuffle_gt() const {
    return is_set_with_env_override(m_enable_bit_shuffle_gt, ENABLE_BIT_SHUFFLE_GT);
  }
  inline bool enable_lz4_compression_gt() const {
    return is_set_with_env_override(m_enable_lz4_compression_gt, ENABLE_LZ4_GT);
  }
  inline bool fail_if_updating() const {
    return m_fail_if_updating;
  }
  inline bool consolidate_tiledb_array_after_load() const {
    return m_consolidate_tiledb_array_after_load;
  }
  inline bool discard_missing_GTs() const {
    return m_discard_missing_GTs;
  }
  inline bool no_mandatory_VCF_fields() const {
    return m_no_mandatory_VCF_fields;
  }
  inline bool treat_deletions_as_intervals() const {
    return m_treat_deletions_as_intervals;
  }
  inline bool produce_tiledb_array() const {
    return m_produce_tiledb_array;
  }
  inline bool produce_combined_vcf() const {
    return m_produce_combined_vcf;
  }
  inline bool discard_vcf_index() const {
    return  m_discard_vcf_index;
  }
  inline int get_num_parallel_vcf_files() const {
    return m_num_parallel_vcf_files;
  }
  inline bool is_row_based_partitioning() const {
    return m_row_based_partitioning;
  }

 protected:
  bool m_standalone_converter_process;
  bool m_treat_deletions_as_intervals;
  bool m_produce_combined_vcf;
  bool m_produce_tiledb_array;
  bool m_compress_tiledb_array;
  bool m_disable_synced_writes;
  bool m_delete_and_create_tiledb_array;
  bool m_row_based_partitioning;
  //do ping-pong buffering
  bool m_do_ping_pong_buffering;
  //Offload VCF output processing to another thread
  bool m_offload_vcf_output_processing;
  //Ignore cells that do not belong to this partition
  bool m_ignore_cells_not_in_partition;
  //Flag that controls whether the VCF indexes should be discarded to reduce memory consumption
  bool m_discard_vcf_index;
  unsigned m_num_entries_in_circular_buffer;
  //#VCF files to open/process in parallel
  int m_num_parallel_vcf_files;
  int m_num_converter_processes;
  int64_t m_per_partition_size;
  int64_t m_max_size_per_callset;
  //max #rows - defining domain of the array
  int64_t m_max_num_rows_in_array;
  //segment size for TileDB array
  size_t m_segment_size;
  //TileDB array #cells/tile
  size_t m_num_cells_per_tile;
  //TileDB compression level
  int m_tiledb_compression_level;
  //flag to disallow TileDB pre compression filter Delta Encoding for offsets to fields
  bool m_disable_delta_encode_offsets;
  //flag to disallow TileDB pre compression filter Delta Encoding for coordinates
  bool m_disable_delta_encode_coords;
  //flag to allow TileDB pre compression filter Bit Shuffle for GT fields
  bool m_enable_bit_shuffle_gt;
  //flag to allow TileDB LZ4 compression for GT fields
  bool m_enable_lz4_compression_gt;
  //flag that causes the loader to fail if this is an update (rather than a fresh load)
  bool m_fail_if_updating;
  //consolidate TileDB array after load - merges fragments
  bool m_consolidate_tiledb_array_after_load;
  //Discard entries with ./. or .|. as the GT field
  bool m_discard_missing_GTs;
  //The array will NOT contain mandatory VCF fields (ref, alt, qual, filter)
  //if this flag is enabled
  bool m_no_mandatory_VCF_fields;
 protected:
  void fix_callset_row_idx_bounds(const int rank);
};

#endif
