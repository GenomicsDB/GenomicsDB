/**
 * @file genomicsdb.h
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2019-2020 Omics Data Automation, Inc.
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
 * Header file to an experimental general GenomicsDB query interface.
 *
 **/

#ifndef GENOMICSDB_H
#define GENOMICSDB_H

#include "genomicsdb_exception.h"
#include "variant_query_config.h"
#include "tiledb.h"
#include "tiledb_utils.h"
#include "genomicsdb_logger.h"

#include <map>
#include <set>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdint.h>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <vector>
#include <mpi.h>
#include <functional>

// Override project visibility set to hidden for api
#if (defined __GNUC__ && __GNUC__ >= 4) || defined __INTEL_COMPILER
#  define GENOMICSDB_EXPORT __attribute__((visibility("default")))
#else
#  define GENOMICSDB_EXPORT
#endif

GENOMICSDB_EXPORT std::string genomicsdb_version();

typedef std::pair<uint64_t, uint64_t> interval_t;

typedef struct genomic_interval_t {
  std::string contig_name;
  interval_t interval;
  genomic_interval_t(std::string contig_name, interval_t interval) {
    this->contig_name = contig_name;
    this->interval = interval;
  }
} genomic_interval_t;

typedef struct genomic_field_type_t {
  enum genomicsdb_basic_type {INT32, FLOAT, CHAR, STRING, UNKNOWN};
  genomicsdb_basic_type type_idx = CHAR;

  //std::type_index type_idx = std::type_index(typeid(char));
  bool is_fixed_num_elements = true;
  size_t num_elements;
  size_t num_dimensions;
  bool contains_phase_info;
  genomic_field_type_t(std::type_index type_idx, bool is_fixed_num_elements, size_t num_elements, size_t num_dimensions, bool contains_phase_info) {
    this->type_idx = type_index_to_basic_type(type_idx);
    this->is_fixed_num_elements = is_fixed_num_elements;
    this->num_elements = num_elements;
    this->num_dimensions = num_dimensions;
    this->contains_phase_info = contains_phase_info;
  }
  genomic_field_type_t(genomicsdb_basic_type type_idx, bool is_fixed_num_elements, size_t num_elements, size_t num_dimensions, bool contains_phase_info) {
    this->type_idx = type_idx;
    this->is_fixed_num_elements = is_fixed_num_elements;
    this->num_elements = num_elements;
    this->num_dimensions = num_dimensions;
    this->contains_phase_info = contains_phase_info;
  }
  inline bool is_int() const {
    return type_idx == INT32;
  }
  inline bool is_float() const {
    return type_idx == FLOAT;
  }
  inline bool is_char() const {
    return type_idx == CHAR && is_fixed_num_elements;
  }
  inline bool is_string() const {
    return type_idx == CHAR && !is_fixed_num_elements;
  }
  inline bool is_cppstring() const {
    return type_idx == STRING;
  }
  inline bool contains_phase_information() const {
    return contains_phase_info;
  }
  static inline genomicsdb_basic_type type_index_to_basic_type(std::type_index ti) {
    if(ti == std::type_index(typeid(int)))
      return INT32;
    if(ti == std::type_index(typeid(float)))
      return FLOAT;
    if(ti == std::type_index(typeid(char)))
      return CHAR;
    return UNKNOWN;
  }
} genomic_field_type_t;

typedef struct genomic_field_t {
  std::string name;
  const void* ptr;
  size_t num_elements;
  genomic_field_t(std::string name, const void* ptr, size_t num_elements) {
    this->name = name;
    this->ptr = ptr;
    this->num_elements = num_elements;
  }
  size_t get_num_elements() const {
    return num_elements;
  }
  inline void check_offset(uint64_t offset) const {
     if (offset >= num_elements) {
       throw GenomicsDBException("Genomic Field="+name+" offset="+std::to_string(offset)+" greater than number of elements");
    }
  }
  inline int int_value_at(uint64_t offset) const {
    check_offset(offset);
    return *(reinterpret_cast<int *>(const_cast<void *>(ptr)) + offset);
  }
  inline float float_value_at(uint64_t offset) const {
    check_offset(offset);
    return *(reinterpret_cast<float *>(const_cast<void *>(ptr)) + offset);
  }
  inline char char_value_at(uint64_t offset) const {
    check_offset(offset);
    return *(reinterpret_cast<char *>(const_cast<void *>(ptr)) + offset);
  }
  inline std::string str_value() const {
    return std::string(reinterpret_cast<const char *>(ptr)).substr(0, num_elements);
  }
  inline std::string cpp_str_value_at(uint64_t offset) const {
    check_offset(offset);
    return *(reinterpret_cast<std::string *>(const_cast<void *>(ptr)) + offset);
  }
  GENOMICSDB_EXPORT std::string recombine_ALT_value(const genomic_field_type_t& field_type, std::string separator=", ") const;
  GENOMICSDB_EXPORT std::string combine_GT_vector(const genomic_field_type_t& field_type) const;
  std::string to_string(uint64_t offset, const genomic_field_type_t& field_type) const {
    if (field_type.is_int()) {
      return std::to_string(int_value_at(offset));
    } else if (field_type.is_float()) {
      return std::to_string(float_value_at(offset));
    } else if (field_type.is_char()) {
      return std::to_string(char_value_at(offset));
    } else if(field_type.is_cppstring()) {
      return cpp_str_value_at(offset);
    } else {
      return "";
    }
  }
  std::string to_string(const genomic_field_type_t& field_type, std::string separator=", ") const {
    if (field_type.is_string()) {
      if (name.compare("ALT") == 0) {
        return recombine_ALT_value(field_type);
      }
      return str_value();
    } else if (num_elements == 1) {
      return to_string(0, field_type);
    } else if (name.compare("GT") == 0) {
      return combine_GT_vector(field_type);
    } else {
      std::string output;
      for (auto i=0ul; i<num_elements; i++) {
        output = output + to_string(i, field_type);
        if (i<num_elements-1) {
          output = output + separator;
        }
      }
      return "[" + output + "]";
    }
  }
} genomic_field_t;


// Opaque data structure to Variant - see src/main/cpp/include/genomicsdb/variant.h
// Similar to GAVariant in GA4GH API
typedef struct genomicsdb_variant_t genomicsdb_variant_t;

// Opaque data structure to VariantCall
// Similar to GACall in GA4GH API. Stores info about 1 CallSet/row for a given position
typedef struct genomicsdb_variant_call_t genomicsdb_variant_call_t;

// Data structure to query a subset of ranges, can represent either column or row ranges
// TODO: This should change to uint64_t
typedef std::vector<std::pair<int64_t, int64_t>> genomicsdb_ranges_t;

// Default Segment Size in bytes = 10MB
#define DEFAULT_SEGMENT_SIZE 10u*1024u*1024u
#define SCAN_FULL {{0, INT64_MAX-1}}
#define ALL_ATTRIBUTES {}

template<typename T>
class GenomicsDBResults {
 public:
  GenomicsDBResults(std::vector<T>* results, std::map<std::string, genomic_field_type_t> genomic_field_types)
      : m_results(results), m_current_pos(0),
        m_genomic_field_types(std::make_shared<std::map<std::string, genomic_field_type_t>>(std::move(genomic_field_types))) {};
  ~GenomicsDBResults() { free(); };
  const std::shared_ptr<std::map<std::string, genomic_field_type_t>> get_genomic_field_types() const {
    return m_genomic_field_types;
  }
  const genomic_field_type_t get_genomic_field_type(const std::string& name) const {
    if (m_genomic_field_types->find(name) != m_genomic_field_types->end()) {
      return m_genomic_field_types->at(name);
    } else {
      throw GenomicsDBException("Genomic Field="+name+" does not seem to have an associated type");
    }
  }
  GENOMICSDB_EXPORT std::size_t size() const noexcept;
  GENOMICSDB_EXPORT const T* at(std::size_t pos);
  inline const T* next() { return at(m_current_pos++); };
 private:
  GENOMICSDB_EXPORT void free();
  std::vector<T>* m_results;
  std::size_t m_current_pos;
  std::shared_ptr<std::map<std::string, genomic_field_type_t>> m_genomic_field_types;
};

// Specializations for the GenomicsDBResults Template
typedef GenomicsDBResults<genomicsdb_variant_t> GenomicsDBVariants;
typedef GenomicsDBResults<genomicsdb_variant_call_t> GenomicsDBVariantCalls;

class GENOMICSDB_EXPORT GenomicsDBVariantCallProcessor {
 public:
  GenomicsDBVariantCallProcessor() {};
  void initialize(std::map<std::string, genomic_field_type_t> genomic_field_types) {
    m_genomic_field_types = std::make_shared<std::map<std::string, genomic_field_type_t>>(std::move(genomic_field_types));
  }
  const std::shared_ptr<std::map<std::string, genomic_field_type_t>> get_genomic_field_types() const {
    return m_genomic_field_types;
  }
  const genomic_field_type_t get_genomic_field_type(const std::string& name) const {
    if (m_genomic_field_types->find(name) != m_genomic_field_types->end()) {
      return m_genomic_field_types->at(name);
    } else {
      throw GenomicsDBException("Genomic Field="+name+" does not seem to have an associated type");
    }
  }
  virtual void process(const interval_t& interval);
  virtual void process(const std::string& sample_name,
                       const int64_t* coordinates,
                       const genomic_interval_t& genomic_interval,
                       const std::vector<genomic_field_t>& genomic_fields);
 private:
  std::shared_ptr<std::map<std::string, genomic_field_type_t>> m_genomic_field_types;
};

class Variant;

class GENOMICSDB_EXPORT GenomicsDBVariantProcessor {
 public:
  GenomicsDBVariantProcessor() {};
  void initialize(std::map<std::string, genomic_field_type_t> genomic_field_types) {
     m_genomic_field_types = std::make_shared<std::map<std::string, genomic_field_type_t>>(std::move(genomic_field_types));
  }
  const std::shared_ptr<std::map<std::string, genomic_field_type_t>> get_genomic_field_types() const {
    return m_genomic_field_types;
  }
  const genomic_field_type_t get_genomic_field_type(const std::string& name) const {
    if (m_genomic_field_types->find(name) != m_genomic_field_types->end()) {
      return m_genomic_field_types->at(name);
    } else {
      throw GenomicsDBException("Genomic Field="+name+" does not seem to have an associated type");
    }
  }
  virtual void process(const Variant& variant) = 0;
  virtual void process(const std::vector<Variant>& variants);
 private:
  std::shared_ptr<std::map<std::string, genomic_field_type_t>> m_genomic_field_types;
};

class GENOMICSDB_EXPORT GenomicsDBPlinkProcessor : public GenomicsDBVariantProcessor {
  public:
    // the formats variable encodes which formats to produce by the values of specific bits: 0 - bgen, 1 - bed, 2 - tped, e.g. formats == 0b110 encodes bed and tped
    GenomicsDBPlinkProcessor(VariantQueryConfig* qc,
                             VidMapper* vid_mapper,
                             const std::string& array,
                             unsigned char formats = 7,
                             int compression = 1,
                             bool verbose = false,
                             double progress_interval = -1,
                             std::string prefix = "output",
                             std::string fam_list = "",
                             int rank = 0);

    ~GenomicsDBPlinkProcessor() {
      TileDBUtils::finalize_codec(codec);
    }

    virtual void process(const interval_t& interval);
    virtual void process(const std::string& sample_name,
                         const int64_t* coordinates,
                         const genomic_interval_t& genomic_interval,
                         const std::vector<genomic_field_t>& genomic_fields,
                         const bool phased);
    virtual void process(const Variant& variant) override;
    void advance_state();
    const char BGEN_MASK = 1;
    const char  BED_MASK = 2;
    const char TPED_MASK = 4;
  private:
    const std::string& array;
    VidMapper* vid_mapper;
    bool make_bgen, make_tped, make_bed;
    int compression = 0; // 0 for none, 1 for zlib, 2 for zstd
    bool verbose = false;
    // flattened coordinate to place in sorted map, phased status of column for bgen purposes (entire column considered unphased if any are unphased)
    std::map<uint64_t, std::pair<int, bool>> variant_map;
    double progress_interval;
    std::string fam_list;
    std::string prefix;
    VariantQueryConfig* query_config;
    // row to place in sorted map and sample name
    std::map<uint64_t, std::pair<int, std::string>> sample_map;
    bool sample_map_initialized = false;
    // fam is identical to tfam, used with bed, tped respectively
    std::fstream tped_file, fam_file, bim_file, bed_file, bgen_file;
    int state = 0;
    int last_sample = -1;
    int num_variants = 0;
    int last_coord = -1;
    int last_alleles = -1;
    bool last_phased;
    int rank;
    int total_rows = 0;
    int total_cols = 0;

    size_t tm = 0;
    void progress_bar(const int64_t* coords) {
      size_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      if (now - progress_interval > tm) {
        // flatten column and row ranges
        int row = 0, col = 0;
        for(auto& a : query_config->get_query_row_ranges(rank)) {
          if(coords[0] <= a.second) {
            row += coords[0] - a.first + 1;
            break;
          }
          else {
            row += a.second - a.first + 1;
          }
        }
        for(auto& b : query_config->get_query_column_ranges(rank)) {
          if(coords[1] <= b.second) {
            col += coords[1] - b.first + 1;
            break;
          }
          else {
            col += b.second - b.first + 1;
          }
        }

        long num = (long)col * total_rows + row + ((bool)state)*((long)total_rows * total_cols);
        long den = (long)total_rows * total_cols * 2;

        logger.info("Plink progress {} / {} = {:.2f}%", num, den, 100*(double)num/den);

        tm = now;
      }
    }

    // populate variant data block in bgen file for variant with given rsid spanning genomic_interval
    // vec is REF and ALT combined (REF is first entry)
    void bgen_variant_data_block(const std::string& rsid, const genomic_interval_t& genomic_interval, const std::vector<std::string>& vec, bool phased);

    // return vector of tokens that were separated by sep in str
    std::vector<std::string> split(std::string str, std::string sep = ",") {
      std::vector<std::string> retval;
      int index;
      if(str.length() >= 2) {
        if(str[0] == '[') {
           str = str.substr(1, str.length() - 2);
        }
      }
      while((index = str.find(sep)) != std::string::npos) {
        retval.push_back(str.substr(0, index));
        str.erase(0, index + 1);
      }
      retval.push_back(str);
      return retval;
    };

    // BED variables/functions
    char bed_buf = 0;
    char bed_buf_state = 0;
    void flush_to_bed() {
      if(bed_buf_state) {
        bed_file.write(&bed_buf, 1);
        bed_buf = 0;
        bed_buf_state = 0;
      }
    }

    void write_to_bed(char x) {
      bed_buf += x << bed_buf_state * 2;
      ++bed_buf_state %= 4;    
      if(!bed_buf_state){
        bed_file.write(&bed_buf, 1);
        bed_buf = 0;
      }
    }
    // BGEN variables
    char min_ploidy, max_ploidy;
    int32_t bgen_gt_size;
    uint32_t samples_in_column = 0;
    void *codec;
    std::string codec_buf;

    // FIXME: hard coded for B = 8
    // callback expects GT vector
    void bgen_enumerate_phased(int ploidy, int alleles, std::function<void(const std::vector<int>&, int)> callback, bool drop_last = true) {
      int size = 0;
      int ind = -1;
      int limit = alleles - drop_last;
      for(int i = 0; i < ploidy; i++) {
        for(int j = 0; j < limit; j++) {
          callback({i, j}, ++ind);
        }
      }
    }

    // get_probabilitiles expects allele counts
    void bgen_enumerate_unphased(int ploidy, int alleles, std::function<void(const std::vector<int>&, size_t)> callback, bool drop_last = true) {
      int size = 0;
      int ind = -1;
      std::vector<int> allele_counts(alleles);

      std::function<void(int, int)> enumerate_unphased;
      enumerate_unphased = [&] (int used, int depth) {
        int limit = ploidy - used - ((depth == alleles - 1) && drop_last); // if the highest depth (rightmost) do not iterate to highest possible count in order to drop last ((0, 0, ..., X) is the last)
        if(depth) {
          for(int i = 0; i <= limit; i++) { 
            allele_counts[depth] = i;
            enumerate_unphased(used + i, depth - 1);
          }
        }
        else {
          allele_counts[depth] = ploidy - used;
          ++ind;
          callback(allele_counts, (size_t)ind);
        }
      };

      enumerate_unphased(0, alleles - 1);    
    }

    void bgen_empty_cell(int ploidy, int alleles, bool phased) {
      auto write_zero = [&] (const std::vector<int>& v, int) {
        char z = 0;
        codec_buf.push_back(0);
      };
      if(phased) {
        bgen_enumerate_phased(ploidy, alleles, write_zero);
      }
      else {
        bgen_enumerate_unphased(ploidy, alleles, write_zero);
      }
    }

    // fill in size, min/max ploidy of last column
    void bgen_finish_gt() {
      // write min ploidy
      codec_buf[7] = min_ploidy;      

      // write max ploidy
      codec_buf[8] = max_ploidy;

      size_t uncompressed_size = codec_buf.size(), data_size = codec_buf.size();

      if(compression) {
        char* data;
        TileDBUtils::compress(codec, (unsigned char*)codec_buf.c_str(), codec_buf.length(), (void**)&data, data_size);
        size_t total_size = data_size + 4;
        bgen_file.write((char*)&total_size, 4); // BGEN: size of previous gt probability data plus D field
        bgen_file.write((char*)&uncompressed_size, 4);
        bgen_file.write(data, data_size);
      }
      else {
        bgen_file.write((char*)&uncompressed_size, 4); // BGEN: size of previous gt probability data plus D field
        bgen_file.write(codec_buf.c_str(), codec_buf.length());
      }

      codec_buf.clear();
      bgen_gt_size = 0;
      min_ploidy = 64;
      max_ploidy = -1;
    }

    // locations in file
    int bgen_gt_size_offset;
    int bgen_min_ploidy_offset;
    int bgen_max_ploidy_offset;
    int bgen_ploidy_info_offset;
    int bgen_probability_offset;
};

// Forward Declarations for keeping Variant* classes opaque
class Variant;
class VariantCall;
class VariantQueryConfig;

/**
 * Experimental Query Interface to GenomicsDB for Arrays partitioned by columns
 * Concurrency support is provided via query json files for now - see 
 *     https://github.com/GenomicsDB/GenomicsDB/wiki/Querying-GenomicsDB#json-configuration-file-for-a-query
 *     https://github.com/GenomicsDB/GenomicsDB/wiki/MPI-with-GenomicsDB
 */
class GenomicsDB {
 public:
  enum GENOMICSDB_EXPORT query_config_type_t { JSON_FILE=0, JSON_STRING=1, PROTOBUF_BINARY_STRING=2 };

  /**
   * Constructor to the GenomicsDB Query API
   *   workspace
   *   callset_mapping_file
   *   vid_mapping_file
   *   reference_genome
   *   attributes, optional 
   *   segment_size, optional
   * Throws GenomicsDBException
   */
  GENOMICSDB_EXPORT GenomicsDB(const std::string& workspace,
             const std::string& callset_mapping_file,
             const std::string& vid_mapping_file,
             const std::string& reference_genome,
             const std::vector<std::string>attributes = ALL_ATTRIBUTES,
             const uint64_t segment_size = DEFAULT_SEGMENT_SIZE);

  /**
   * Constructor to the GenomicsDB Query API with configuration json files
   *   query_configuration - describe the query configuration in either a JSON file or string
   *   query_configuration_type - type of query configuration, could be a JSON_FILE or JSON_STRING
   *   loader_config_json_file, optional - describe the loader configuration in a JSON file.
   *           If a configuration key exists in both the query and the loader configuration, the query
   *           configuration takes precedence
   *   concurrency_rank, optional - if greater than 0, 
   *           the constraints(workspace, array, column and row ranges) are surmised
   *           using the rank as an index into their corresponding vectors
   * Throws GenomicsDBException
   */
  GENOMICSDB_EXPORT GenomicsDB(const std::string& query_configuration,
                               const query_config_type_t query_configuration_type = JSON_FILE,
                               const std::string& loader_configuration_json_file = std::string(),
                               const int concurrency_rank=0);
  /**
   * Destructor
   */
  GENOMICSDB_EXPORT ~GenomicsDB();

  /**
   * Query GenomicsDB array for variants constrained by column and row ranges.
   * Variants are similar to GAVariant in GA4GH API
   *   array
   *   column_ranges, optional
   *   row_ranges, optional
   */
  GENOMICSDB_EXPORT GenomicsDBVariants query_variants(const std::string& array,
                                                      genomicsdb_ranges_t column_ranges=SCAN_FULL,
                                                      genomicsdb_ranges_t row_ranges={});

  /**
   * Query using set configuration for variants. Useful when using parallelism paradigms(MPI, Intel TBB)
   * Variants are similar to GAVariant in GA4GH API
   */
  GENOMICSDB_EXPORT GenomicsDBVariants  query_variants();

  /**
   * Query the array for variant calls constrained by the column and row ranges.
   * Variant Calls are similar to GACall in GA4GH API.
   *   array
   *   column_ranges, optional
   *   row_ranges, optional
   */
  GENOMICSDB_EXPORT GenomicsDBVariantCalls query_variant_calls(const std::string& array,
                                                               genomicsdb_ranges_t column_ranges=SCAN_FULL,
                                                               genomicsdb_ranges_t row_ranges={});

  /**
   * Query the array for variant calls constrained by the column and row ranges.
   * Variant Calls are similar to GACall in GA4GH API.
   *   array
   *   column_ranges, optional
   *   row_ranges, optional
   */
  GENOMICSDB_EXPORT GenomicsDBVariantCalls query_variant_calls(GenomicsDBVariantCallProcessor& processor,
                                                               const std::string& array,
                                                               genomicsdb_ranges_t column_ranges=SCAN_FULL,
                                                               genomicsdb_ranges_t row_ranges={});
 
  /**
   * Query using set configuration for variant calls. Useful when using parallelism paradigms(MPI, Intel TBB)
   * Variant Calls are similar to GACall in GA4GH API.
   */
  GENOMICSDB_EXPORT GenomicsDBVariantCalls query_variant_calls();

  /**
   * Query using set configuration for variant calls. Useful when using parallelism paradigms(MPI, Intel TBB)
   * Variant Calls are similar to GACall in GA4GH API.
   */
  GENOMICSDB_EXPORT GenomicsDBVariantCalls query_variant_calls(GenomicsDBVariantCallProcessor& processor);

  GENOMICSDB_EXPORT void generate_vcf(const std::string& array,
                                 genomicsdb_ranges_t column_ranges,
                                 genomicsdb_ranges_t row_ranges,
                                 const std::string& output = "",
                                 const std::string& output_format = "",
                                 bool overwrite = false);

  GENOMICSDB_EXPORT void generate_vcf(const std::string& output = "",
                                      const std::string& output_format = "",
                                      bool overwrite = false);

  /**
   * Query by column and row ranges, use results to generate plink .ped and .map files
   * The two files will be named <output_prefix>.ped and <output_prefix>.map
   */
  GENOMICSDB_EXPORT void generate_plink(const std::string& array,
                                        VariantQueryConfig* query_config,
                                        unsigned char format = 7,
                                        int compression = 1,
                                        bool one_pass = false,
                                        bool verbose = false,
                                        double progress_interval = -1,
                                        const std::string& output_prefix = "output",
                                        const std::string& fam_list = "");

  /**
   * Utility template functions to extract information from Variant and VariantCall classes
   */
  GENOMICSDB_EXPORT interval_t get_interval(const genomicsdb_variant_t* variant);
  GENOMICSDB_EXPORT interval_t get_interval(const genomicsdb_variant_call_t* variant_call);

  GENOMICSDB_EXPORT genomic_interval_t get_genomic_interval(const genomicsdb_variant_t* variant);
  GENOMICSDB_EXPORT genomic_interval_t get_genomic_interval(const genomicsdb_variant_call_t* variant_call);

  GENOMICSDB_EXPORT std::vector<genomic_field_t> get_genomic_fields(const std::string& array, const genomicsdb_variant_t* variant);
  GENOMICSDB_EXPORT std::vector<genomic_field_t> get_genomic_fields(const std::string& array, const genomicsdb_variant_call_t* variant_call);

  GENOMICSDB_EXPORT GenomicsDBVariantCalls get_variant_calls(const std::string& array, const genomicsdb_variant_t* variant);

  GENOMICSDB_EXPORT int64_t get_row(const genomicsdb_variant_call_t* variant_call);
  
 private:
  std::vector<Variant>*  query_variants(const std::string& array,
                                        VariantQueryConfig *query_config);
  void query_variants(const std::string& array,
                      VariantQueryConfig* query_config,
                      GenomicsDBVariantProcessor& proc);
  std::vector<VariantCall>* query_variant_calls(const std::string& array,
                                                VariantQueryConfig *query_config,
                                                GenomicsDBVariantCallProcessor& processor);
  void generate_vcf(const std::string& array,
                    VariantQueryConfig* query_config,
                    const std::string& output,
                    const std::string& output_format,
                    bool overwrite);

  /**
   * Utility function to set the type of the ALT field to string
   */ 
  void change_alt_to_string(std::map<std::string, genomic_field_type_t>& types);

  VariantQueryConfig* get_query_config_for(const std::string& array);

  void* m_storage_manager = nullptr; // Pointer to VariantStorageManager instance
  void* m_vid_mapper = nullptr;      // Pointer to VidMapper instance
  void* m_query_config = nullptr;    // Pointer to a base VariantQueryConfig instance
  void *m_annotation_service = nullptr; // Pointer to AnnotationService instance

  int m_concurrency_rank = 0;

  //TODO: Get VariantFields to have names instead of indices at the point of building the data structure, so
  //      we don't have to maintain m_query_configs_map
  // Associate array names with VariantQueryConfig 
  std::map<std::string, VariantQueryConfig> m_query_configs_map; 
};

// genomicsdb_variant_t specialization of GenomicsDBResults template
template<>
GENOMICSDB_EXPORT std::size_t GenomicsDBResults<genomicsdb_variant_t>::size() const noexcept;
template<>
GENOMICSDB_EXPORT const genomicsdb_variant_t* GenomicsDBResults<genomicsdb_variant_t>::at(std::size_t pos);
template<>
GENOMICSDB_EXPORT void GenomicsDBResults<genomicsdb_variant_t>::free();

// genomicsdb_variant_call_t specialization for the GenomicsDBResults template
template<>
GENOMICSDB_EXPORT std::size_t GenomicsDBResults<genomicsdb_variant_call_t>::size() const noexcept;
template<>
GENOMICSDB_EXPORT const genomicsdb_variant_call_t* GenomicsDBResults<genomicsdb_variant_call_t>::at(std::size_t pos);
template<>
GENOMICSDB_EXPORT void GenomicsDBResults<genomicsdb_variant_call_t>::free();


#endif /* GENOMICSDB_H */
