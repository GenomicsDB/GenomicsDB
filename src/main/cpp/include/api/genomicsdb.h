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

#include <map>
#include <memory>
#include <sstream>
#include <stdint.h>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <vector>

// Override project visibility set to hidden for api
#if (defined __GNUC__ && __GNUC__ >= 4) || defined __INTEL_COMPILER
#  define GENOMICSDB_EXPORT __attribute__((visibility("default")))
#else
#  define GENOMICSDB_EXPORT
#endif

GENOMICSDB_EXPORT std::string genomicsdb_version();

#define PHASED_ALLELE_SEPARATOR '|';
#define UNPHASED_ALLELE_SEPARATOR '/';
#define ALLELE_SEPARATOR '|'
#define NON_REF_VARIANT "&"
#define NON_REF "<NON_REF>"

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
  std::type_index type_idx = std::type_index(typeid(char));
  bool is_fixed_num_elements = true;
  size_t num_elements;
  size_t num_dimensions;
  bool contains_phase_info;
  genomic_field_type_t(std::type_index type_idx, bool is_fixed_num_elements, size_t num_elements, size_t num_dimensions, bool contains_phase_info) {
    this->type_idx = type_idx;
    this->is_fixed_num_elements = is_fixed_num_elements;
    this->num_elements = num_elements;
    this->num_dimensions = num_dimensions;
    this->contains_phase_info = contains_phase_info;
  }
  inline bool is_int() const {
    return type_idx == std::type_index(typeid(int));
  }
  inline bool is_float() const {
    return type_idx == std::type_index(typeid(float));
  }
  inline bool is_char() const {
    return (type_idx == std::type_index(typeid(char))) && is_fixed_num_elements;
  }
  inline bool is_string() const {
    return (type_idx == std::type_index(typeid(char))) && !is_fixed_num_elements;
  }
  inline bool contains_phase_information() const {
    return contains_phase_info;
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
  std::string recombine_ALT_value(char delim=ALLELE_SEPARATOR, std::string separator=", ") const {
    std::stringstream ss(str_value());
    std::string output;
    std::string item;
    while (std::getline(ss, item, delim)){
      if (item.compare(NON_REF_VARIANT) == 0) {
        output.empty()?output=NON_REF:output=output+separator+NON_REF;
      } else {
        output.empty()?output=item:output=output+separator+item;
      }
    }
    return "[" + output + "]";
  }
  std::string combine_GT_vector(const genomic_field_type_t& field_type) const {
    std::string output;
    if (field_type.contains_phase_information()) {
      for (auto i=0ul; i<num_elements; i++) {
        output = output + to_string(i, field_type);
        if (i+1<num_elements) {
          if (to_string(i+1, field_type).compare("0") == 0) {
            output = output + UNPHASED_ALLELE_SEPARATOR;
          } else {
            output = output + PHASED_ALLELE_SEPARATOR;
          }
          i++;
        }
      }
    } else {
      for (auto i=0ul; i<num_elements; i++) {
        output = output + to_string(i, field_type);
        if (i+1<num_elements) {
          output = output + UNPHASED_ALLELE_SEPARATOR;
        }
      }
    }
    return output;
  }
  std::string to_string(uint64_t offset, const genomic_field_type_t& field_type) const {
    if (field_type.is_int()) {
      return std::to_string(int_value_at(offset));
    } else if (field_type.is_float()) {
      return std::to_string(float_value_at(offset));
    } else if (field_type.is_char()) {
      return std::to_string(char_value_at(offset));
    } else {
      return "";
    }
  }
  std::string to_string(const genomic_field_type_t& field_type, std::string separator=", ") const {
    if (field_type.is_string()) {
      if (name.compare("ALT") == 0) {
        return recombine_ALT_value();
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
  std::vector<VariantCall>* query_variant_calls(const std::string& array,
                                                VariantQueryConfig *query_config,
                                                GenomicsDBVariantCallProcessor& processor);
  void generate_vcf(const std::string& array,
                    VariantQueryConfig* query_config,
                    const std::string& output,
                    const std::string& output_format,
                    bool overwrite);

  VariantQueryConfig* get_query_config_for(const std::string& array);

  void* m_storage_manager = nullptr; // Pointer to VariantStorageManager instance
  void* m_vid_mapper = nullptr;      // Pointer to VidMapper instance
  void* m_query_config = nullptr;    // Pointer to a base VariantQueryConfig instance

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
