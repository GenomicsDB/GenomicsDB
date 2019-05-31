/**
 * @file genomicsdb.cc
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
 * Implementation of the query interface to GenomicsDB
 *
 **/

#include "genomicsdb.h"

#include <iostream>
#include <string>
#include "query_variants.h"
#include "timer.h"
#include "variant_field_data.h"
#include "variant_query_config.h"
#include "vid_mapper_pb.h"

#define TO_VARIANT_QUERY_CONFIG(X) (reinterpret_cast<VariantQueryConfig *>(static_cast<void *>(X)))
#define TO_VARIANT_STORAGE_MANAGER(X) (reinterpret_cast<VariantStorageManager *>(static_cast<void *>(X)))
#define TO_VID_MAPPER(X) (reinterpret_cast<VidMapper *>(static_cast<void *>(X)))

#define TO_VARIANT(X) (reinterpret_cast<const Variant *>(static_cast<const void *>(X)))
#define TO_VARIANT_CALL(X) (reinterpret_cast<const VariantCall *>(static_cast<const void *>(X)))

#define TO_GENOMICSDB_VARIANT_T_VECTOR(X) (reinterpret_cast<std::vector<genomicsdb_variant_t> *>(static_cast<void *>(X)))
#define TO_GENOMICSDB_VARIANT_CALL_T_VECTOR(X) (reinterpret_cast<std::vector<genomicsdb_variant_call_t> *>(static_cast<void *>(X))

std::string genomicsdb_version() {
  return GENOMICSDB_VERSION;
}

#define VERIFY(X) if(!(X)) throw GenomicsDBException(#X);

inline void check(const std::string& workspace,
                  const uint64_t segment_size,
                  const std::string& callset_mapping_file,
                  const std::string& vid_mapping_file,
                  const std::string& reference_genome) {
  VERIFY(!workspace.empty() && "Workspace specified cannot be empty");
  VERIFY(segment_size && "Segment size specified has to be greater than 0");
  VERIFY(!callset_mapping_file.empty() && "Callset Mapping File cannot be empty");
  VERIFY(!vid_mapping_file.empty() && "Vid Mapping File cannot be empty");
  VERIFY(!reference_genome.empty() && "Reference Genome cannot be empty for querying");
}

GenomicsDB::GenomicsDB(const std::string& workspace,
                       const std::string& callset_mapping_file,
                       const std::string& vid_mapping_file,
                       const std::string& reference_genome,
                       const std::vector<std::string>attributes,
                       const uint64_t segment_size) {
  check(workspace, segment_size, callset_mapping_file, vid_mapping_file, reference_genome);
    
  // Create base query configuration
  m_query_config = new VariantQueryConfig();

  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);

  query_config->set_workspace(workspace);
  query_config->set_callset_mapping_file(callset_mapping_file);
  query_config->set_vid_mapping_file(vid_mapping_file);
  query_config->set_reference_genome(reference_genome);
  query_config->set_attributes_to_query(attributes);

  // Create storage manager
  m_storage_manager = new VariantStorageManager(workspace, segment_size);

  // Read in vid mapping information
  m_vid_mapper = new FileBasedVidMapper(vid_mapping_file);
}

GenomicsDB::GenomicsDB(const std::string& query_configuration_json_file,
                       const std::string& loader_configuration_json_file,
                       const int concurrency_rank) : m_concurrency_rank(concurrency_rank) {
  // Create base query configuration
  m_query_config = new VariantQueryConfig();

  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);

  GenomicsDBImportConfig loader_config;
  if (!loader_configuration_json_file.empty()) {
    loader_config.read_from_file(loader_configuration_json_file, concurrency_rank);
    query_config->update_from_loader(loader_config, concurrency_rank);
  }

  query_config->read_from_file(query_configuration_json_file, concurrency_rank);

  check(query_config->get_workspace(concurrency_rank), query_config->get_segment_size(),
        query_config->get_callset_mapping_file(), query_config->get_vid_mapping_file(),
        query_config->get_reference_genome());

  // Discard intervals not part of this partition
  query_config->subset_query_column_ranges_based_on_partition(loader_config, concurrency_rank);

  // Create storage manager
  m_storage_manager = new VariantStorageManager(query_config->get_workspace(concurrency_rank), query_config->get_segment_size());

  // Read in vid mapping information
  m_vid_mapper = new FileBasedVidMapper(query_config->get_vid_mapping_file());
}

GenomicsDB::~GenomicsDB() {
  // TODO: delete variant_query_config per array
  delete TO_VARIANT_QUERY_CONFIG(m_query_config);
  delete TO_VARIANT_STORAGE_MANAGER(m_query_config);
  delete TO_VID_MAPPER(m_vid_mapper);
}

std::vector<genomicsdb_variant_t>* GenomicsDB::query_variants(const std::string& array,
                                                              const genomicsdb_ranges_t column_ranges,
                                                              const genomicsdb_ranges_t row_ranges) {
  // Create Variant Config for given concurrency_rank
  VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  VariantQueryConfig query_config(*base_query_config);
  query_config.set_array_name(array);
  query_config.set_query_column_ranges(column_ranges);
  query_config.set_query_row_ranges(row_ranges);
 
  query_config.validate(m_concurrency_rank);

  return TO_GENOMICSDB_VARIANT_T_VECTOR(query_variants(array, &query_config));
}

std::vector<genomicsdb_variant_t>* GenomicsDB::query_variants() {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);
  return TO_GENOMICSDB_VARIANT_T_VECTOR(query_variants(array, query_config));
}


GenomicsDBVariants GenomicsDB::query_variants1(const std::string& array,
                                           genomicsdb_ranges_t column_ranges,
                                           genomicsdb_ranges_t row_ranges) {
    // Create Variant Config for given concurrency_rank
  VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  VariantQueryConfig query_config(*base_query_config);
  query_config.set_array_name(array);
  query_config.set_query_column_ranges(column_ranges);
  query_config.set_query_row_ranges(row_ranges);

  return GenomicsDBResults<genomicsdb_variant_t>(TO_GENOMICSDB_VARIANT_T_VECTOR(query_variants(array, &query_config)));
}


std::vector<Variant>* GenomicsDB::query_variants(const std::string& array, VariantQueryConfig* query_config) {
#if(0)
  auto query_timer = Timer();
  query_timer.start();
#endif

  query_config->validate();

  // Get a Variant Query Processor
  VariantQueryProcessor *query_processor = new VariantQueryProcessor(
      TO_VARIANT_STORAGE_MANAGER(m_storage_manager), array, query_config->get_vid_mapper());
  query_processor->do_query_bookkeeping(query_processor->get_array_schema(),
                                        *query_config, query_config->get_vid_mapper(), true);

  // Perform Query over all the intervals
  std::vector<Variant>* pvariants = new std::vector<Variant>;
  //  std::vector<Variant> variants;
  for (auto i=0u; i<query_config->get_num_column_intervals(); i++) {
    query_processor->gt_get_column_interval(query_processor->get_array_descriptor(),
                                            *query_config, i, *pvariants);
  }

#if(1)
  print_variants(*pvariants, "default", *query_config);
#endif

  free(query_processor);

#if(0)
  query_timer.stop();
  query_timer.print();
#endif

  // TODO: Should not need this map, but required for printing out variant field names. 
  m_query_configs_map.emplace(array, *query_config);

  return pvariants;
}

std::vector<genomicsdb_variant_call_t>* GenomicsDB::query_variant_calls(const std::string& array,
                                                                        genomicsdb_ranges_t column_ranges,
                                                                        genomicsdb_ranges_t row_ranges) {
  // NYI
  std::cerr << "query_variant_calls NOT YET IMPLEMENTED" << std::endl;
  return nullptr;
}
  
std::vector<genomicsdb_variant_call_t>* GenomicsDB::query_variant_calls() {
  // NYI
  std::cerr << "query_variant_calls NOT YET IMPLEMENTED" << std::endl;
  return nullptr;
}

/*GenomicsDBVariantCalls GenomicsDB::query_variant_calls1(const std::string& array,
                                                        genomicsdb_ranges_t column_ranges,
                                                        genomicsdb_ranges_t row_ranges) {
  //
  }*/

std::vector<VariantCall>* query_variant_calls(VariantQueryConfig *query_config) {
  // NYI
  std::cerr << "query_variant_calls NOT YET IMPLEMENTED" << std::endl;
  return nullptr;
}

// Template to get the mapped interval from the GenomicsDB array for the Variant(Call)
template<class VariantOrVariantCall>
interval_t get_interval_for(const VariantOrVariantCall* variant_or_variant_call) {
  return std::make_pair<uint64_t, uint64_t>(variant_or_variant_call->get_column_begin(), variant_or_variant_call->get_column_end());
}

interval_t GenomicsDB::get_interval(const genomicsdb_variant_t* variant) {
  return get_interval_for(TO_VARIANT(variant));
}

interval_t GenomicsDB::get_interval(const genomicsdb_variant_call_t* variant_call) {
  return get_interval_for(TO_VARIANT_CALL(variant_call));
}

// Template to get the genomic interval from the GenomicsDB array associated with the Variant(Call)
template<class VariantOrVariantCall>
genomic_interval_t get_genomic_interval_for(const VariantOrVariantCall* variant_or_variant_call, VidMapper *vid_mapper) {
  std::string contig_name;
  int64_t contig_position;
  vid_mapper->get_contig_location(variant_or_variant_call->get_column_begin(), contig_name,
                                                   contig_position);
  contig_position++;
  return genomic_interval_t(std::move(contig_name),
                            std::make_pair(contig_position,
                                           contig_position+(variant_or_variant_call->get_column_end()
                                                            -variant_or_variant_call->get_column_begin())));
}

genomic_interval_t GenomicsDB::get_genomic_interval(const genomicsdb_variant_t* variant) {
  return get_genomic_interval_for(TO_VARIANT(variant), TO_VID_MAPPER(m_vid_mapper));
}

genomic_interval_t GenomicsDB::get_genomic_interval(const genomicsdb_variant_call_t* variant_call) {
  return get_genomic_interval_for(TO_VARIANT_CALL(variant_call), TO_VID_MAPPER(m_vid_mapper));
}

// Template to get genomic fields from the GenomicsDB array associated with the Variant(Call)
template<class VariantOrVariantCall>
const std::vector<std::unique_ptr<VariantFieldBase>>& get_fields(const VariantOrVariantCall* variant_or_variant_call);
template<>
inline const std::vector<std::unique_ptr<VariantFieldBase>>& get_fields(const Variant* variant) {
  return variant->get_common_fields();
}
template<>
inline const std::vector<std::unique_ptr<VariantFieldBase>>& get_fields(const VariantCall* variant_call) {
  return variant_call->get_all_fields();
}

template<class VariantOrVariantCall>
std::string get_query_attribute_name(const VariantOrVariantCall* variant_or_variant_call, VariantQueryConfig *query_config, uint index);
template<>
inline std::string get_query_attribute_name(const Variant* variant, VariantQueryConfig *query_config, uint index) {
  return query_config->get_query_attribute_name(variant->get_query_idx_for_common_field(index));
}
template<>
inline std::string get_query_attribute_name(const VariantCall* variant_call, VariantQueryConfig *query_config, uint index) {
  return query_config->get_query_attribute_name(index);
}

template<class VariantOrVariantCall>
std::vector<genomic_field_t> get_genomic_fields_for(const std::string& array, const VariantOrVariantCall* variant_or_variant_call, VariantQueryConfig* query_config) {
  std::vector<genomic_field_t> fields;
  auto index = 0u;
  for (const auto& field: get_fields(variant_or_variant_call)) {
    if (field->is_valid()) {
      std::stringstream field_value_stream;
      field->print(field_value_stream);
      genomic_field_t field_vec(get_query_attribute_name(variant_or_variant_call, query_config, index),
                                field_value_stream.str());
      fields.push_back(field_vec);
   }
    index++;
  }

  return fields;
}

VariantQueryConfig* GenomicsDB::get_query_config_for(const std::string& array) {
  auto query_config_for_array = m_query_configs_map.find(array);
  if (query_config_for_array != m_query_configs_map.end()) {
    return &query_config_for_array->second;
  } else {
    return TO_VARIANT_QUERY_CONFIG(m_query_config);
  }
}

std::vector<genomic_field_t> GenomicsDB::get_genomic_fields(const std::string& array, const genomicsdb_variant_t* variant) {
  return get_genomic_fields_for(array, TO_VARIANT(variant), get_query_config_for(array));
}

std::vector<genomic_field_t> GenomicsDB::get_genomic_fields(const std::string& array, const genomicsdb_variant_call_t* variant_call) {
  return get_genomic_fields_for(array, TO_VARIANT_CALL(variant_call), get_query_config_for(array));
}

/*GenomicsDBVariantCalls get_variant_calls(const Variant* variant) {
  
  }

#define TO_XXX(X) (reinterpret_cast<std::vector<genomicsdb_variant_call_t> *>(X))
GenomicsDBVariantCalls get_variant_calls(const genomicsdb_variant_t* variant) {
  const int a = 3;
  int *b = (int *)(&a);
  return get_variant_calls(TO_VARIANT(variant));
}
*/
#define TO_VARIANT_VECTOR(X) (reinterpret_cast<std::vector<Variant> *>(static_cast<void *>(X)))
#define TO_VARIANT_CALL_VECTOR(X) (reinterpret_cast<std::vector<VariantCall> *>(static_cast<void *>(X)))

#define TO_GENOMICSDB_T(X) (reinterpret_cast<const T *>(static_cast<const void *>(X)))

// Navigate GenomicsDBResults
template<class T>
std::size_t GenomicsDBResults<T>::size() const noexcept {
  return m_results->size();
}

template<class T>
const T* GenomicsDBResults<T>::at(std::size_t pos) {
  if (pos >= size()) {
    return nullptr;
  } else {
    return TO_GENOMICSDB_T(&m_results[pos]);
  }
}

template<class T>
void GenomicsDBResults<T>::free() {
  delete m_results;
}
