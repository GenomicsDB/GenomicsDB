/**
 * @file genomicsdb.cc
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
 * Implementation of the query interface to GenomicsDB
 *
 **/

#include "genomicsdb.h"

#include <iostream>
#include <string>

#include "broad_combined_gvcf.h"
#include "query_variants.h"
#include "tiledb_utils.h"
#include "timer.h"
#include "variant_field_data.h"
#include "variant_operations.h"
#include "variant_query_config.h"
#include "vcf_adapter.h"
#include "vid_mapper_pb.h"
#include "genomicsdb_logger.h"

#define TO_VARIANT_QUERY_CONFIG(X) (reinterpret_cast<VariantQueryConfig *>(static_cast<void *>(X)))
#define TO_VARIANT_STORAGE_MANAGER(X) (reinterpret_cast<VariantStorageManager *>(static_cast<void *>(X)))
#define TO_VID_MAPPER(X) (reinterpret_cast<VidMapper *>(static_cast<void *>(X)))

#define TO_VARIANT(X) (reinterpret_cast<const Variant *>(static_cast<const void *>(X)))
#define TO_VARIANT_CALL(X) (reinterpret_cast<const VariantCall *>(static_cast<const void *>(X)))

#define TO_GENOMICSDB_VARIANT(X) (reinterpret_cast<const genomicsdb_variant_t *>(static_cast<const void *>(X)))
#define TO_GENOMICSDB_VARIANT_CALL(X) (reinterpret_cast<const genomicsdb_variant_call_t *>(static_cast<const void *>(X)))

#define TO_VARIANT_VECTOR(X) (reinterpret_cast<std::vector<Variant> *>(static_cast<void *>(X)))
#define TO_VARIANT_CALL_VECTOR(X) (reinterpret_cast<std::vector<VariantCall> *>(static_cast<void *>(X)))

#define TO_GENOMICSDB_VARIANT_VECTOR(X) (reinterpret_cast<std::vector<genomicsdb_variant_t> *>(static_cast<void *>(X)))
#define TO_GENOMICSDB_VARIANT_CALL_VECTOR(X) (reinterpret_cast<std::vector<genomicsdb_variant_call_t> *>(static_cast<void *>(X)))


std::string genomicsdb_version() {
  return GENOMICSDB_VERSION;
}

#define VERIFY(X) if(!(X)) throw GenomicsDBException(#X);

void check(const std::string& workspace,
                  const uint64_t segment_size,
                  const std::string& reference_genome) {
  VERIFY(!workspace.empty() && "Workspace specified cannot be empty");
  VERIFY(segment_size && "Segment size specified has to be greater than 0");
  VERIFY(!reference_genome.empty() && "Reference Genome cannot be empty for querying");
}

void check(const std::string& workspace,
                  const uint64_t segment_size,
                  const std::string& callset_mapping_file,
                  const std::string& vid_mapping_file,
                  const std::string& reference_genome) {
  check(workspace, segment_size, reference_genome);
  VERIFY(!callset_mapping_file.empty() && "Callset Mapping File cannot be empty");
  VERIFY(!vid_mapping_file.empty() && "Vid Mapping File cannot be empty");
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
  TO_VID_MAPPER(m_vid_mapper)->parse_callsets_json(callset_mapping_file, true);
}

GenomicsDB::GenomicsDB(const std::string& query_configuration,
                       const query_config_type_t query_configuration_type,
                       const std::string& loader_configuration_json_file,
                       const int concurrency_rank) : m_concurrency_rank(concurrency_rank) {
  VERIFY(!query_configuration.empty() && "Specified query configuration cannot be empty");

  // Create base query configuration
  m_query_config = new VariantQueryConfig();

  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);

  GenomicsDBImportConfig loader_config;
  if (!loader_configuration_json_file.empty()) {
    loader_config.read_from_file(loader_configuration_json_file, concurrency_rank);
    query_config->update_from_loader(loader_config, concurrency_rank);
  }

  switch (query_configuration_type) {
    case JSON_FILE:
      query_config->read_from_file(query_configuration, concurrency_rank); break;
    case JSON_STRING:
      query_config->read_from_JSON_string(query_configuration, concurrency_rank); break;
    case PROTOBUF_BINARY_STRING:
      query_config->read_from_PB_binary_string(query_configuration, concurrency_rank); break;
    default:
      throw GenomicsDBException("Unsupported query configuration type specified to the GenomicsDB constructor");
  }

  check(query_config->get_workspace(concurrency_rank),
        query_config->get_segment_size(),
        query_config->get_reference_genome());

  // Discard intervals not part of this partition
  query_config->subset_query_column_ranges_based_on_partition(loader_config, concurrency_rank);

  // Create storage manager
  m_storage_manager = new VariantStorageManager(query_config->get_workspace(concurrency_rank), query_config->get_segment_size(), query_config->enable_shared_posixfs_optimizations());
}

GenomicsDB::~GenomicsDB() {
  // TODO: delete variant_query_config per array
  if (m_vid_mapper != nullptr) {
    delete TO_VID_MAPPER(m_vid_mapper);
  }
  if (m_query_config != nullptr) {
    delete TO_VARIANT_QUERY_CONFIG(m_query_config);
  }
  if (m_storage_manager != nullptr) {
    delete TO_VARIANT_STORAGE_MANAGER(m_storage_manager);
  }
}

std::map<std::string, genomic_field_type_t> create_genomic_field_types(const VariantQueryConfig &query_config) {
  std::map<std::string, genomic_field_type_t> genomic_field_types;
  for (auto i=0u; i<query_config.get_num_queried_attributes(); i++) {
    const std::string attribute_name = query_config.get_query_attribute_name(i);
    const FieldInfo* field_info = query_config.get_vid_mapper().get_field_info(attribute_name);
    if (field_info) {
      const FieldElementTypeDescriptor descr = field_info->get_vcf_type();
      assert(descr.get_num_elements_in_tuple() > 0);
      const std::type_index type_index = descr.get_tuple_element_type_index(0);
      const FieldLengthDescriptor length_descr = field_info->m_length_descriptor;
      genomic_field_types.insert(std::make_pair(
          attribute_name,
          genomic_field_type_t(type_index,
                               length_descr.is_fixed_length_field(),
                               length_descr.is_fixed_length_field()?length_descr.get_num_elements():0,
                               length_descr.get_num_dimensions(),
                               length_descr.is_length_ploidy_dependent()&&length_descr.contains_phase_information())));
    }
  }
  return genomic_field_types;
}

GenomicsDBVariants GenomicsDB::query_variants(const std::string& array,
                                               genomicsdb_ranges_t column_ranges,
                                               genomicsdb_ranges_t row_ranges) {
  
  // Create Variant Config for given concurrency_rank
  VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  VariantQueryConfig query_config(*base_query_config);
  query_config.set_array_name(array);
  if (column_ranges.size() > 0) {
    query_config.set_query_column_ranges(column_ranges);
  }
  if (row_ranges.size() > 0) {
    query_config.set_query_row_ranges(row_ranges);
  }

  query_config.validate();

  return GenomicsDBVariants(TO_GENOMICSDB_VARIANT_VECTOR(query_variants(array, &query_config)),
                            create_genomic_field_types(query_config));
}

GenomicsDBVariants GenomicsDB::query_variants() {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);
  return GenomicsDBVariants(TO_GENOMICSDB_VARIANT_VECTOR(query_variants(array, query_config)),
                            create_genomic_field_types(*query_config));
}


std::vector<Variant>* GenomicsDB::query_variants(const std::string& array, VariantQueryConfig* query_config) {
#if(0)
  auto query_timer = Timer();
  query_timer.start();
#endif

  // Get a Variant Query Processor
  VariantQueryProcessor *query_processor = new VariantQueryProcessor(
      TO_VARIANT_STORAGE_MANAGER(m_storage_manager), array, query_config->get_vid_mapper());
  query_processor->do_query_bookkeeping(query_processor->get_array_schema(),
                                        *query_config, query_config->get_vid_mapper(), true);

  // Perform Query over all the intervals
  std::vector<Variant>* pvariants = new std::vector<Variant>;

  for (auto i=0u; i<query_config->get_num_column_intervals(); i++) {
    query_processor->gt_get_column_interval(query_processor->get_array_descriptor(),
                                            *query_config, i, *pvariants);
  }

#if(DEBUG)
  //print_variants(*pvariants, "default", *query_config);
#endif

  delete query_processor;

#if(0)
  query_timer.stop();
  query_timer.print();
#endif

  // TODO: Should not need this map, but required for associating variant field names with indices in Variant
  m_query_configs_map.emplace(array, *query_config);
  
  return pvariants;
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls(const std::string& array,
                                                       genomicsdb_ranges_t column_ranges,
                                                       genomicsdb_ranges_t row_ranges) {
  GenomicsDBVariantCallProcessor processor;
  return query_variant_calls(processor, array, column_ranges, row_ranges);
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls(GenomicsDBVariantCallProcessor& processor,
                                                       const std::string& array,
                                                       genomicsdb_ranges_t column_ranges,
                                                       genomicsdb_ranges_t row_ranges) {
  // Create Variant Config for given concurrency_rank
  VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  VariantQueryConfig query_config(*base_query_config);
  query_config.set_array_name(array);
  query_config.set_query_column_ranges(column_ranges);
  query_config.set_query_row_ranges(row_ranges);

  query_config.validate();

  return GenomicsDBVariantCalls(TO_GENOMICSDB_VARIANT_CALL_VECTOR(query_variant_calls(array, &query_config, processor)),
                                create_genomic_field_types(query_config));
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls() {
  GenomicsDBVariantCallProcessor processor;
  return query_variant_calls(processor);
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls(GenomicsDBVariantCallProcessor& processor) {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);
  return GenomicsDBVariantCalls(TO_GENOMICSDB_VARIANT_CALL_VECTOR(query_variant_calls(array, query_config, processor)),
                                create_genomic_field_types(*query_config));
}

VidMapper* get_vid_mapper(void* vid_mapper, void* query_config) {
  if (vid_mapper)  {
    return TO_VID_MAPPER(vid_mapper);
  } else {
    return const_cast<VidMapper *>(&(TO_VARIANT_QUERY_CONFIG(query_config)->get_vid_mapper()));
  }
}

class GatherVariantCalls : public SingleCellOperatorBase {
 public:
  GatherVariantCalls(GenomicsDBVariantCallProcessor& variant_call_processor, const VariantQueryConfig& query_config) :
      SingleCellOperatorBase(), m_variant_call_processor(variant_call_processor), m_vid_mapper(query_config.get_vid_mapper()) {
    initialize(query_config);
  }
  void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
  void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
                                const VariantArraySchema& schema);
 private:
  void initialize(const VariantQueryConfig& query_config);
  GenomicsDBVariantCallProcessor& m_variant_call_processor;
  const VidMapper& m_vid_mapper;
  std::shared_ptr<std::map<std::string, genomic_field_type_t>> m_genomic_field_types;
};

void GatherVariantCalls::initialize(const VariantQueryConfig& query_config) {
  m_variant_call_processor.initialize(create_genomic_field_types(query_config));
};

void GatherVariantCalls::operate(VariantCall& variant_call,
                                 const VariantQueryConfig& query_config,
                                 const VariantArraySchema& schema) {
  // m_variant_calls.push_back(std::move(variant_call));
  std::cout << "TBD: In GatherVariantCalls::operate()" << std::endl;
}

void GatherVariantCalls::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell,
                                                   const VariantQueryConfig& query_config,
                                                   const VariantArraySchema& schema) {
  if (cell.at_new_query_column_interval()) {
    auto curr_query_column_interval_idx = cell.get_current_query_column_interval_idx();
    auto begin = std::max<int64_t>(query_config.get_column_begin(curr_query_column_interval_idx), 0);
    auto end = std::min<int64_t>(query_config.get_column_end(curr_query_column_interval_idx), INT64_MAX-1);
    // TODO: Change Column/RowRange to be uint64_t intervals instead of int64_t, so we don't have to typecast
    m_variant_call_processor.process(std::make_pair<uint64_t, uint64_t>((uint64_t)begin, (uint64_t)end));
  }

  auto coords = cell.get_coordinates();

  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto end_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX);
  auto end_position = *(reinterpret_cast<const int64_t*>(cell.get_field_ptr_for_query_idx(end_query_idx)));

  std::string contig_name;
  int64_t contig_position;
  if (!m_vid_mapper.get_contig_location(coords[1], contig_name, contig_position)) {
    std::cerr << "Could not find genomic interval associated with Variant(Call) at "
              << coords[1] << std::endl;
    return;
  }

  contig_position++;
  genomic_interval_t genomic_interval(std::move(contig_name),
                                      std::make_pair(contig_position, contig_position+end_position-coords[1]));

  std::vector<genomic_field_t> genomic_fields;
  // Ignore first field as it is "END"
  for (auto i=1u; i<query_config.get_num_queried_attributes(); i++) {
    if (cell.is_valid(i)) {
      genomic_field_t field_vec(query_config.get_query_attribute_name(i),
                                cell.get_field_ptr_for_query_idx(i),
                                cell.get_field_length(i));
      genomic_fields.push_back(field_vec);
    }
  }

  std::string sample_name;
  if (!m_vid_mapper.get_callset_name(coords[0], sample_name)) {
    sample_name = "NONE";
  }
   
  m_variant_call_processor.process(sample_name, coords, genomic_interval, genomic_fields);
}

#if(DEBUG)
void print_variant_calls(const VariantQueryConfig& query_config,
                         const VariantQueryProcessor& query_processor,
                         const VidMapper& vid_mapper) {
  std::string indent_prefix = "    ";
  std::cout << "{\n";
  //variant_calls is an array of dictionaries
  std::cout << indent_prefix << "\"variant_calls\": [\n";
  VariantCallPrintOperator printer(std::cout, indent_prefix+indent_prefix, &vid_mapper);
  query_processor.iterate_over_cells(query_processor.get_array_descriptor(), query_config, printer, true);
  std::cout << "\n" << indent_prefix << "]\n";
  std::cout << "}\n";
}
#endif

std::vector<VariantCall>* GenomicsDB::query_variant_calls(const std::string& array, VariantQueryConfig *query_config, GenomicsDBVariantCallProcessor& processor) {
#if(0)
  auto query_timer = Timer();
  query_timer.start();
#endif

  // Get a Variant Query Processor
  VariantQueryProcessor *query_processor = new VariantQueryProcessor(
      TO_VARIANT_STORAGE_MANAGER(m_storage_manager), array, query_config->get_vid_mapper());
  query_processor->do_query_bookkeeping(query_processor->get_array_schema(),
                                        *query_config, query_config->get_vid_mapper(), true);

#if(DEBUG)
  //print_variant_calls(*query_config, *query_processor, query_config->get_vid_mapper());
#endif

  // Perform Query over all the intervals
  std::vector<VariantCall> *pvariant_calls = new std::vector<VariantCall>;
  VidMapper* vid_mapper;
  if (m_vid_mapper)  {
    vid_mapper = TO_VID_MAPPER(m_vid_mapper);
  } else {
    vid_mapper = const_cast<VidMapper *>(&query_config->get_vid_mapper());
  }
  GatherVariantCalls gather_variant_calls(processor, *query_config);
  query_processor->iterate_over_cells(query_processor->get_array_descriptor(), *query_config, gather_variant_calls, true);

  delete query_processor;

#if(0)
  query_timer.stop();
  query_timer.print();
#endif

  return pvariant_calls;
}

void GenomicsDB::generate_vcf(const std::string& array,
                              genomicsdb_ranges_t column_ranges,
                              genomicsdb_ranges_t row_ranges,
                              const std::string& output,
                              const std::string& output_format,
                              bool overwrite) {
  // Create Variant Config for given concurrency_rank
  VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  VariantQueryConfig query_config(*base_query_config);
  query_config.set_array_name(array);
  if (column_ranges.size() > 0) {
    query_config.set_query_column_ranges(column_ranges);
  }
  if (row_ranges.size() > 0) {
    query_config.set_query_row_ranges(row_ranges);
  }

  query_config.validate();

  generate_vcf(array, &query_config, output, output_format, overwrite); 
}

void GenomicsDB::generate_vcf(const std::string& output, const std::string& output_format, bool overwrite) {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);
  generate_vcf(array, query_config, output, output_format, overwrite);
}

void GenomicsDB::generate_vcf(const std::string& array, VariantQueryConfig* query_config, const std::string& output, const std::string& output_format, bool overwrite) {
  if (output.length()) {
    query_config->set_vcf_output_filename(output);
  }
  if (output_format.length()) {
    query_config->set_vcf_output_format(output_format);
  }
  query_config->set_index_output_VCF(true);

  VERIFY(query_config->get_vcf_output_filename().length() && "VCF output filename not specified");
  if (!overwrite && TileDBUtils::is_file(query_config->get_vcf_output_filename())) {
      throw GenomicsDBException("VCF output file exists and overwrite set to false");
  }

  // Get a Variant Query Processor
  VariantQueryProcessor *query_processor = new VariantQueryProcessor(
      TO_VARIANT_STORAGE_MANAGER(m_storage_manager), array, query_config->get_vid_mapper());
  query_processor->do_query_bookkeeping(query_processor->get_array_schema(),
                                        *query_config, query_config->get_vid_mapper(), true);

  // Get a VCF Adapter
  VCFAdapter vcf_adapter;
  vcf_adapter.initialize(*query_config);
  auto *gvcf_operator = new BroadCombinedGVCFOperator(vcf_adapter, query_config->get_vid_mapper(), *query_config);

  for (auto i=0u; i<query_config->get_num_column_intervals(); i++) {
    query_processor->scan_and_operate(query_processor->get_array_descriptor(),
                                      *query_config, *gvcf_operator, i, true);
  }

  delete gvcf_operator;
  delete query_processor;
}

void GenomicsDB::generate_ped_map(const std::string& array,
                                  VariantQueryConfig* query_config,
                                  const std::string& output_prefix,
                                  double progress_interval,
                                  const std::string& fam_list) {

  std::cout << "FAM LIST " << fam_list << std::endl;

  GenomicsDBPedMapProcessor proc(query_config, progress_interval, fam_list); 

  query_variant_calls(array, query_config, (GenomicsDBVariantCallProcessor&)proc);
  proc.advance_state();
  query_variant_calls(array, query_config, (GenomicsDBVariantCallProcessor&)proc);
  proc.advance_state();
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
  if (vid_mapper->get_contig_location(variant_or_variant_call->get_column_begin(), contig_name,
                                      contig_position)) {
    contig_position++;
    return genomic_interval_t(std::move(contig_name),
                              std::make_pair(contig_position,
                                             contig_position
                                             +(variant_or_variant_call->get_column_end()
                                               -variant_or_variant_call->get_column_begin())));
  } else {
    // TODO: Extend api to return an genomic interval representing a NULL instead of throwing an exception
    throw GenomicsDBException("Could not find genomic interval associated with Variant(Call) at "
                             + std::to_string(variant_or_variant_call->get_column_begin()));
  }
}

genomic_interval_t GenomicsDB::get_genomic_interval(const genomicsdb_variant_t* variant) {
  return get_genomic_interval_for(TO_VARIANT(variant), get_vid_mapper(m_vid_mapper, m_query_config));
}

genomic_interval_t GenomicsDB::get_genomic_interval(const genomicsdb_variant_call_t* variant_call) {
  return get_genomic_interval_for(TO_VARIANT_CALL(variant_call), get_vid_mapper(m_vid_mapper, m_query_config));
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
    if (field && field->is_valid()) {
      std::string name = get_query_attribute_name(variant_or_variant_call, query_config, index);
      void* ptr = const_cast<void *>(field->get_raw_pointer());
      if (name.compare("ALT") == 0) {
        void *p;
        unsigned size;
        bool allocated;
        std::type_index type_idx= field->get_C_pointers(size, &p, allocated);
        if (allocated) {
          if (type_idx == std::type_index(typeid(int))) {
            ptr = *(reinterpret_cast<int **>(p));
          } else if (type_idx == std::type_index(typeid(float))) {
            ptr = *(reinterpret_cast<float **>(p));
          } else  if (type_idx == std::type_index(typeid(char))) {
            ptr = *(reinterpret_cast<char **>(p));
          }
        }
      }
      genomic_field_t field_vec(name,
                                ptr,
                                field->length());
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

GenomicsDBVariantCalls GenomicsDB::get_variant_calls(const std::string& array, const genomicsdb_variant_t* variant) {
  std::vector<VariantCall>* variant_calls = const_cast<std::vector<VariantCall>*>(&(TO_VARIANT(variant)->get_calls()));
  return GenomicsDBResults<genomicsdb_variant_call_t>(TO_GENOMICSDB_VARIANT_CALL_VECTOR(variant_calls),
                                                      create_genomic_field_types(*get_query_config_for(array)));
}

int64_t GenomicsDB::get_row(const genomicsdb_variant_call_t* variant_call) {
  return TO_VARIANT_CALL(variant_call)->get_row_idx();
}

// Navigate GenomicsDBResults

template<>
std::size_t GenomicsDBResults<genomicsdb_variant_t>::size() const noexcept {
  return TO_VARIANT_VECTOR(m_results)->size();
}

template<>
const genomicsdb_variant_t* GenomicsDBResults<genomicsdb_variant_t>::at(std::size_t pos) {
  if (pos >= size()) {
    return nullptr;
  } else {
    Variant *variant = TO_VARIANT_VECTOR(m_results)->data() + pos;
    return TO_GENOMICSDB_VARIANT(variant);
  }
}

template<>
void GenomicsDBResults<genomicsdb_variant_t>::free() {
  if (m_results != nullptr) {
    delete TO_VARIANT_VECTOR(m_results);
  }
  m_results = nullptr;
}

template<>
std::size_t GenomicsDBResults<genomicsdb_variant_call_t>::size() const noexcept {
  return TO_VARIANT_CALL_VECTOR(m_results)->size();
}

template<>
const genomicsdb_variant_call_t* GenomicsDBResults<genomicsdb_variant_call_t>::at(std::size_t pos) {
  if (pos >= size()) {
    return nullptr;
  } else {
    VariantCall *variant = TO_VARIANT_CALL_VECTOR(m_results)->data() + pos;
    return TO_GENOMICSDB_VARIANT_CALL(variant);
  }
}

template<>
void GenomicsDBResults<genomicsdb_variant_call_t>::free() {
  // NOP: free'ing is done at the Variant(genomicsdb_variant_t) level
}

void GenomicsDBVariantCallProcessor::process(const std::string& sample_name,
                                             const int64_t* coords,
                                             const genomic_interval_t& genomic_interval,
                                             const std::vector<genomic_field_t>& genomic_fields) {
  std::cout << "\t sample=" << sample_name << "\n";
  std::cout << "\t row=" << coords[0] << " position=" << coords[1]
            << "\n\t genomic_interval=" << genomic_interval.contig_name
            << ":" << genomic_interval.interval.first << "," << genomic_interval.interval.second << "\n";
  std::cout << "\t genomic_fields\n";
  for(auto genomic_field: genomic_fields) {
    std::cout << "\t\t" << genomic_field.name << ":" << genomic_field.to_string(get_genomic_field_type(genomic_field.name));
  }
  std::cout << std::endl;
}

void GenomicsDBVariantCallProcessor::process(const interval_t& interval) {
  std::cout << "----------------\nInterval:[" << interval.first << "," << interval.second << "]\n\n";
}

void GenomicsDBPedMapProcessor::process(const std::string& sample_name,
                                        const int64_t* coords,
                                        const genomic_interval_t& genomic_interval,
                                        const std::vector<genomic_field_t>& genomic_fields) {
  static size_t tm = 0;
  auto progress_bar = [&] () {
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
      
      logger.info("Query progress: r:{} c:{}, {} / {} = {:.2f}%", row, col, (long)col * total_rows + row, (long)total_rows * total_cols, 100 * ((double)col * total_rows + row)/(total_rows * total_cols));

      tm = now;
    }
  };

  if(progress_interval > 0) {
    progress_bar();
  }

  // ===============================================
  /*if(state) {
    std::cout << "\t sample=" << sample_name << "\n";
  std::cout << "\t row=" << coords[0] << " position=" << coords[1]
            << "\n\t genomic_interval=" << genomic_interval.contig_name
            << ":" << genomic_interval.interval.first << "," << genomic_interval.interval.second << "\n";
  std::cout << "\t genomic_fields\n";
  for(auto genomic_field: genomic_fields) {
    std::cout << "\t\t" << genomic_field.name << ":" << genomic_field.to_string(get_genomic_field_type(genomic_field.name));
  }
  std::cout << std::endl;
  }*/
  // ===============================================
  std::string ref_string, alt_string, gt_string, id_string;
  for(auto& f : genomic_fields) {
    if(f.name == "ALT") {
      std::string combined_alt = f.recombine_ALT_value();
      if(combined_alt.size()) {
        alt_string = combined_alt.substr(1, combined_alt.length() - 2);
      }
    }
    if(f.name == "REF") {
      ref_string = f.to_string(get_genomic_field_type(f.name));
    }
    if(f.name == "GT") {
      gt_string = f.to_string(get_genomic_field_type(f.name));
    }
    if(f.name == "ID") {
      id_string = f.to_string(get_genomic_field_type(f.name));
    }
  }

  if(!alt_string.size()) {
    logger.error("No ALT field for sample {}", sample_name);
    exit(1);
  }
  if(!ref_string.size()) {
    logger.error("No REF field for sample {}", sample_name);
    exit(1);
  }
  if(!gt_string.size()) {
    logger.error("No GT field for sample {}", sample_name);
    exit(1);
  }

  bool new_col = (last_coord != coords[1]);

  if(!state) {
    sample_map.insert(std::make_pair(sample_name, std::make_pair(-1, "")));
    variant_map.insert(std::make_pair(coords[1], std::make_pair(-1, ref_string)));
    if(new_col) {
      if(id_string.size()) {
        map_file << genomic_interval.contig_name << " " << id_string << " 0 " << genomic_interval.interval.first << std::endl;
      }
      else {
        map_file << genomic_interval.contig_name << " " << genomic_interval.contig_name << ":" << genomic_interval.interval.first << " 0 " << genomic_interval.interval.first << std::endl;
      }
    }
    last_coord = coords[1];
    return;
  }

  int sind = sample_map[sample_name].first;
  int vind = variant_map[coords[1]].first;

  if(state == 1) {
    // if skipped some samples backfill with ref
    /*int lo = sind + 1;
    if(vind > last_variant && last_variant != -1){
      lo = 0;
      for(int i = last_sample + 1; i < sample_map.size(); i++) {
        alt_file << variant_map[last_coord].second << " " << variant_map[last_coord].second << " ";
        alt_file_entries += 2;
      }
    }
    for(int i = lo; i < sind; i++) {
      alt_file << variant_map[coords[1]].second << " " << variant_map[coords[1]].second << " ";
      alt_file_entries += 2;
    }*/

    int to_fill = (vind - last_variant) * sample_map.size() + sind - (last_sample + 1);
    for(int i = 0; i < to_fill; i++) {
      alt_file << "0 0 ";
    }

    std::vector<std::string> vec = {ref_string};
    int index;
    while((index = alt_string.find(", ")) != std::string::npos) {
      vec.push_back(alt_string.substr(0, index));
      alt_string.erase(0, index + 2);
    }
    vec.push_back(alt_string);

    std::vector<int> gt_vec;
    while((index = gt_string.find("|")) != std::string::npos) {
      gt_vec.push_back(std::stoi(gt_string.substr(0, index)));
      gt_string.erase(0, index + 1);
    }
    gt_vec.push_back(std::stoi(gt_string));

    if(gt_vec.size() < 2) {
      logger.error("Samples must be at least diploid");
      exit(1);
    }

    alt_file << vec[gt_vec[0]] << " " << vec[gt_vec[1]] << " ";
    alt_file_entries += 2;
    last_sample = sind;
    last_variant = vind;
    last_coord = coords[1];
  }
}

void GenomicsDBPedMapProcessor::process(const interval_t& interval) {
  GenomicsDBVariantCallProcessor::process(interval);
}

void GenomicsDBPedMapProcessor::advance_state() {
  if(!state) {
    // associate samples with sorted position
    int i = -1;
    for(auto& a : sample_map) {
      ++i;
      a.second.first = (uint64_t)i;
    }
    // associate variants with sorted position
    i = -1;
    for(auto& b : variant_map) {
      ++i;
      b.second.first = (uint64_t)i;
    }
    last_sample = -1;
    last_variant = 0;
    last_coord = -1;
  
    // Find samples coincident with entries in fam files (if specified) and associate informatoin with sample name
    std::cout << "fam_list " << fam_list << std::endl;
    if(fam_list.length()) {
      std::ifstream fam_list_file(fam_list);
      std::string fname;
      std::cout << "reading " << fam_list << std::endl;
      while(std::getline(fam_list_file, fname)) {
        std::ifstream file(fname);
        std::string entry;
        std::cout << "fname is " << fname << std::endl;
        while(std::getline(file, entry)) {
          std::cout << "entry is " << entry << std::endl;
          std::string fid, wfid, fthid, mthid;
          char sex, pt;
          std::stringstream(entry) >> fid >> wfid >> fthid >> mthid >> sex >> pt;
          if(sample_map.count(wfid)) {
            sample_map[wfid].second = entry;
            std::cout << "sample_map[" << wfid << "] = " << sample_map[wfid].second << std::endl;
          }
        }
      }
    }
  }

  if(state == 1) {
    std::cout << "advance state " << std::endl;
    // if skipped some samples at end, fill with reample_map.insert(std::make_pair(sample_name, -1));
    for(int i = alt_file_entries; i < sample_map.size() * variant_map.size() * 2; i+=2) {
      //alt_file << variant_map[last_coord].second << " " << variant_map[last_coord].second << " ";
      alt_file << "0 0 ";
    }

    std::cout << "done with backfill" << std::endl;

    alt_file.clear();
    alt_file.seekg(0);

    std::cout << sample_map.size() << " samples " << std::endl;
    std::cout << variant_map.size() << " variants " << std::endl;

    int j = -1;
    for(auto& s : sample_map) {
      ++j;
      std::cout << "processing " << s.first << ", " << j << " / " << sample_map.size() << std::endl;
      // first 6 columns of ped: family id, within family id, within family id of father, mother, sex code, phenotype value
      if(s.second.second.length()) {
        ped_file << s.second.second;
      }
      else {
        ped_file << s.first << " " << s.first << " 0 0 0 0";
      }

      int index = s.second.first;
      for(int i = 0; i < variant_map.size(); i++) {
        if(!(i % 100)) {
          std::cout << "variant " << i << " / " << variant_map.size() << std::endl;
        }

        alt_file.clear();
        alt_file.seekg(0);

        std::string s1, s2;
        for(int j = 0; j <= index; j++){
          alt_file >> s1 >> s2;
        }
    
        ped_file << " " << s1 << " " << s2;
        index += sample_map.size();
      }
      ped_file << std::endl;
    }
  }
  state++;
}
