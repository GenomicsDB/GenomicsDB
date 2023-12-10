/**
 * @file genomicsdb.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2019-2020,2022 Omics Data Automation, Inc.
 * Copyright (c) 2023 dātma, inc™
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

#include "annotation_service.h"
#include "broad_combined_gvcf.h"
#include "query_variants.h"
#include "tiledb_utils.h"
#include "timer.h"
#include "variant_field_data.h"
#include "variant_operations.h"
#include "variant_query_config.h"
#include "vcf_adapter.h"
#include "vid_mapper_pb.h"

#define TO_VARIANT_QUERY_CONFIG(X) (reinterpret_cast<VariantQueryConfig *>(static_cast<void *>(X)))
#define TO_VARIANT_STORAGE_MANAGER(X) (reinterpret_cast<VariantStorageManager *>(static_cast<void *>(X)))
#define TO_VID_MAPPER(X) (reinterpret_cast<VidMapper *>(static_cast<void *>(X)))
#define TO_ANNOTATION_SERVICE(X) (reinterpret_cast<AnnotationService *>(static_cast<void *>(X)))

#define TO_VARIANT(X) (reinterpret_cast<const Variant *>(static_cast<const void *>(X)))
#define TO_VARIANT_CALL(X) (reinterpret_cast<const VariantCall *>(static_cast<const void *>(X)))

#define TO_GENOMICSDB_VARIANT(X) (reinterpret_cast<const genomicsdb_variant_t *>(static_cast<const void *>(X)))
#define TO_GENOMICSDB_VARIANT_CALL(X) (reinterpret_cast<const genomicsdb_variant_call_t *>(static_cast<const void *>(X)))

#define TO_VARIANT_VECTOR(X) (reinterpret_cast<std::vector<Variant> *>(static_cast<void *>(X)))
#define TO_VARIANT_CALL_VECTOR(X) (reinterpret_cast<std::vector<VariantCall> *>(static_cast<void *>(X)))

#define TO_GENOMICSDB_VARIANT_VECTOR(X) (reinterpret_cast<std::vector<genomicsdb_variant_t> *>(static_cast<void *>(X)))
#define TO_GENOMICSDB_VARIANT_CALL_VECTOR(X) (reinterpret_cast<std::vector<genomicsdb_variant_call_t> *>(static_cast<void *>(X)))


#define VERIFY(X) if(!(X)) throw GenomicsDBException(#X);

void check(const std::string& workspace,
           const uint64_t segment_size) {
  VERIFY(!workspace.empty() && "Workspace specified cannot be empty");
  VERIFY(segment_size && "Segment size specified has to be greater than 0");
}

void check(const std::string& workspace,
           const uint64_t segment_size,
           const std::string& callset_mapping_file,
           const std::string& vid_mapping_file) {
  check(workspace, segment_size);
  VERIFY(!callset_mapping_file.empty() && "Callset Mapping File cannot be empty");
  VERIFY(!vid_mapping_file.empty() && "Vid Mapping File cannot be empty");
}


GenomicsDB::GenomicsDB(const std::string& workspace,
                       const std::string& callset_mapping_file,
                       const std::string& vid_mapping_file,
                       const std::vector<std::string>attributes,
                       const uint64_t segment_size) {
  check(workspace, segment_size, callset_mapping_file, vid_mapping_file);

  // Create base query configuration
  m_query_config = new VariantQueryConfig();

  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);

  query_config->set_workspace(workspace);
  query_config->set_callset_mapping_file(callset_mapping_file);
  query_config->set_vid_mapping_file(vid_mapping_file);
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
    case PROTOBUF_BINARY_STRING: {
      query_config->read_from_PB_binary_string(query_configuration, concurrency_rank);
      if (query_config->has_annotation_sources()) {
        // Determine which chromosomes are included in the query and pass the list to the annotation service
        // so it knows not to bother opening chromosome-specific vcfs that aren't relevant.
        // ColumnRange is std::pair<int64_t, int64_t> see genomicsdb.h:71
        ColumnRange column_partition = query_config->get_column_partition(concurrency_rank);
        // ContigIntervalTuple is std::tuple<std::string, int64_t, int64_t> see vid_mapper.cc:33
        std::vector<ContigIntervalTuple> contig_intervals = query_config->get_vid_mapper().get_contig_intervals_for_column_partition(column_partition.first, column_partition.second, true);
        std::set<std::string> contigs;
        for (auto interval : contig_intervals) {
          contigs.insert(std::get<0>(interval));
        }
        // Create an annotationService class
        m_annotation_service = new AnnotationService(query_configuration, contigs);
        AnnotationService* annotation_service = TO_ANNOTATION_SERVICE(m_annotation_service);
      }
      break;
    }
    default:
      throw GenomicsDBException("Unsupported query configuration type specified to the GenomicsDB constructor");
  }

  check(query_config->get_workspace(concurrency_rank),
        query_config->get_segment_size());

  // Discard intervals not part of this partition
  if (!loader_configuration_json_file.empty()) {
    query_config->subset_query_column_ranges_based_on_partition(loader_config, concurrency_rank);
  }

  // Create storage manager
  m_storage_manager = new VariantStorageManager(query_config->get_workspace(concurrency_rank), query_config->get_segment_size(), query_config->enable_shared_posixfs_optimizations());
}

GenomicsDB::~GenomicsDB() {
  // TODO: delete variant_query_config per array

  if (m_annotation_service != nullptr) {
  	delete TO_ANNOTATION_SERVICE(m_annotation_service);
  }
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

std::map<std::string, genomic_field_type_t> create_genomic_field_types(const VariantQueryConfig &query_config,
                                                                       void *annotation_service, bool change_alt_to_string=false) {
  std::map<std::string, genomic_field_type_t> genomic_field_types;
  for (auto i=0u; i<query_config.get_num_queried_attributes(); i++) {
    const std::string attribute_name = query_config.get_query_attribute_name(i);
    const FieldInfo* field_info = query_config.get_vid_mapper().get_field_info(attribute_name);
    if (field_info) {
      const FieldElementTypeDescriptor descr = field_info->get_vcf_type();
      assert(descr.get_num_elements_in_tuple() > 0);
      const std::type_index type_index = descr.get_tuple_element_type_index(0);
      const FieldLengthDescriptor length_descr = field_info->m_length_descriptor;
      if (change_alt_to_string && attribute_name == "ALT") {
        // ALT pointer will be to the first std::string in a std::vector owned by an instance of VariantFieldALTData
        genomic_field_types.insert(std::make_pair(
            attribute_name,
            genomic_field_type_t(genomic_field_type_t::genomicsdb_basic_type::STRING,
                                 length_descr.is_fixed_length_field(),
                                 length_descr.is_fixed_length_field()?length_descr.get_num_elements():0,
                                 length_descr.get_num_dimensions(),
                                 length_descr.is_length_ploidy_dependent()&&length_descr.contains_phase_information())));
      } else {
        genomic_field_types.insert(std::make_pair(
            attribute_name,
            genomic_field_type_t(type_index,
                                 length_descr.is_fixed_length_field(),
                                 length_descr.is_fixed_length_field()?length_descr.get_num_elements():0,
                                 length_descr.get_num_dimensions(),
                                 length_descr.is_length_ploidy_dependent()&&length_descr.contains_phase_information())));
      }
    }
  }
  if (annotation_service) {
    for(auto annotation_source: TO_ANNOTATION_SERVICE(annotation_service)->get_annotation_sources()) {
      auto field_types = annotation_source.field_types();
      genomic_field_types.insert(field_types.begin(), field_types.end());
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
  if (column_ranges.size() == 0) {
    query_config.set_query_column_ranges(SCAN_FULL);
  } else {
    query_config.set_query_column_ranges(column_ranges);
  }
  if (row_ranges.size() > 0) {
    query_config.set_query_row_ranges(row_ranges);
  }

  query_config.validate(m_concurrency_rank);

  return GenomicsDBVariants(TO_GENOMICSDB_VARIANT_VECTOR(query_variants(array, &query_config)),
                            create_genomic_field_types(query_config, m_annotation_service, true));
}

GenomicsDBVariants GenomicsDB::query_variants() {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);

  return GenomicsDBVariants(TO_GENOMICSDB_VARIANT_VECTOR(query_variants(array, query_config)),
                            create_genomic_field_types(*query_config, m_annotation_service, true));
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

void GenomicsDB::query_variants(const std::string& array, VariantQueryConfig* query_config, GenomicsDBVariantProcessor& proc) {
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
    proc.process(*pvariants);
    pvariants->clear();
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
  if (column_ranges.size() == 0) {
    query_config.set_query_column_ranges(SCAN_FULL);
  } else {
    query_config.set_query_column_ranges(column_ranges);
  }
  if (row_ranges.size() > 0) {
    query_config.set_query_row_ranges(row_ranges);
  }

  query_config.validate(m_concurrency_rank);

  return GenomicsDBVariantCalls(TO_GENOMICSDB_VARIANT_CALL_VECTOR(query_variant_calls(array, &query_config, processor)),
                                create_genomic_field_types(query_config, m_annotation_service));
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls() {
  GenomicsDBVariantCallProcessor processor;
  return query_variant_calls(processor, "", NONE);
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls(GenomicsDBVariantCallProcessor& processor,
                                                       const std::string& query_configuration,
                                                       const query_config_type_t query_configuration_type) {
  try {
    // Create Variant Config for given concurrency_rank
    VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
    VariantQueryConfig query_config(*base_query_config);

    if (query_configuration.size() > 0) {
      if (query_configuration_type != PROTOBUF_BINARY_STRING) {
        logger.fatal(GenomicsDBException(), "Unsupported query configuration type={} specified to query_variant_calls()",
                     query_configuration_type);
      }

      genomicsdb_pb::QueryConfiguration query_config_pb;
      auto status = query_config_pb.ParseFromString(query_configuration);
      if(!status || !query_config_pb.IsInitialized()) {
        logger.fatal(GenomicsDBException("Could not parse query_configuration. Only protobuf QueryConfiguration binary strings accepted as input argument"));
      }
      query_config.subset_query_from_PB(&query_config_pb, m_concurrency_rank);
    }

    if (!query_config.has_array_name(m_concurrency_rank)) {
      auto array_names = TileDBUtils::get_array_names(query_config.get_workspace(m_concurrency_rank));
      if (array_names.size() == 1) {
        query_config.set_array_name(array_names[0]);
      } else {
        // TODO: Should figure out the array names based on contig intervals as done on the JAVA side
        logger.fatal(GenomicsDBConfigException("Query configuration must either have array_name set or should be a single array in the workspace for now"));
      }
    }

    query_config.validate(m_concurrency_rank);
    
    return GenomicsDBVariantCalls(TO_GENOMICSDB_VARIANT_CALL_VECTOR(
        query_variant_calls(query_config.get_array_name(m_concurrency_rank), &query_config, processor)),
                                  create_genomic_field_types(query_config, m_annotation_service));
  } catch(const GenomicsDBIteratorException& e) {
    throw GenomicsDBException(e.what());
  }
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
  GatherVariantCalls(GenomicsDBVariantCallProcessor& variant_call_processor, const VariantQueryConfig& query_config, void* annotation_service) :
      SingleCellOperatorBase(), m_variant_call_processor(variant_call_processor), m_vid_mapper(query_config.get_vid_mapper()) {
    this->m_annotation_service = annotation_service;
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
  void* m_annotation_service = NULL;
};

void GatherVariantCalls::initialize(const VariantQueryConfig& query_config) {
  std::map<std::string, genomic_field_type_t> genomic_field_types = create_genomic_field_types(query_config, m_annotation_service);
  m_variant_call_processor.initialize(genomic_field_types);
};

void GatherVariantCalls::operate(VariantCall& variant_call,
                                 const VariantQueryConfig& query_config,
                                 const VariantArraySchema& schema) {
  // m_variant_calls.push_back(std::move(variant_call));
  logger.info("TBD: In GatherVariantCalls::operate()");
}

void GatherVariantCalls::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell,
                                                   const VariantQueryConfig& query_config,
                                                   const VariantArraySchema& schema) {
  if (cell.at_new_query_column_interval()) {
    auto curr_query_column_interval_idx = cell.get_current_query_column_interval_idx();
    auto begin = query_config.get_num_column_intervals()>0?query_config.get_column_begin(curr_query_column_interval_idx):0;
    auto end = query_config.get_num_column_intervals()>0?query_config.get_column_end(curr_query_column_interval_idx):INT64_MAX-1;
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
    logger.warn("Could not find genomic interval associated with Variant(Call) at {}", coords[1]);
    return;
  }

  contig_position++;
  genomic_interval_t genomic_interval(std::move(contig_name),
                                      std::make_pair(contig_position, contig_position+end_position-coords[1]));

  std::vector<genomic_field_t> genomic_fields;
  
  // Ignore first field as it is "END"
  for (auto i=1u; i<query_config.get_num_queried_attributes(); i++) {
    if (cell.is_valid(i)) {
      genomic_field_t field(query_config.get_query_attribute_name(i),
                            cell.get_field_ptr_for_query_idx(i),
                            cell.get_field_length(i));
      genomic_fields.push_back(field);
    }
  }

  if (m_annotation_service) {
    std::string alt_value;
    std::string ref_value;
    for (auto field: genomic_fields) {
      if (field.name == "REF") {
        ref_value = field.str_value();
      } else if (field.name == "ALT" && field.get_num_elements() >= 1) {
        alt_value = field.str_value().substr(0, field.str_value().find_first_of("|"));
      }
    }

    if (!alt_value.empty()) {
      AnnotationService* annotation_service = TO_ANNOTATION_SERVICE(m_annotation_service);
      annotation_service->annotate(genomic_interval, ref_value, alt_value, genomic_fields);
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

#if(0)
  print_variant_calls(*query_config, *query_processor, query_config->get_vid_mapper());
#endif

  // Perform Query over all the intervals
  std::vector<VariantCall> *pvariant_calls = new std::vector<VariantCall>;

  GatherVariantCalls gather_variant_calls(processor, *query_config, m_annotation_service);
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
                              const std::string& reference_genome,
                              const std::string& vcf_header,
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
  query_config.set_reference_genome(reference_genome);
  query_config.set_vcf_header_filename(vcf_header);

  query_config.validate(m_concurrency_rank);

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
      genomic_field_t field_vec(name,
                                ptr,
                                field->length());
      fields.push_back(field_vec);
    }
    index++;
  }
  return fields;
}

template
std::vector<genomic_field_t> get_genomic_fields_for<Variant>(const std::string& array, const Variant* variant_or_variant_call, VariantQueryConfig* query_config);
template
std::vector<genomic_field_t> get_genomic_fields_for<VariantCall>(const std::string& array, const VariantCall* variant_or_variant_call, VariantQueryConfig* query_config);

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
                                                      create_genomic_field_types(*get_query_config_for(array), m_annotation_service, true));
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

#define PIPED_SEP '|'
#define SLASHED_SEP '/'
std::string_view get_segment(std::string_view str, int segment_pos) {
  std::string::size_type pos;
  std::string::size_type next_pos = -1;
  for (int j=0; j<segment_pos; j++) {
    pos = next_pos + 1;
    next_pos = str.find(PIPED_SEP);
    if (next_pos == std::string::npos) {
      if (j != segment_pos-1) {
        return ".";
      } else {
        return str.substr(pos, str.length());
      }
    }
  }
  return str.substr(pos, next_pos);
}

std::string resolve(std::vector<int32_t>& gt, std::string_view ref, std::string_view alt) {
  std::string resolved_gt;
  for (auto i=0ul; i<gt.size(); i++) {
    if (i%2) { // Phase
      if (gt[i]) resolved_gt += PIPED_SEP; else resolved_gt += SLASHED_SEP;
    } else if (gt[i] > 0) { // ALT
      resolved_gt += get_segment(alt, gt[i]);
    } else if (gt[i] == 0) { // REF
      resolved_gt += ref;
    } else { // UNKNOWN
      resolved_gt += ".";
    }
  }
  return resolved_gt;
}

const std::string resolve_gt(const std::vector<genomic_field_t>& genomic_fields)  {
  std::string ref;
  std::string alt;
  std::vector<int32_t> gt_vec;
  for (auto genomic_field: genomic_fields) {
    if (genomic_field.name == "REF") {
      ref = genomic_field.str_value();
    } else if (genomic_field.name == "ALT") {
      alt = genomic_field.str_value();
    } else if (genomic_field.name == "GT") {
      for (auto i=0u; i<genomic_field.num_elements; i++) {
        gt_vec.push_back(genomic_field.int_value_at(i));
      }
    }
  }
  return resolve(gt_vec, ref, alt);
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

void GenomicsDBVariantProcessor::process(const std::vector<Variant>& variants) {
  for(const auto& variant : variants) {
    process(variant);
  }
}
