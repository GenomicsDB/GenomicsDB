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
    case PROTOBUF_BINARY_STRING: {
      query_config->read_from_PB_binary_string(query_configuration, concurrency_rank);

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

      // Create an annotationService class.
      m_annotation_service = new AnnotationService(query_configuration, contigs);
      AnnotationService* annotation_service = TO_ANNOTATION_SERVICE(m_annotation_service);
      break;
    }
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
                                       void *annotation_service) {
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
  if (column_ranges.size() > 0) {
    query_config.set_query_column_ranges(column_ranges);
  }
  if (row_ranges.size() > 0) {
    query_config.set_query_row_ranges(row_ranges);
  }

  query_config.validate();

  // ALT pointer will be to the first std::string in a std::vector owned by an instance of VariantFieldALTData
  auto field_types = create_genomic_field_types(query_config, m_annotation_service);

  auto it = field_types.find("ALT");
  if(it != field_types.end()) {
    it->second.type_idx = genomic_field_type_t::genomicsdb_basic_type::STRING;
  }

  return GenomicsDBVariants(TO_GENOMICSDB_VARIANT_VECTOR(query_variants(array, &query_config)), field_types);
}

GenomicsDBVariants GenomicsDB::query_variants() {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);

  // ALT pointer will be to the first std::string in a std::vector owned by an instance of VariantFieldALTData
  auto field_types = create_genomic_field_types(*query_config, m_annotation_service);

  auto it = field_types.find("ALT");
  if(it != field_types.end()) {
    it->second.type_idx = genomic_field_type_t::genomicsdb_basic_type::STRING;
  }

  return GenomicsDBVariants(TO_GENOMICSDB_VARIANT_VECTOR(query_variants(array, query_config)), field_types);
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
  query_config.set_query_column_ranges(column_ranges);
  query_config.set_query_row_ranges(row_ranges);

  query_config.validate();

  return GenomicsDBVariantCalls(TO_GENOMICSDB_VARIANT_CALL_VECTOR(query_variant_calls(array, &query_config, processor)),
                                create_genomic_field_types(query_config, m_annotation_service));
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls() {
  GenomicsDBVariantCallProcessor processor;
  return query_variant_calls(processor);
}

GenomicsDBVariantCalls GenomicsDB::query_variant_calls(GenomicsDBVariantCallProcessor& processor) {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  const std::string& array = query_config->get_array_name(m_concurrency_rank);
  return GenomicsDBVariantCalls(TO_GENOMICSDB_VARIANT_CALL_VECTOR(query_variant_calls(array, query_config, processor)),
                                create_genomic_field_types(*query_config, m_annotation_service));
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
    logger.warn("Could not find genomic interval associated with Variant(Call) at {}", coords[1]);
    return;
  }

  contig_position++;
  genomic_interval_t genomic_interval(std::move(contig_name),
                                      std::make_pair(contig_position, contig_position+end_position-coords[1]));

  std::vector<genomic_field_t> genomic_fields;
  std::string ref_value;
  std::string alt_value;
  // Ignore first field as it is "END"
  for (auto i=1u; i<query_config.get_num_queried_attributes(); i++) {
    if (cell.is_valid(i)) {
      genomic_field_t field_vec(query_config.get_query_attribute_name(i),
                                cell.get_field_ptr_for_query_idx(i),
                                cell.get_field_length(i));
      genomic_fields.push_back(field_vec);

      if(query_config.get_query_attribute_name(i) == "REF") {
        ref_value = field_vec.str_value();
      } else if(query_config.get_query_attribute_name(i) == "ALT" && field_vec.get_num_elements() >= 1) {
        alt_value = field_vec.str_value().substr(0, field_vec.str_value().find_first_of("|"));
      }
    }
  }

  if(!alt_value.empty() && m_annotation_service) {
    AnnotationService* annotation_service = TO_ANNOTATION_SERVICE(m_annotation_service);
    annotation_service->annotate(genomic_interval, ref_value, alt_value, genomic_fields);
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

template<class VariantOrVariantCall>
std::vector<genomic_field_t> get_genomic_fields_for(const std::string& array, const VariantOrVariantCall* variant_or_variant_call, VariantQueryConfig* query_config);

void GenomicsDB::generate_plink(const std::string& array,
                                VariantQueryConfig* query_config,
                                unsigned char format,
                                int compression,
                                bool one_pass,
                                bool verbose,
                                double progress_interval,
                                const std::string& output_prefix,
                                const std::string& fam_list) {
  if(!m_vid_mapper) {
    logger.error("No valid VidMapper, PLINK generation cancelled");
    return;
  }

  GenomicsDBPlinkProcessor proc(query_config, static_cast<VidMapper*>(m_vid_mapper), array, format, compression, verbose, progress_interval, output_prefix, fam_list, m_concurrency_rank);
  auto types = create_genomic_field_types(*query_config, m_annotation_service);
  auto it = types.find("ALT");
  if(it != types.end()) {
    it->second.type_idx = genomic_field_type_t::genomicsdb_basic_type::STRING;
  }

  proc.initialize(types);

  if(!one_pass) { // if one_pass is true, skip first pass and get participating samples from callset (will include samples without data)
    query_variants(array, query_config, proc);
  }
  proc.advance_state();
  query_variants(array, query_config, proc);
  proc.advance_state();

  VidMapper* vid_mapper = static_cast<VidMapper*>(m_vid_mapper);
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

  auto types = create_genomic_field_types(*get_query_config_for(array), m_annotation_service);
  auto it = types.find("ALT"); // REMOVE?
  if(it != types.end()) {
    it->second.type_idx = genomic_field_type_t::genomicsdb_basic_type::STRING;
  }

  return GenomicsDBResults<genomicsdb_variant_call_t>(TO_GENOMICSDB_VARIANT_CALL_VECTOR(variant_calls), types);
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

void GenomicsDBVariantProcessor::process(const std::vector<Variant>& variants) {
  for(const auto& variant : variants) {
    process(variant);
  }
}

GenomicsDBPlinkProcessor::GenomicsDBPlinkProcessor(VariantQueryConfig* qc,
                                                   VidMapper* vid_mapper,
                                                   const std::string& array,
                                                   unsigned char formats,
                                                   int compression,
                                                   bool verbose,
                                                   double progress_interval,
                                                   std::string prefix,
                                                   std::string fam_list,
                                                   int rank
                                                  ) : query_config(qc), vid_mapper(vid_mapper), array(array), compression(compression), verbose(verbose), progress_interval(progress_interval), prefix(prefix), fam_list(fam_list), rank(rank) {
  make_bgen = (bool)(formats & BGEN_MASK);
  make_bed = (bool)(formats & BED_MASK);
  make_tped = (bool)(formats & TPED_MASK);

  // For use with BGEN compression
  if(compression == 1) {
    TileDBUtils::create_codec(&codec, TILEDB_GZIP, Z_DEFAULT_COMPRESSION);
  }
  else {
    TileDBUtils::create_codec(&codec, TILEDB_ZSTD, 9);
  }

  // open various files
  if(make_tped) {
    tped_file.open(prefix + ".tped", std::ios::out);
  }
  if(make_bed) {
    bed_file.open(prefix + ".bed", std::ios::out | std::ios::binary);
    bim_file.open(prefix + ".bim", std::ios::out);
  }
  if(make_tped || make_bed) {
    fam_file.open(prefix + ".fam", std::ios::out);
  }
  if(make_bgen) {
    bgen_file.open(prefix + ".bgen", std::ios::out | std::ios::binary);
  }

  if(make_bed) {
    char magic_numbers[] = {0x6c, 0x1b, 0x01};
    bed_file.write(magic_numbers, 3); // BED: magic numbers
  }

  if(make_bgen) {
    int32_t zero = 0;
    int32_t offset = 20; // BGEN: offset of first variant data block relative to 5th byte. Always 20 here because free data area left empty
    bgen_file.write((char*)&offset, 4); // BGEN: first 4 bytes has offset
    bgen_file.write((char*)&offset, 4); // BGEN: beginning of header, next 4 bytes same in this case
    bgen_file.write((char*)&zero, 4); // BGEN: 4 bytes number of variants (M), filled in later
    bgen_file.write((char*)&zero, 4); // BGEN: 4 bytes number of samples (N), filled in later

    char bgen_magic_numbers[] = {'b', 'g', 'e', 'n'};
    bgen_file.write(bgen_magic_numbers, 4); // BGEN: 4 bytes bgen magic number

    int32_t flags = 0b10000000000000000000000000001000; // BGEN: layout 2, sample identifier block present
    flags = flags | compression;
    bgen_file.write((char*)&flags, 4); // BGEN: 4 bytes flags, end of header
  }

  // for use with progress bar
  for(auto& a : query_config->get_query_row_ranges(rank)) {
    total_rows += a.second - a.first + 1;
  }
  for(auto& b : query_config->get_query_column_ranges(rank)) {
    total_cols += b.second - b.first + 1;
  }
}

void GenomicsDBPlinkProcessor::bgen_variant_data_block(const std::string& rsid, const genomic_interval_t& genomic_interval, const std::vector<std::string>& vec, bool phased) {
  // BGEN: fill in genotype block size and min/max ploidy from previous iteration
  if(last_sample != -1) { // no need on first column
    if(samples_in_column < sample_map.size()) {
      min_ploidy = (min_ploidy > 2) ? 2 : min_ploidy;
      max_ploidy = (max_ploidy < 2) ? 2 : max_ploidy;
    }
    samples_in_column = 0;
    bgen_finish_gt();
  }
  min_ploidy = 64;
  max_ploidy = -1;

  // BGEN: variant data blocks
  int16_t zero = 0;
  int16_t one = 1;
  bgen_file.write((char*)&one, 2); // BGEN: 2 byte length of variant identifier, not stored in GenomicsDB so using dummy
  bgen_file.write((char*)&one, 1); // BGEN: dummy variant id
  int16_t rsid_len = rsid.length();
  bgen_file.write((char*)&rsid_len, 2); // BGEN: 2 byte length of rsid
  bgen_file.write(rsid.c_str(), rsid_len); // BGEN: rsid
  std::string chrom = genomic_interval.contig_name;
  int16_t chrom_len = chrom.length();
  bgen_file.write((char*)&chrom_len, 2); // BGEN: 2 byte chrom length
  bgen_file.write(chrom.c_str(), chrom_len); // BGEN: chrom
  uint32_t variant_pos = genomic_interval.interval.first;
  bgen_file.write((char*)&variant_pos, 4); // BGEN: 4 byte variant position
  int16_t K = vec.size();
  bgen_file.write((char*)&K, 2);// BGEN: 2 byte K for number of alleles
  // write alleles and lengths
  int32_t len;
  for(auto& a : vec) { // iterate over alleles
    len = a.length();
    bgen_file.write((char*)&len, 4); // BGEN: 4 byte length of allele
    bgen_file.write(a.c_str(), len); // BGEN: allele
  }

  // BGEN: genotype data block (layout 2, uncompressed)
  bgen_gt_size_offset = bgen_file.tellp();
  bgen_gt_size = 0;
  int32_t fourB_zero = 0;

  // BGEN: preallocate probability data storage
  int32_t N = sample_map.size();
  codec_buf.append((char*)&N, 4); // BGEN: 4 byte N
  codec_buf.append((char*)&K, 2); // BGEN: 2 byte K
  bgen_min_ploidy_offset = bgen_file.tellp();
  codec_buf.append((char*)&fourB_zero, 1); // BGEN: 1 byte min ploidy (placeholder)
  bgen_max_ploidy_offset = bgen_file.tellp();
  codec_buf.append((char*)&fourB_zero, 1); // BGEN: 1 byte max ploidy (placeholder)
  bgen_ploidy_info_offset = bgen_file.tellp();
  char default_sample_info = 0b10000010; // default missingness and ploidy information: set to missing/diploid unspecified
  for(int j = 0; j < N; j++) {
    codec_buf.append(&default_sample_info, 1); // BGEN: default sample information within this variant: because missing is set to 1, no need to backfill skipped cells
  }
  codec_buf.append((char*)&phased, 1); // BGEN: 1 byte phased
  char B = 8; // precision at one byte to avoid alignment difficulties
  codec_buf.append(&B, 1); // BGEN: 1 byte unsigned bits of precision
  bgen_probability_offset = bgen_file.tellp();
}

void GenomicsDBPlinkProcessor::process(const Variant& variant) {
  auto& calls = variant.get_calls();

  bool phased = true;
  std::string gt_string;

  for(auto& c : calls) {
    auto fields = get_genomic_fields_for(array, &c, query_config);
    
    for(auto& f : fields) {
      if(f.name == "GT") {
        gt_string = f.to_string(get_genomic_field_type(f.name));

        if(gt_string.find('/') != std::string::npos && gt_string.find('.') == std::string::npos) { // unphased if contains a slash or is ./. (might be able to use GL/PL, which are unphased probabilities)
          phased = false;
        }
        break;
      }
    }
  }

  for(auto& c : calls) {
    auto fields = get_genomic_fields_for(array, &c, query_config);

    int64_t coords[] = {(int64_t)c.get_row_idx(), (int64_t)c.get_column_begin()};
    int64_t end_position = c.get_column_end();

    std::string contig_name;
    int64_t contig_position;
    if (!vid_mapper->get_contig_location(coords[1], contig_name, contig_position)) {
      std::cerr << "Could not find genomic interval associated with Variant(Call) at "
        << coords[1] << std::endl;
      continue;
    }

    contig_position++;
    genomic_interval_t genomic_interval(std::move(contig_name),
                                        std::make_pair(contig_position, contig_position+end_position-coords[1]));

    std::string sample_name;
    if (!vid_mapper->get_callset_name(coords[0], sample_name)) {
      sample_name = "NONE";
    }

    process(sample_name, coords, genomic_interval, fields, phased);
  }
}

void GenomicsDBPlinkProcessor::process(const std::string& sample_name,
                                       const int64_t* coords,
                                       const genomic_interval_t& genomic_interval,
                                       const std::vector<genomic_field_t>& genomic_fields,
                                       const bool phased) {
  last_phased = phased;

  if(progress_interval > 0) {
    progress_bar(coords);
  }

  std::string ref_string, alt_string, gt_string, id_string, pl_string, gl_string, pq_string;
  for(auto& f : genomic_fields) {
    if(f.name == "ALT") {
      std::string combined_alt = f.recombine_ALT_value(get_genomic_field_type(f.name));
      if(combined_alt.size()) {
        alt_string = combined_alt.substr(1, combined_alt.length() - 2);
      }
    }
    else if(f.name == "REF") {
      ref_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "GT") {
      gt_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "ID") {
      id_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "PL") {
      pl_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "GL") {
      gl_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "PQ") {
      pq_string = f.to_string(get_genomic_field_type(f.name));
    }
  }

  if(!alt_string.size()) {
    logger.error("No ALT field for sample: {} row: {} column: {}", sample_name, coords[0], coords[1]);
    exit(1);
  }
  if(!ref_string.size()) {
    logger.error("No REF field for sample: {} row: {} column: {}", sample_name, coords[0], coords[1]);
    exit(1);
  }
  if(!gt_string.size()) {
    logger.error("No GT field for sample: {} row: {} column: {}", sample_name, coords[0], coords[1]);
    exit(1);
  }

  // collect alleles
  std::vector<std::string> vec = {ref_string};
  size_t index;
  while((index = alt_string.find(", ")) != std::string::npos) {
    vec.push_back(alt_string.substr(0, index));
    alt_string.erase(0, index + 2);
  }
  vec.push_back(alt_string);

  std::vector<int> gt_vec;
  auto iter = gt_string.begin();

  // parse GT
  try {
    while((iter = find_if(gt_string.begin(), gt_string.end(), [] (char c) { return c == '|' || c == '/'; })) != gt_string.end()) {
      index = iter - gt_string.begin();
      gt_vec.push_back(std::stoi(gt_string.substr(0, index)));
      gt_string.erase(0, index + 1);
    }
    gt_vec.push_back(std::stoi(gt_string));
  }
  catch (...) { // probably ./., treat as missing cell
    return;
  }

  for(auto g : gt_vec) {
    if(g < 0 || g >= vec.size()) {
      if(verbose) {
        logger.error("GT field for sample: {} row: {} column: {},  contains index {}, which is out of bounds (1 ref allele, {} alt allele(s))", sample_name, coords[0], coords[1], g, vec.size() - 1);
       }
      return;
    }
  }

  int16_t ploidy = gt_vec.size();

  if(!state) {
    sample_map_initialized = true; // possible to skip first pass, sample map will be populated from callset

    sample_map.insert(std::make_pair(coords[0], std::make_pair(-1, sample_name)));
    last_coord = coords[1];
    return;
  }

  if(ploidy != 2 && (make_bed || make_tped)) { 
    logger.error("The tped/bed formats do not support ploidies other than 2.");
    make_bed = 0;
    make_tped = 0;
    if(make_bgen) {
      logger.info("Continuing bgen generation");
    } 
  }

  if(state == 1) {
    std::string rsid;
    std::string rsid_row;
    if(id_string.size()) {
      rsid = id_string;
    }
    else {
      rsid = genomic_interval.contig_name + ":" + std::to_string(genomic_interval.interval.first);
    }
    rsid_row = genomic_interval.contig_name + " " + rsid + " 0 " + std::to_string(genomic_interval.interval.first);

    int sind = sample_map[coords[0]].first;

    // backfill if needed
    bool backfill = coords[1] - last_coord && last_coord != -1;
    int add_to_prev = backfill ? sample_map.size() - (last_sample + 1) : 0;
    int add_to_current = backfill ? sind : sind - (last_sample + 1);

    if(last_coord != coords[1]) { // first in variant
      num_variants++;

      for(int i = 0; i < add_to_prev; i++) { // backfill samples missing from previous variant
        if(make_tped) {
          tped_file << " 0 0";
        }
        if(make_bed) {
          write_to_bed(1);
        }
        if(make_bgen) {
          // BGEN: backfill probability data
          bgen_empty_cell(2, last_alleles, phased);
        }
      }
      if(make_bed) {
        flush_to_bed();
      }

      if(make_tped) {
        if(last_sample != -1) { // first line should not have newline
          tped_file << std::endl;
        }
        tped_file << rsid_row;
      }

      if(make_bgen) {
        bgen_variant_data_block(rsid, genomic_interval, vec, phased);
      }
    }
    else {
      samples_in_column++;
    }

    for(int i = 0; i < add_to_current; i++) { // backfill samples missing from current variant
      if(make_bed) {
        write_to_bed(1);
      }

      if(make_tped) {
        tped_file << " 0 0";
      }

      if(make_bgen) {
        // BGEN: backfill probability data
        bgen_empty_cell(2, vec.size(), phased);
      }
    }

    // safe to update now that backfilling is over
    last_alleles = vec.size();

    if(make_bed) {
      // 0 is homozygous for first allele
      // 1 is missing genotype
      // 2 is heterozygous
      // 3 is homozygous for second allele
      char gt_code;

      // FIXME cover with tests
      // FIXME will never get code 1 (call skipped if no GT) or 3 (no clear concept of a first vs second allele if they match)
      gt_code = char((gt_vec[0] != gt_vec[1]) * 2);

      if(last_coord != coords[1]) {
        bim_file << rsid_row;
        bim_file << " " << vec[gt_vec[0]] << " " << vec[gt_vec[1]] << std::endl;
      }
      write_to_bed(gt_code);
    }

    if(make_tped) {
      tped_file << " " << vec[gt_vec[0]] << " " << vec[gt_vec[1]];
    }

    last_sample = sample_map[coords[0]].first;
    //last_variant = vind;
    last_coord = coords[1];

    // convert PL to BGEN format
    std::vector<double> pl_vec;
    bool pl_dot = false;

    try {
      for(auto& tok : split(pl_string)) {
        if(tok == ".") pl_dot = true;
        int val = std::stoi(tok);
        if(val >= 0) {
          pl_vec.push_back(val);
        }
        else {
          throw std::runtime_error("PL value is negative");
        }
      }
    }
    catch(...) {
      pl_vec.clear();
      if(!pl_dot && verbose) {
        logger.error("PL field for sample: {} row: {} column: {} contains a negative value or otherwise did not parse, full string: {}", sample_name, coords[0], coords[1], pl_string);
      }
    }

    std::vector<char> pl_probs;

    double pl_total = 0;
    double epsilon = .05;
    for(auto& a : pl_vec) {
      double prob = std::pow(10, double(a)/-10);
      pl_total += prob;
      pl_probs.push_back(char(std::numeric_limits<unsigned char>::max() * prob));
    }

    // parse GL
    std::vector<double> gl_vec;
    bool gl_dot = false;

    try {
      for(auto& tok : split(gl_string)) {
        if(tok == ".") gl_dot = true;
        double val = std::stod(tok);
        if(val <= 0) {
          gl_vec.push_back(val);
        }
        else {
          throw std::runtime_error("GL value is positive");
        }
      }
    }
    catch(...) {
      gl_vec.clear();
      if(!gl_dot && verbose) {
        logger.error("GL field for sample: {} row: {} column: {} contains a strictly positive value or otherwise did not parse, full string {}", sample_name, coords[0], coords[1], gl_string);
      }
    }

    std::vector<char> gl_probs;

    double gl_total = 0;
    for(auto& a : gl_vec) {
      double prob = std::pow(10, a);
      gl_total += prob;
      gl_probs.push_back(char(std::numeric_limits<unsigned char>::max() * prob));
    }

    std::vector<char> probs;

    // subtract 1 representing reference
    auto num_genotypes = KnownFieldInfo::get_number_of_genotypes(vec.size() - 1, ploidy);

    if(gl_probs.size()) { // prefer gl as it is more precise (pl is rounded)
      if(gl_probs.size() == num_genotypes) {
        if(std::abs(1 - gl_total) > epsilon && verbose) {
          logger.warn("GL probabilities at sample: {} row: {} column: {} sum to {}, expected near 1. Generated BGEN may be invalid", sample_name, coords[0], coords[1], gl_total);
        }
        probs = gl_probs;
      }
      else {
        if(verbose) {
          logger.error("GL length at sample: {} row: {} column: {} is {}, expected {}. Defaulting to using GT to construct probabilities", sample_name, coords[0], coords[1], gl_vec.size(), num_genotypes);
        }
      }
    }
    else if(pl_probs.size()) {
      if(pl_probs.size() == num_genotypes) {
        if(std::abs(1 - pl_total) > epsilon && verbose) {
          logger.warn("PL probabilities at sample: {} row: {} column: {} sum to {}, expected near 1. Generated BGEN may be invalid", sample_name, coords[0], coords[1], pl_total);
        }
        probs = pl_probs;
      }
      else {
        if(verbose) {
          logger.error("PL length at sample: {} row: {} column: {} is {}, expected {}. Defaulting to using GT to construct probabilities", sample_name, coords[0], coords[1], pl_vec.size(), num_genotypes);
        }
      }
    }

    double pq;
    if(pq_string.length()) {
      try {
        pq = std::pow(10, (double)std::stoi(pq_string)/-10);
      }
      catch(...) {
        pq_string.clear();
      }
    }

    auto write_phased_probability = [&] (const std::vector<int>& v, size_t ind) {
      char p = gt_vec[v[0]] == v[1] ? -1 : 0;
      codec_buf.push_back(p);
    };

    auto write_unphased_probability = [&] (const std::vector<int>& v, size_t ind) {
      char p;

      if(!probs.size()) {
        std::vector<int> counts(vec.size(), 0);
        for(auto& g : gt_vec) {
          counts[g]++;
        }
        p = counts == v ? -1 : 0;
      }
      else {
        if(ind < probs.size()) {
          p = probs[ind];
        }
        else {
          // NOTE: this should never be triggered, GL/PL length is checked above
          logger.error("BGEN generation error: GL/PL probabilies only have {} term(s), halting BGEN generation", probs.size());
          make_bgen = 0;
        }
      }
      codec_buf.push_back(p);
    };

    if(make_bgen) {
      // store sample as not missing/ploidy info
      if(ploidy < 64) {
         min_ploidy = ploidy < min_ploidy ? ploidy : min_ploidy;
        max_ploidy = ploidy > max_ploidy ? ploidy : max_ploidy;
        char p = ploidy;
        codec_buf[8 + sample_map[coords[0]].first] = p;

        if(phased) { // phased
          bgen_enumerate_phased(ploidy, vec.size(), write_phased_probability);
        }
         else { // unphased
          bgen_enumerate_unphased(ploidy, vec.size(), write_unphased_probability);
        }
      }
    }
  }
}

void GenomicsDBPlinkProcessor::process(const interval_t& interval) {
  return;
}

void GenomicsDBPlinkProcessor::advance_state() {
  if(!state) {
    if(!sample_map_initialized) { // populate from callset
      std::string str;
      for(auto& row : query_config->get_rows_to_query()) {
        vid_mapper->get_callset_name(row, str);
        sample_map.insert(std::make_pair(row, std::make_pair(-1, str)));
      }
    }

    // associate samples with sorted position
    int i = -1;
    for(auto& a : sample_map) {
      ++i;
      a.second.first = (uint64_t)i;
    }
    // associate variants with sorted position
    i = -1;
    last_sample = -1;
    last_coord = -1;
  
    if(make_bed || make_tped) {
      // Find samples coincident with entries in fam files (if specified) and associate information with sample name
      std::map<std::string, std::string> fam_entries;
    
      if(fam_list.length()) {
        std::ifstream fam_list_file(fam_list);
        std::string fname;
        while(std::getline(fam_list_file, fname)) {
          std::ifstream file(fname);
          std::string entry;
          while(std::getline(file, entry)) {
            std::string fid, wfid, fthid, mthid;
            char sex, pt;
            std::stringstream(entry) >> fid >> wfid >> fthid >> mthid >> sex >> pt;
            fam_entries.insert(std::make_pair(wfid, entry));
          }
        }
      }
      for(auto& s : sample_map) {
        if(fam_entries.count(s.second.second)) {
          fam_file << fam_entries[s.second.second] << std::endl;
        }
        else {
          fam_file << s.second.second << " " << s.second.second << " 0 0 0 0" << std::endl;
        }
      }
    }

    if(make_bgen) {
      // BGEN: fill in N in header
      bgen_file.seekp(12);
      int32_t N = sample_map.size();
      bgen_file.write((char*)&N, 4); // BGEN: 4 byte N
      // Fill in M at end, as first pass, which counts variants, may be skipped

      // BGEN: write sample identifier block
      bgen_file.seekp(24);
      int32_t lsi = 0;
      bgen_file.write((char*)&lsi, 4); // BGEN: 4 byte total length of sample identifier block, filled in last
      bgen_file.write((char*)&N, 4); // BGEN: 4 byte N, deliberately duplicated here as per spec
      lsi = 8; // total length includes 8 bytes metainfo

      int16_t len;
      for(auto& s : sample_map) { // BGEN: write each sample id, can potentially be combined with above foreach loop
        len = s.second.second.length();
        lsi += len + 2; // total length increased by length of identifier and 2 byte length field
       
        bgen_file.write((char*)&len, 2);
        bgen_file.write(s.second.second.c_str(), len);
      }

      bgen_file.seekp(24);
      bgen_file.write((char*)&lsi, 4); // BGEN: 4 byte total length of sample identifier block, now with correct value
       bgen_file.seekp(0);
      int32_t offset = lsi + 20;
      bgen_file.write((char*)&offset, 4);
      // BGEN: update initial offset to include size of sample block
      bgen_file.seekp(24 + lsi); // seek to end of sample identifier block
    }
  }

  if(state == 1) {
    // if skipped some samples at end, fill with reample_map.insert(std::make_pair(sample_name, -1));
    int add_to_current = sample_map.size() - (last_sample + 1);

    for(int i = 0; i < add_to_current; i++) {
      if(make_tped) {
        tped_file << " 0 0";
      }
      if(make_bed) {
        write_to_bed(1);
      }
      if(make_bgen) {
        // BGEN: backfill last variant
        bgen_empty_cell(2, last_alleles, last_phased);
      }
    }
    if(make_bed) {
      flush_to_bed();
    }

    if(make_tped) {
      tped_file << std::endl;
    }

    if(sample_map.size() && make_bgen) {
      bgen_finish_gt();
    }

    if(make_bgen) {
      // BGEN: fill in M in header
      bgen_file.seekp(8);
      int32_t M = num_variants;
      bgen_file.write((char*)&M, 4); // BGEN: 4 byte M
    }
  }
  state++;
}
