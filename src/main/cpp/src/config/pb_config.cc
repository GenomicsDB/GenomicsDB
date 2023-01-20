/**
 * The MIT License (MIT)
 * Copyright (c) 2019, 2023 Omics Data Automation, Inc
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

//Enable asserts
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <zlib.h>
#include <google/protobuf/util/json_util.h>
#include <google/protobuf/util/type_resolver.h>
#include <google/protobuf/util/type_resolver_util.h>

#include "tiledb_utils.h"
#include "json_config.h"
#include "genomicsdb_config_base.h"
#include "genomicsdb_logger.h"
#include "genomicsdb_status.h"

#include "genomicsdb_coordinates.pb.h"
#include "genomicsdb_import_config.pb.h"
#include "genomicsdb_export_config.pb.h"
#include "genomicsdb_vid_mapping.pb.h"
#include "genomicsdb_callsets_mapping.pb.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBConfigException(#X);

void GenomicsDBConfigBase::read_from_PB_binary_string(const std::string& str, const int rank) {
  genomicsdb_pb::ExportConfiguration export_config;
  auto status = export_config.ParseFromString(str);
  if(!status || !export_config.IsInitialized()) {
    throw GenomicsDBConfigException("Could not deserialize PB binary string for export");
  }
  read_from_PB(&export_config, rank);
}

int GenomicsDBImportConfig::read_from_PB_file(const std::string& filename, const int rank) {
  genomicsdb_pb::ImportConfiguration import_config;
  int status = get_pb_from_json_file(&import_config, filename);
  if (status == GENOMICSDB_OK) {
    read_from_PB(&import_config, rank);
  }
  return status;
}

ContigInfo GenomicsDBImportConfig::get_contig_info(const ContigPosition* contig_position) {
  ContigInfo contig_info;
  if (!m_vid_mapper.get_contig_info(contig_position->contig(), contig_info)) {
    logger.fatal(GenomicsDBConfigException(logger.format("Could not locate contig({}) position({}) in the vid mapper",
                                                         contig_position->contig(), contig_position->position())));
  }
  return std::move(contig_info);
}


ColumnRange GenomicsDBImportConfig::extract_contig_interval(const ContigPosition* start,
                                                            const ContigPosition* end) {
  ContigInfo contig_info = get_contig_info(start);
  ColumnRange column_range;
  if (end) {
    if (start->contig() != end->contig()) {
      logger.fatal(GenomicsDBConfigException(logger.format("Start contig({}) and End contig({}) in protobuf are different",
                                                           start->contig(), end->contig())));
    }
    column_range = verify_contig_position_and_get_tiledb_column_interval(contig_info,
                                                                         start->position(),
                                                                         end->position());
  } else {
    column_range = verify_contig_position_and_get_tiledb_column_interval(contig_info,
                                                                         start->position(),
                                                                         INT64_MAX-1);
  }
  return column_range;
}

void GenomicsDBImportConfig::read_from_PB(const genomicsdb_pb::ImportConfiguration* import_config, const int rank) {
  base_read_from_PB(import_config);

  // Get column partitions
  if (import_config->column_partitions_size() == 0) {
    logger.fatal(GenomicsDBConfigException("Loader PB should have at least one partition specified"));
  }
  m_column_partitions_specified = true;

  m_array_names.resize(import_config->column_partitions_size());
  m_sorted_column_partitions.resize(import_config->column_partitions_size());
  m_column_ranges.resize(import_config->column_partitions_size());
  
  auto partitions = import_config->column_partitions();
  m_single_workspace_path = true;
  m_single_array_name = partitions.size() == 1;
  m_column_partitions_specified = true;

  // Store the begins for every partition for lookup later
  std::unordered_map<int64_t, unsigned> begin_to_idx;
  auto partition_idx = 0u;
  for (auto partition: partitions) {
    if (partition.has_workspace()) {
      if (m_workspaces.size() == 0) {
          m_workspaces.emplace_back(partition.workspace());
      } else {
        m_workspaces.resize(import_config->column_partitions_size());
        m_workspaces[partition_idx] = partition.workspace();
        m_single_workspace_path = false;
      }
    }
    if (partition.generate_array_name_from_partition_bounds()) {
      logger.fatal(GenomicsDBConfigException("Generate array name from partition bounds not yet supported"));
    } else if (!partition.has_array_name()) {
      logger.fatal(GenomicsDBConfigException("Currently every partition should have an array name specified"));
    }
    m_array_names[partition_idx] = partition.array_name();
    
    m_column_ranges[partition_idx].resize(1);      //only 1 std::pair
    if (partition.begin().has_contig_position()) {
      auto* start = &partition.begin().contig_position();
      auto* end = partition.has_end()?&partition.end().contig_position():nullptr;
      m_column_ranges[partition_idx][0] = std::move(extract_contig_interval(start, end));
    } else {
      m_column_ranges[partition_idx][0].first = partition.begin().tiledb_column();
      if (partition.has_end()) {
        m_column_ranges[partition_idx][0].second = partition.end().tiledb_column();
      } else {
        m_column_ranges[partition_idx][0].second = INT64_MAX-1;
      }
    }  
    if (m_column_ranges[partition_idx][0].first > m_column_ranges[partition_idx][0].second) {
      std::swap<int64_t>(m_column_ranges[partition_idx][0].first, m_column_ranges[partition_idx][0].second);
    }
    m_sorted_column_partitions[partition_idx].first = m_column_ranges[partition_idx][0].first;
    m_sorted_column_partitions[partition_idx].second = m_column_ranges[partition_idx][0].second;
    begin_to_idx[m_column_ranges[partition_idx][0].first] = partition_idx;
    partition_idx++;
  }

  if ((m_single_workspace_path && m_workspaces.size() != 1) ||
      (!m_single_workspace_path && m_workspaces.size() != partitions.size())) {
    logger.fatal(GenomicsDBConfigException("List of workspaces should either be one for a single workspace or have workspaces specified for every partition"));
  }

  // Sort in ascending order and set end value if not set for the partition
  std::sort(m_sorted_column_partitions.begin(), m_sorted_column_partitions.end(), ColumnRangeCompare);
  for (auto i=0ull; i+1u<m_sorted_column_partitions.size(); ++i) {
    if (m_sorted_column_partitions[i].first == m_sorted_column_partitions[i+1u].first) {
      logger.fatal(GenomicsDBConfigException("Column partitions cannot have the same begin value"));
    }
    if (m_sorted_column_partitions[i].second >= m_sorted_column_partitions[i+1u].first) {
      m_sorted_column_partitions[i].second = m_sorted_column_partitions[i+1u].first-1;
    }
    auto idx = begin_to_idx[m_sorted_column_partitions[i].first];
    m_column_ranges[idx][0].second = m_sorted_column_partitions[i].second;
  }
  
  if (import_config->has_size_per_column_partition()) m_per_partition_size = import_config->size_per_column_partition();
  if (import_config->has_treat_deletions_as_intervals()) m_treat_deletions_as_intervals = import_config->treat_deletions_as_intervals();

  // genomicsdb options
  if (import_config->has_fail_if_updating()) m_fail_if_updating = import_config->fail_if_updating();
  if (import_config->has_disable_synced_writes()) m_disable_synced_writes = import_config->disable_synced_writes();
  if (import_config->has_ignore_cells_not_in_partition()) m_ignore_cells_not_in_partition = import_config->ignore_cells_not_in_partition();
    
  // vcf import options
  if (import_config->has_num_parallel_vcf_files()) m_num_parallel_vcf_files = import_config->num_parallel_vcf_files();
  if (import_config->has_do_ping_pong_buffering()) m_do_ping_pong_buffering = import_config->do_ping_pong_buffering();
  if (import_config->has_offload_vcf_output_processing()) m_offload_vcf_output_processing = import_config->offload_vcf_output_processing();

  // vcf generate options
  if (import_config->has_produce_combined_vcf()) {
    m_produce_combined_vcf = import_config->produce_combined_vcf();
    m_vcf_output_filename = "-"; // stdout, TODO: hardcoded for now, protobuf has it as part of partition
    if (import_config->has_reference_genome()) {
      m_reference_genome = import_config->reference_genome();
    } else {
      logger.fatal(GenomicsDBConfigException("Should specify a reference genome when produce_conbined_vcf is set to true"));
    }
    if (import_config->has_vcf_header_filename()) m_vcf_header_filename = import_config->vcf_header_filename();
    if (import_config->has_discard_vcf_index()) m_discard_vcf_index = import_config->discard_vcf_index();
  }

  // callset processing for sample rows to import
  if (import_config->has_lb_callset_row_idx()) m_lb_callset_row_idx = import_config->lb_callset_row_idx();
  if (import_config->has_ub_callset_row_idx()) m_ub_callset_row_idx = import_config->ub_callset_row_idx();
  fix_callset_row_idx_bounds(rank);

  // tiledb options
  if (import_config->has_produce_tiledb_array()) m_produce_tiledb_array = import_config->produce_tiledb_array();
  if (import_config->has_delete_and_create_tiledb_array()) m_delete_and_create_tiledb_array = import_config->delete_and_create_tiledb_array();
  if (import_config->has_segment_size()) m_segment_size = import_config->segment_size();
  if (import_config->has_num_cells_per_tile()) m_num_cells_per_tile = import_config->num_cells_per_tile();
  if (import_config->has_consolidate_tiledb_array_after_load()) m_consolidate_tiledb_array_after_load = import_config->consolidate_tiledb_array_after_load();

  // tiledb compression options
  if (import_config->has_compress_tiledb_array()) m_compress_tiledb_array = import_config->compress_tiledb_array();
  if (import_config->has_tiledb_compression_type()) m_tiledb_compression_type = import_config->tiledb_compression_type();
  if (import_config->has_tiledb_compression_level()) m_tiledb_compression_level = import_config->tiledb_compression_level();

  // other tiledb options
  if (import_config->has_enable_shared_posixfs_optimizations()) m_enable_shared_posixfs_optimizations = import_config->enable_shared_posixfs_optimizations();
  if (import_config->has_disable_delta_encode_for_offsets()) m_disable_delta_encode_offsets = import_config->disable_delta_encode_for_offsets();
  if (import_config->has_disable_delta_encode_for_coords()) m_disable_delta_encode_coords = import_config->disable_delta_encode_for_coords();
  if (import_config->has_enable_bit_shuffle_gt()) m_enable_bit_shuffle_gt = import_config->enable_bit_shuffle_gt();
  if (import_config->has_enable_lz4_compression_gt()) m_enable_lz4_compression_gt = import_config->enable_lz4_compression_gt();
}

void GenomicsDBConfigBase::read_from_PB(const genomicsdb_pb::ExportConfiguration* export_config, const int rank) {
  base_read_from_PB(export_config);
  
  m_workspaces.clear();
  m_workspaces.emplace_back(export_config->workspace());
  m_single_workspace_path = true;
  m_array_names.clear();
  if(!export_config->has_array_name())
    throw GenomicsDBConfigException("PB export config object must have array_name set");
  m_array_names.emplace_back(export_config->array_name());
  m_single_array_name = true;
  if(export_config->has_scan_full() && export_config->scan_full())
    scan_whole_array();
  m_attributes.resize(export_config->attributes_size());
  for(auto i=0;i<export_config->attributes_size();++i)
    m_attributes[i] = std::move(std::string(export_config->attributes(i)));
  if(export_config->has_query_filter())
    m_query_filter = export_config->query_filter();
  if(export_config->has_vcf_header_filename())
    m_vcf_header_filename = export_config->vcf_header_filename();
  m_vcf_output_filename = export_config->has_vcf_output_filename()
    ? export_config->vcf_output_filename() : "-"; //else stdout
  m_vcf_output_format = export_config->has_vcf_output_format()
    ? export_config->vcf_output_format() : "";
  if(!m_scan_whole_array
      && export_config->query_column_ranges_size() == 0
      && export_config->query_contig_intervals_size() == 0
      && export_config->query_row_ranges_size() == 0
      && export_config->query_sample_names_size() == 0)
    throw GenomicsDBConfigException("Must have one of \"scan_full\" or \"query_column_ranges\" or \
	\"query_contig_intervals\" or \"query_row_ranges\" \
	or \"query_sample_names\" in the Protobuf query config object");
  if(!m_scan_whole_array) {
    if(export_config->query_column_ranges_size() > 0 && export_config->query_contig_intervals_size() > 0)
      throw GenomicsDBConfigException("Protobuf query object cannot have both query_column_ranges and query_contig_intervals at the same time");
    if(export_config->query_column_ranges_size() > 0) {
      if(export_config->query_column_ranges_size() == 1)
	m_single_query_column_ranges_vector = true;
      m_column_ranges.resize(export_config->query_column_ranges_size());
      for(auto i=0;i<export_config->query_column_ranges_size();++i) {
	const auto& inner_list_wrapper = export_config->query_column_ranges(i);
	m_column_ranges[i].resize(inner_list_wrapper.column_or_interval_list_size());
	for(auto j=0;j<inner_list_wrapper.column_or_interval_list_size();++j) {
	  const auto& gdb_column_or_interval = inner_list_wrapper.column_or_interval_list(j);
	  if(gdb_column_or_interval.has_column()) {
	    const auto& gdb_column = gdb_column_or_interval.column();
	    if(gdb_column.has_tiledb_column()) {
	      m_column_ranges[i][j].first = gdb_column.tiledb_column();
	      m_column_ranges[i][j].second = gdb_column.tiledb_column();
	    }
	    else {
	      if(!gdb_column.has_contig_position())
		throw GenomicsDBConfigException("Protobuf object GenomicsDBColumn must have either tiledb_column \
		    or contig_position defined");
	      const auto& contig_position = gdb_column.contig_position();
	      assert(m_vid_mapper.is_initialized());
	      ContigInfo contig_info;
	      if (!m_vid_mapper.get_contig_info(contig_position.contig(), contig_info))
		throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : "
		    + contig_position.contig());
	      m_column_ranges[i][j] = verify_contig_position_and_get_tiledb_column_interval(contig_info,
		  contig_position.position(), contig_position.position());
	    }
	  }
	  else {
	    if(!gdb_column_or_interval.has_column_interval())
	      throw GenomicsDBConfigException("Protobuf object GenomicsDBColumnOrInterval must have either column \
		  or column_interval set");
	    const auto& gdb_column_interval = gdb_column_or_interval.column_interval();
	    if(gdb_column_interval.has_tiledb_column_interval()) {
	      m_column_ranges[i][j].first = gdb_column_interval.tiledb_column_interval().begin();
	      m_column_ranges[i][j].second = gdb_column_interval.tiledb_column_interval().end();
	    }
	    else {
	      if(!gdb_column_interval.has_contig_interval())
		throw GenomicsDBConfigException("GenomicsDBColumnInterval Protobuf object must have tiledb_column_interval\
		    or contig_interval");
	      const auto& contig_interval = gdb_column_interval.contig_interval();
	      assert(m_vid_mapper.is_initialized());
	      ContigInfo contig_info;
	      if (!m_vid_mapper.get_contig_info(contig_interval.contig(), contig_info))
		throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : "
		    + contig_interval.contig());
	      m_column_ranges[i][j] = verify_contig_position_and_get_tiledb_column_interval(contig_info,
		  contig_interval.begin(), contig_interval.end());
	    }
	  }
	  if(m_column_ranges[i][j].first > m_column_ranges[i][j].second)
	    std::swap(m_column_ranges[i][j].first, m_column_ranges[i][j].second);
	}
      }
    }
    if(export_config->query_contig_intervals_size() > 0) {
      m_single_query_column_ranges_vector = true;
      m_column_ranges.resize(1u);
      m_column_ranges[0].resize(export_config->query_contig_intervals_size());
      for(auto i=0;i<export_config->query_contig_intervals_size();++i) {
	const auto& contig_interval = export_config->query_contig_intervals(i);
	ContigInfo contig_info;
	if (!m_vid_mapper.get_contig_info(contig_interval.contig(), contig_info))
	  throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : "
	      + contig_interval.contig());
	if(!contig_interval.has_begin() && contig_interval.has_end())
	  throw GenomicsDBConfigException("ContigInterval must have begin if end is specified");
	auto begin = contig_interval.has_begin() ? contig_interval.begin() : 1ll;
	auto end = contig_interval.has_end() ? contig_interval.end()
	  : (contig_interval.has_begin() ? contig_interval.begin() : contig_info.m_length);
	m_column_ranges[0][i] = verify_contig_position_and_get_tiledb_column_interval(contig_info, begin, end);
      }
    }

    if(export_config->query_row_ranges_size() > 0 && export_config->query_sample_names_size() > 0)
      throw GenomicsDBConfigException("Protobuf query object cannot have both query_row_ranges and query_sample_names at the same time");
    if(export_config->query_sample_names_size() > 0) {
      m_single_query_row_ranges_vector = true;
      m_row_ranges.resize(1u);
      m_row_ranges[0].resize(export_config->query_sample_names_size());
      for(auto i=0;i<export_config->query_sample_names_size();++i) {
	int64_t row_idx = -1;
	auto status = m_vid_mapper.get_tiledb_row_idx(row_idx, export_config->query_sample_names(i));
	if(!status)
	  throw GenomicsDBConfigException(std::string("Unknown sample name ") + export_config->query_sample_names(i));
	m_row_ranges[0][i].first = row_idx;
	m_row_ranges[0][i].second = row_idx;
      }
    }
    if(export_config->query_row_ranges_size() > 0) {
      m_row_ranges.resize(export_config->query_row_ranges_size());
      m_single_query_row_ranges_vector = (m_row_ranges.size() == 1u);
      for(auto i=0;i<export_config->query_row_ranges_size();++i) {
	auto& row_range_vec = m_row_ranges[i];
	const auto& pb_row_range_list = export_config->query_row_ranges(i);
	row_range_vec.resize(pb_row_range_list.range_list_size());
	for(auto j=0;j<pb_row_range_list.range_list_size();++j) {
	  row_range_vec[j].first = pb_row_range_list.range_list(j).low();
	  row_range_vec[j].second = pb_row_range_list.range_list(j).high();
	  if(row_range_vec[j].second < row_range_vec[j].first)
	    std::swap(row_range_vec[j].first, row_range_vec[j].second);
	}
      }
    }
  }
  if(export_config->has_reference_genome())
    m_reference_genome = export_config->reference_genome();
  //Limit on max #alt alleles so that PL fields get re-computed
  m_max_diploid_alt_alleles_that_can_be_genotyped = export_config->has_max_diploid_alt_alleles_that_can_be_genotyped()
    ? export_config->max_diploid_alt_alleles_that_can_be_genotyped() : MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;
  m_max_genotype_count = export_config->has_max_genotype_count()
    ? export_config->max_genotype_count() : MAX_GENOTYPE_COUNT;
  m_combined_vcf_records_buffer_size_limit = export_config->has_combined_vcf_records_buffer_size_limit()
    ? export_config->combined_vcf_records_buffer_size_limit() : DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE;
  //Cannot be 0
  m_combined_vcf_records_buffer_size_limit = std::max<size_t>(1ull, m_combined_vcf_records_buffer_size_limit);
  //GATK CombineGVCF does not produce GT field by default - option to produce GT
  m_produce_GT_field = export_config->has_produce_gt_field() ? export_config->produce_gt_field() : false; 
  //GATK CombineGVCF does not produce FILTER field by default - option to produce FILTER
  m_produce_FILTER_field = export_config->has_produce_filter_field() ? export_config->produce_filter_field() : false; 
  //index output VCF file
  m_index_output_VCF = export_config->has_index_output_vcf() ? export_config->index_output_vcf() : false;
  //sites-only query - doesn't produce any of the FORMAT fields
  m_sites_only_query = export_config->has_sites_only_query() && export_config->sites_only_query();
  //when producing GT, use the min PL value GT for spanning deletions
  m_produce_GT_with_min_PL_value_for_spanning_deletions =
    export_config->has_produce_gt_with_min_pl_value_for_spanning_deletions()
    ? export_config->produce_gt_with_min_pl_value_for_spanning_deletions() : false;
  //Bypass the first intersecting interval phase and use just the simple traversal mode for a fast
  //fuzzy search if that is acceptable
  m_bypass_intersecting_intervals_phase =
      export_config->has_bypass_intersecting_intervals_phase() ? export_config->bypass_intersecting_intervals_phase() : false;
}

template<typename T>
void GenomicsDBConfigBase::base_read_from_PB(const T* pb_config) {
  VERIFY_OR_THROW(pb_config && pb_config->IsInitialized());
  read_and_initialize_vid_and_callset_mapping_if_available(pb_config);
  VERIFY_OR_THROW(m_vid_mapper.is_initialized() && m_vid_mapper.is_callset_mapping_initialized());
  
  if (pb_config->has_segment_size()) {
    m_segment_size = pb_config->segment_size();
  }
  
  //Enable Shared PosixFS Optimizations in TileDB
  if (pb_config->has_enable_shared_posixfs_optimizations()){
    m_enable_shared_posixfs_optimizations = pb_config->enable_shared_posixfs_optimizations();
  }
  
  m_workspaces.clear();
  m_single_workspace_path = true;
  m_array_names.clear();
}

template<typename T>
void GenomicsDBConfigBase::read_and_initialize_vid_and_callset_mapping_if_available(
    const T* pb_config) {
  if (pb_config->has_vid_mapping_file()) {
    VidMappingPB vidmap_pb;
    auto& vidmap_file = pb_config->vid_mapping_file();
    if (get_pb_from_json_file(&vidmap_pb, vidmap_file) == GENOMICSDB_OK) {
      m_vid_mapper.parse_vidmap_protobuf(&vidmap_pb);
    } else {
      logger.warn("Could not deserialize vid mapping file {} as protobuf. Trying to parse as a regular JSON file instead", vidmap_file);
      m_vid_mapping_file = vidmap_file;
      m_vid_mapper = std::move(FileBasedVidMapper(m_vid_mapping_file));
    }
  } else if(pb_config->has_vid_mapping()) {
    m_vid_mapper.parse_vidmap_protobuf(&(pb_config->vid_mapping()));
  }

  if (pb_config->has_callset_mapping_file()) {
    CallsetMappingPB callset_mapping_pb;
    auto& callset_mapping_file = pb_config->callset_mapping_file();
    if (get_pb_from_json_file(&callset_mapping_pb, callset_mapping_file) == GENOMICSDB_OK) {
      m_vid_mapper.parse_callset_protobuf(&callset_mapping_pb);
    } else {
      logger.warn("Could not deserialize callset mapping file {} as protobuf. Trying to parse as a regular JSON file instead", callset_mapping_file);
      m_callset_mapping_file = callset_mapping_file;
      m_vid_mapper.parse_callsets_json(callset_mapping_file, true);
    }
  } else if (pb_config->has_callset_mapping()) {
    m_vid_mapper.parse_callset_protobuf(&(pb_config->callset_mapping()));
  }

  if (!m_vid_mapper.is_callset_mapping_initialized()) {
    logger.fatal(GenomicsDBConfigException("Could not initialize callset mapping with either protobuf or regular JSON files"));
  }
  if (!m_vid_mapper.is_initialized()) {
    logger.fatal(GenomicsDBConfigException("Could not initialize vid mapping with either protobuf or regular JSON files"));
  }
}

int GenomicsDBConfigBase::get_pb_from_json_file(
        google::protobuf::Message *pb_config,
        const std::string& json_file) {
  char *buffer = 0;
  size_t buffer_length;
  if (TileDBUtils::read_entire_file(json_file, (void **)&buffer, &buffer_length) != TILEDB_OK || !buffer || buffer_length == 0) {
    free(buffer);
    logger.fatal(GenomicsDBConfigException("Could not read json file " + json_file));
  }
  google::protobuf::StringPiece pb_string(buffer, buffer_length);
  google::protobuf::util::JsonParseOptions json_options;
  json_options.ignore_unknown_fields = true;
  auto status = google::protobuf::util::JsonStringToMessage(pb_string, pb_config, json_options);
  free(buffer);
  if (!status.ok()) {
    logger.debug("json file {} does not seem to be serialized protobuf - error:{}", json_file, status.ToString());
    return GENOMICSDB_ERR;
  } else if (!pb_config->IsInitialized()) {
    logger.debug("Could not deserialize protobuf from json file {}", json_file);
    return GENOMICSDB_ERR;
  }
  return GENOMICSDB_OK;
}
