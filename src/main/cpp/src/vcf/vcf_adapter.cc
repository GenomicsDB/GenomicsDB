/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 * Copyright (c) 2020-2021 Omics Data Automation, Inc.
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

#ifdef HTSDIR

#include "vcf_adapter.h"

#include "genomicsdb_logger.h"
#include "hfile_genomicsdb.h"
#include "htslib/tbx.h"
#include "tiledb_utils.h"
#include "vid_mapper.h"

//ReferenceGenomeInfo functions
void ReferenceGenomeInfo::initialize(const std::string& reference_genome) {
  if (!reference_genome.empty()) {
    if (!TileDBUtils::is_file(reference_genome)) {
      throw VCFAdapterException(logger.format("Specified reference genome {} not found", reference_genome));
    }
    m_reference_faidx = fai_load(reference_genome.c_str());
    if (!m_reference_faidx) {
      throw VCFAdapterException(logger.format("Issues with opening reference genome {}", reference_genome));
    }
    m_reference_last_seq_read = "";
    //buffer
    m_buffer.resize(32768u+8);     //32KB
  }
}

char ReferenceGenomeInfo::get_reference_base_at_position(const char* contig, const int64_t pos) {
  //See if pos is within the last buffer read
  if (strcmp(m_reference_last_seq_read.c_str(), contig) == 0 && m_reference_last_read_pos <= pos) {
    int offset = pos - m_reference_last_read_pos;
    if (offset < m_reference_num_bases_read)
      return m_buffer[offset];
  }
  hts_pos_t length = 0;
  faidx_fetch_seq_into_buffer(m_reference_faidx, contig, pos, pos+m_buffer.size()-8u, &(m_buffer[0]), &length);
  assert(length > 0);
  m_reference_last_seq_read = contig;
  m_reference_last_read_pos = pos;
  m_reference_num_bases_read = length;
  return m_buffer[0];
}

//VCFAdapter functions
bool VCFAdapter::add_field_to_hdr_if_missing(bcf_hdr_t* hdr, const VidMapper* id_mapper, const std::string& field_name, int field_type_idx) {
  auto field_info_ptr = id_mapper->get_field_info(field_name);
  auto is_multid_vector_or_tuple_element_field = (
        (field_info_ptr != 0)
        && ((field_info_ptr->get_genomicsdb_type().get_num_elements_in_tuple() > 1u)
            || (field_info_ptr->m_length_descriptor.get_num_dimensions() > 1u)));
  //VCF header structure in htslib is a mess
  auto field_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, field_name.c_str());
  //field name might exist in the hdr hash, but hrec might be null
  auto field_exists_in_vcf_hdr = (field_idx >= 0)
                                 && (bcf_hdr_idinfo_exists(hdr, field_type_idx, field_idx))
                                 && ((bcf_hdr_id2hrec(hdr, BCF_DT_ID, field_type_idx, field_idx)) != 0);
  //For multi-d vector fields and fields whose elements are tuples, sometimes the header
  //information is completely wrong. It should be of type String.
  //Do this by removing and adding the field again
  //Correctness depends on the newly added hrec getting the exact same idx as
  //the deleted one. In debug mode, verify this
  auto old_field_idx_before_deletion = -1;
  std::string description_value;
  if (is_multid_vector_or_tuple_element_field && field_exists_in_vcf_hdr) {
    //(hdr->id[BCF_DT_ID][field_idx].val->info[field_type_idx])
    //= field_type_idx      //BCF_HL_*
    //| (BCF_HT_STR << 4)   //String
    //| (BCF_VL_FIXED << 8) //Fixed (since Number=1)
    //| (1 << 12)
    //;
    auto hrec = (bcf_hdr_id2hrec(hdr, BCF_DT_ID, field_type_idx, field_idx));
    auto description_idx = bcf_hrec_find_key(hrec, "Description");
    if (description_idx >= 0)
      description_value = hrec->vals[description_idx];
    bcf_hdr_remove(hdr, field_type_idx, field_name.c_str());
    if (bcf_hdr_sync(hdr)) {
      logger.fatal(VCFAdapterException("Possible realloc() failure from bcf_hdr_sync() while adding missing field to hdr"));
    }
    field_exists_in_vcf_hdr = false;
    old_field_idx_before_deletion = field_idx;
  }
  //Field not found
  if (!field_exists_in_vcf_hdr) {
    std::string header_line = "##";
    switch (field_type_idx) {
    case BCF_HL_FLT:
      header_line += "FILTER";
      break;
    case BCF_HL_INFO:
      header_line += "INFO";
      break;
    case BCF_HL_FMT:
      header_line += "FORMAT";
      break;
    default:
      throw VCFAdapterException(std::string("Unknown field type ")+std::to_string(field_type_idx));
      break;
    }
    header_line += "=<ID="+field_name;
    if (field_type_idx != BCF_HL_FLT) {
      //GT is weird
      if (field_type_idx == BCF_HL_FMT && field_name == "GT")
        header_line += ",Number=1,Type=String,Description=\"Genotype\"";
      else {
        assert(field_info_ptr);
        header_line += ",Number=";
        if (field_info_ptr->get_vcf_type().get_tuple_element_bcf_ht_type(0u) == BCF_HT_FLAG)
          header_line += '0';
        else {
          //Must be string with Number=1 in the hdr
          if (is_multid_vector_or_tuple_element_field) {
            assert(field_info_ptr->get_vcf_type().get_tuple_element_bcf_ht_type(0u) == BCF_HT_STR);
            header_line += "1";
          } else {
            auto& length_descriptor = field_info_ptr->m_length_descriptor;
            switch (length_descriptor.get_length_descriptor(0u)) {
            case BCF_VL_FIXED:
              header_line += std::to_string(length_descriptor.get_num_elements());
              break;
            case BCF_VL_VAR:
              header_line += ".";
              break;
            case BCF_VL_A:
              header_line += "A";
              break;
            case BCF_VL_R:
              header_line += "R";
              break;
            case BCF_VL_G:
              header_line += "G";
              break;
            default:
              throw VCFAdapterException("Unhandled field length descriptor "
                                        +std::to_string(length_descriptor.get_length_descriptor(0u)));
              break;
            }
          }
        }
        header_line += ",Type=";
        if (is_multid_vector_or_tuple_element_field)
          header_line+="String";
        else {
          switch (field_info_ptr->get_vcf_type().get_tuple_element_bcf_ht_type(0u)) {
          case BCF_HT_FLAG:
            header_line += "Flag";
            break;
          case BCF_HT_INT:
          case BCF_HT_INT64:
            header_line += "Integer";
            break;
          case BCF_HT_REAL:
            header_line += "Float";
            break;
          case BCF_HT_CHAR:
          case BCF_HT_STR:
            header_line += "String";
            break;
          default:
            throw VCFAdapterException("Field type "
                                      +std::to_string(field_info_ptr->get_vcf_type().get_tuple_element_bcf_ht_type(0u))
                                      +" not handled");
            break;
          }
        }
      }
    }
    if (old_field_idx_before_deletion >=0 && !description_value.empty())
      header_line += ",Description="+description_value; //looks like description_value has " internally
    else
      header_line += ",Description=\""+field_name+"\"";
    header_line += ">";
    int line_length = 0;
    auto hrec = bcf_hdr_parse_line(hdr, header_line.c_str(), &line_length);
    bcf_hdr_add_hrec(hdr, hrec);
    if (bcf_hdr_sync(hdr)) {
      logger.fatal(VCFAdapterException("Possible realloc() failure from bcf_hdr_sync() while adding missing field to hdr"));
    }
#ifdef DEBUG
    if (old_field_idx_before_deletion >= 0)
      assert(bcf_hdr_id2int(hdr, BCF_DT_ID, field_name.c_str()) == old_field_idx_before_deletion);
#endif
    return true;
  } else {
    const auto* field_info_ptr = id_mapper->get_field_info(field_name);
    assert(field_info_ptr);
    auto field_ht_type = bcf_hdr_id2type(hdr, field_type_idx, field_idx);
    //Don't bother doing any checks for the GT field
    if (field_name != "GT" && field_ht_type != BCF_HT_STR && field_type_idx != BCF_HL_FLT) {
      //Allowed configurations - both the JSON and the header specify that:
      //The field is fixed length and agree on the length OR
      //The field is variable length OR
      //field type is BCF_HT_FLAG and VCF header says length is 0 (VCF spec), vid JSON says that length is 1 OR
      //field is a string and VCF header says length is 0
      if (!((field_info_ptr->m_length_descriptor.is_fixed_length_field()
             && bcf_hdr_id2length(hdr, field_type_idx, field_idx) == BCF_VL_FIXED
             && field_info_ptr->m_length_descriptor.get_num_elements()
             == static_cast<size_t>(bcf_hdr_id2number(hdr, field_type_idx, field_idx)))
            ||
            (!field_info_ptr->m_length_descriptor.is_fixed_length_field()
             && bcf_hdr_id2length(hdr, field_type_idx, field_idx) != BCF_VL_FIXED
            )
            ||
            (field_ht_type == BCF_HT_FLAG
             && field_info_ptr->m_length_descriptor.is_fixed_length_field()
             && field_info_ptr->m_length_descriptor.get_num_elements() == 1
             && bcf_hdr_id2length(hdr, field_type_idx, field_idx) == BCF_VL_FIXED
             && static_cast<int>(bcf_hdr_id2number(hdr, field_type_idx, field_idx)) == 0
            )
            ||
            ((field_ht_type == BCF_HT_CHAR || field_ht_type == BCF_HT_STR)
             && !(field_info_ptr->m_length_descriptor.is_fixed_length_field())
             && field_ht_type == BCF_HT_STR
             && bcf_hdr_id2length(hdr, field_type_idx, field_idx) == BCF_VL_FIXED
             && bcf_hdr_id2number(hdr, field_type_idx, field_idx) == 1)
           )
         )
        throw VCFAdapterException(std::string("Conflicting field length descriptors and/or field lengths in the vid JSON and VCF header for field ")+field_name);
      //Check for compatible field types
      auto compatible_types = std::unordered_map<int, std::unordered_set<int>> {
        { BCF_HT_FLAG, { BCF_HT_FLAG } },
        { BCF_HT_INT, { BCF_HT_INT, BCF_HT_INT64 } },
        { BCF_HT_REAL, { BCF_HT_REAL } },
        { BCF_HT_INT64, { BCF_HT_INT64 } },
        { BCF_HT_VOID, { BCF_HT_VOID } },
        { BCF_HT_CHAR, { BCF_HT_CHAR, BCF_HT_STR } },
        { BCF_HT_STR, { BCF_HT_CHAR, BCF_HT_STR } }
      };
      if (
        compatible_types.find(field_ht_type) != compatible_types.end()
        && compatible_types[field_ht_type].find(field_info_ptr->get_vcf_type().get_tuple_element_bcf_ht_type(0u))
        == compatible_types[field_ht_type].end()
      )
        throw VCFAdapterException(std::string("Conflicting data types in the vid JSON and VCF header for field ")+field_name);
    }
  }
  return false;
}

VCFAdapter::VCFAdapter(bool open_output) {
  m_open_output = open_output;
  clear();
  m_template_vcf_hdr = 0;
  m_output_fptr = 0;
  m_is_bcf = true;
  m_config_base_ptr = 0;
  m_output_VCF_index_type = VCFIndexType::VCF_INDEX_NONE;

  // Register plugin for URI support for htslib via GenomicsDB
  genomicsdb_htslib_plugin_initialize();
}

VCFAdapter::~VCFAdapter() {
  clear();
  if (m_template_vcf_hdr)
    bcf_hdr_destroy(m_template_vcf_hdr);
  if (m_open_output && m_output_fptr) {
    bcf_close(m_output_fptr);
    auto status = 0;
    switch (m_output_VCF_index_type) {
    case VCFIndexType::VCF_INDEX_CSI:
      status = bcf_index_build(m_config_base_ptr->get_vcf_output_filename().c_str(), 14); //bcftools had default 14
      break;
    case VCFIndexType::VCF_INDEX_TBI:
      status = tbx_index_build(m_config_base_ptr->get_vcf_output_filename().c_str(), 0, &tbx_conf_vcf);
      break;
    default:
      break; //do nothing
    }
    if (status != 0) {
      logger.warn("error in creating index for output file {}", m_config_base_ptr->get_vcf_output_filename());
    }
  }
  m_output_fptr = 0;
#ifdef DO_PROFILING
  m_vcf_serialization_timer.print("bcf_t serialization", std::cerr);
#endif
}

void VCFAdapter::clear() {
  m_reference_genome_info.clear();
}

void VCFAdapter::initialize(const GenomicsDBConfigBase& config_base) {
  //Read template header with fields and contigs
  if (config_base.get_vcf_header_filename().empty())
    m_template_vcf_hdr = initialize_default_header();
  else {
    std::string vcf_header_filename = config_base.get_vcf_header_filename();
    auto created_header_file_copy = false;
    // Move  contents to a temporary local file for non-local URIs
    if (TileDBUtils::is_cloud_path(vcf_header_filename)) {
      char tmp_filename[PATH_MAX];
      auto status = TileDBUtils::create_temp_filename(tmp_filename, PATH_MAX);
      if (status != 0)
        throw VCFAdapterException(std::string("Could not create temp filename for VCF template header"));
      status = TileDBUtils::move_across_filesystems(vcf_header_filename, tmp_filename);
      if (status != TILEDB_OK)
        throw VCFAdapterException(std::string("Could not copy contents of VCF header filename ")
                                  +vcf_header_filename+" to temporary file "+tmp_filename);
      vcf_header_filename.assign(tmp_filename);
      created_header_file_copy = true;
    }
    auto* fptr = bcf_open(vcf_header_filename.c_str(), "r");
    if (fptr == 0)
      throw VCFAdapterException(std::string("Could not open template VCF header file ")
                                +vcf_header_filename);
    m_template_vcf_hdr = bcf_hdr_read_required_sample_line(fptr, 0); //sample line is not required
    bcf_close(fptr);
    if (created_header_file_copy)
      unlink(vcf_header_filename.c_str());
  }
  //Output fptr
  auto& output_format = config_base.get_vcf_output_format();
  assert(GenomicsDBConfigBase::m_vcf_output_format_to_is_bcf_flag.find(output_format)
         != GenomicsDBConfigBase::m_vcf_output_format_to_is_bcf_flag.end());
  m_is_bcf = GenomicsDBConfigBase::m_vcf_output_format_to_is_bcf_flag[output_format];
  auto& output_filename = config_base.get_vcf_output_filename();
  m_output_VCF_index_type = VCF_INDEX_NONE;
  if (m_open_output) {
    m_output_fptr = bcf_open(output_filename.c_str(), ("w"+output_format).c_str());
    if (m_output_fptr == 0) {
      logger.fatal(VCFAdapterException(logger.format("Cannot write to output file {}", output_filename)));
    }
    if (config_base.index_output_VCF() && !GenomicsDBConfigBase::output_to_stdout(output_filename)) {
      if (output_format == "z")
        m_output_VCF_index_type = VCFIndexType::VCF_INDEX_TBI;
      else {
        if (output_format == "b")
          m_output_VCF_index_type = VCFIndexType::VCF_INDEX_CSI;
      }
    }
  }
  //Reference genome
  m_reference_genome_info.initialize(config_base.get_reference_genome());
  m_config_base_ptr = &config_base;
}

bcf_hdr_t* VCFAdapter::initialize_default_header() {
  auto hdr = bcf_hdr_init("w");
  bcf_hdr_append(hdr, "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">");
  bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">");
  if (bcf_hdr_sync(hdr)) {
    logger.fatal(VCFAdapterException("Posssible realloc() failure from bcf_hdr_sync() while initializing default header"));
  }
  return hdr;
}

void VCFAdapter::print_header() {
  if (bcf_hdr_write(m_output_fptr, m_template_vcf_hdr)) {
    logger.fatal(VCFAdapterException("bcf_hdr_write() failed while printing header"));
  }
}

void VCFAdapter::handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size) {
  if (bcf_write(m_output_fptr, m_template_vcf_hdr, line)) {
    logger.fatal(VCFAdapterException(logger.format("Failed to write VCF/BCF record at position {}, {}",
                                                   bcf_hdr_id2name(m_template_vcf_hdr, line->rid), line->pos+1)));
  }
}

void VCFAdapter::close_file()
{
  if(m_open_output && m_output_fptr) {
    bcf_close(m_output_fptr);
  }
  m_output_fptr = 0;
}

BufferedVCFAdapter::BufferedVCFAdapter(unsigned num_circular_buffers, unsigned max_num_entries)
  : VCFAdapter(true), CircularBufferController(num_circular_buffers) {
  clear();
  m_line_buffers.resize(num_circular_buffers);
  m_num_valid_entries.resize(num_circular_buffers, 0u);
  m_combined_vcf_records_buffer_sizes.resize(num_circular_buffers, 0ull);
  //Initialize buffers
  for (auto& line_buffer : m_line_buffers)
    resize_line_buffer(line_buffer, max_num_entries);
}

BufferedVCFAdapter::~BufferedVCFAdapter() {
  for (auto& line_buffer : m_line_buffers)
    for (auto& line : line_buffer)
      bcf_destroy(line);
  clear();
}

void BufferedVCFAdapter::clear() {
  m_line_buffers.clear();
  m_num_valid_entries.clear();
  m_combined_vcf_records_buffer_sizes.clear();
}

void BufferedVCFAdapter::handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size) {
  auto write_idx = get_write_idx();
  auto& line_buffer = m_line_buffers[write_idx];
  //Need to resize buffer - non-common case
  if (m_num_valid_entries[write_idx] >= line_buffer.size())
    resize_line_buffer(line_buffer, 2u*line_buffer.size()+1u);
  assert(m_num_valid_entries[write_idx] < line_buffer.size());
  std::swap<bcf1_t*>(line, line_buffer[m_num_valid_entries[write_idx]]);
  ++(m_num_valid_entries[write_idx]);
  m_combined_vcf_records_buffer_sizes[write_idx] += bcf_record_size;
}

void BufferedVCFAdapter::resize_line_buffer(std::vector<bcf1_t*>& line_buffer, unsigned new_size) {
  if (new_size <= line_buffer.size())    //never reduce
    return;
  auto curr_idx = line_buffer.size();
  line_buffer.resize(new_size);
  for (auto i=curr_idx; i<line_buffer.size(); ++i)
    line_buffer[i] = bcf_init();
}

void BufferedVCFAdapter::advance_write_idx() {
  auto write_idx = get_write_idx();
  //Advance if something was written
  if (m_num_valid_entries[write_idx])
    CircularBufferController::advance_write_idx();
}

void BufferedVCFAdapter::do_output() {
  if (get_num_entries_with_valid_data() == 0u)
    return;
  auto read_idx = get_read_idx();
  assert(m_num_valid_entries[read_idx] <= m_line_buffers[read_idx].size());
  for (auto i=0u; i<m_num_valid_entries[read_idx]; ++i) {
    assert(m_line_buffers[read_idx][i]);
    auto write_status = bcf_write(m_output_fptr, m_template_vcf_hdr, m_line_buffers[read_idx][i]);
    if (write_status != 0)
      throw VCFAdapterException(std::string("Failed to write VCF/BCF record at position ")
                                +bcf_hdr_id2name(m_template_vcf_hdr, m_line_buffers[read_idx][i]->rid)+", "
                                +std::to_string(m_line_buffers[read_idx][i]->pos+1));
  }
  m_num_valid_entries[read_idx] = 0u;
  m_combined_vcf_records_buffer_sizes[read_idx] = 0ull;
  advance_read_idx();
}

void VCFSerializedBufferAdapter::print_header() {
  assert(m_rw_buffer);
  auto offset = bcf_hdr_serialize(m_template_vcf_hdr, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
                                  m_is_bcf ? 1u : 0u, m_keep_idx_fields_in_bcf_header ? 1u : 0u);
  //Buffer capacity was too small, resize
  while (offset == m_rw_buffer->m_num_valid_bytes) {
    m_rw_buffer->m_buffer.resize(2u*(m_rw_buffer->m_buffer.size())+1u);
    offset = bcf_hdr_serialize(m_template_vcf_hdr, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
                               m_is_bcf ? 1u : 0u, m_keep_idx_fields_in_bcf_header ? 1u : 0u);
  }
  m_rw_buffer->m_num_valid_bytes = offset;
}

void VCFSerializedBufferAdapter::handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size) {
#ifdef DO_PROFILING
  m_vcf_serialization_timer.start();
#endif
  assert(m_rw_buffer);
  auto offset = bcf_serialize(line, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
                              m_is_bcf ? 1u : 0u, m_template_vcf_hdr, &m_hts_string);
  //Buffer capacity was too small, resize
  while (offset == m_rw_buffer->m_num_valid_bytes) {
    m_rw_buffer->m_buffer.resize(2u*(m_rw_buffer->m_buffer.size())+1u);
    offset = bcf_serialize(line, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
                           m_is_bcf ? 1u : 0u, m_template_vcf_hdr, &m_hts_string);
  }
  m_rw_buffer->m_num_valid_bytes = offset;
#ifdef DO_PROFILING
  m_vcf_serialization_timer.stop();
#endif
}

void VCFSerializedBufferAdapter::initialize(const GenomicsDBConfigBase& config) {
  VCFAdapter::initialize(config);
  if (m_do_output)
    m_write_fptr = GenomicsDBConfigBase::output_to_stdout(config.get_vcf_output_filename())
                   ? stdout
                   : fopen(config.get_vcf_output_filename().c_str(), "w");
}

#endif //ifdef HTSDIR
