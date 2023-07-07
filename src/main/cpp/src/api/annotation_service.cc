/**
 * @file annotation_service.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2021 Omics Data Automation, Inc.
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
 * Implementation of the AnnotationService interface to GenomicsDB
 *
 **/

#include "annotation_service.h"

#include "htslib/vcf.h"

#include "genomicsdb_export_config.pb.h"
#include "hfile_genomicsdb.h"
#include "genomicsdb_logger.h"
#include "tiledb_utils.h"
#include "vid_mapper.h"

int set_length_descriptor(FieldLengthDescriptor& length_descr, const std::string& field_type, char *field_number) {
  auto found_ht_map = VidMapper::m_typename_string_to_bcf_ht_type.find(field_type);
  if (found_ht_map == VidMapper::m_typename_string_to_bcf_ht_type.end()) {
    logger.format("Could not map bcf field type {} to vid_mapper ht_type", field_type);
    return -1;
  }
  auto ht_type = found_ht_map->second;
  char* p_end;
  auto field_number_value = std::strtoull(field_number, &p_end, 0);
  if (field_number == p_end || ht_type == BCF_HT_STR) {
    length_descr.set_length_descriptor(0, BCF_VL_VAR);
  } else if (ht_type == BCF_HT_FLAG) {
    length_descr.set_length_descriptor(0, 1);
  } else {
    length_descr.set_num_elements(0, field_number_value);
  }
  return 0;
}

/**
   Read a VCF INFO field definition from the header and return it as a genomic_field_type_t
*/
genomic_field_type_t annotation_source_t::construct_field_type(bcf_hrec_t* hrec) {
    int type_index = bcf_hrec_find_key(hrec, "Type");
    int number_index = bcf_hrec_find_key(hrec, "Number");

    VERIFY2(type_index >= 0 && number_index >=0, "Could not find either Type or Number associated with hrec");

    std::string field_type = hrec->vals[type_index];
    auto found = VidMapper::m_typename_string_to_type_index.find(field_type);
    VERIFY2(found != VidMapper::m_typename_string_to_type_index.end(),
           logger.format("Could not map bcf field type {} to vid_mapper type index", field_type));
    auto genomic_type_index = found->second;

    FieldLengthDescriptor length_descr;
    char *field_number = hrec->vals[number_index];
    auto found_map = VidMapper::m_length_descriptor_string_to_int.find(field_number);
    if (found_map ==  VidMapper::m_length_descriptor_string_to_int.end()) {
      VERIFY2(set_length_descriptor(length_descr, field_type, field_number) == 0,
             logger.format("Could not set_length_descriptor for field type {}", field_type));
    } else {
      FieldInfo field_info;
      FieldElementTypeDescriptor descr(genomic_type_index, found_map->second);
      field_info.set_type(descr);
      length_descr = field_info.m_length_descriptor;
    }

    return genomic_field_type_t(genomic_type_index,
                                length_descr.is_fixed_length_field(),
                                length_descr.is_fixed_length_field()?length_descr.get_num_elements():0,
                                length_descr.get_num_dimensions(),
                                length_descr.is_length_ploidy_dependent()&&length_descr.contains_phase_information());
}

AnnotationService::AnnotationService(const std::string& export_configuration, std::set<std::string> contigs) {
  genomicsdb_pb::ExportConfiguration export_config;
  export_config.ParseFromString(export_configuration);
  for(auto i=0; i<export_config.annotation_source_size(); ++i) {
    auto filename  = export_config.annotation_source(i).filename();
    if (!TileDBUtils::is_file(filename)) {
      throw GenomicsDBException(logger.format("Annotation data source {} not found", filename));
    }

    // If the vcf file is chromosome-specific and the query doesn't cover those chromosomes
    // then don't bother loading the file.
    if(!export_config.annotation_source(i).file_chromosomes().empty()) {
      bool isVcfChromosomeInQuery = false;
      for(auto j=0; j<export_config.annotation_source(i).file_chromosomes_size(); ++j) {
        std::string vcf_contig = export_config.annotation_source(i).file_chromosomes(j);
        if(contigs.find(vcf_contig) != contigs.end()) {
          // The chromosome-specifc vcf and the query have a chromosome in common
          isVcfChromosomeInQuery = true;
          break;
        }
      }

      if(!isVcfChromosomeInQuery) {
      	// Don't load the chromosome-specific VCF 
        continue;
      }
    }

    genomicsdb_htslib_plugin_initialize();

    std::set<std::string> fields;
    for(auto field : export_config.annotation_source(i).attributes()) {
      fields.insert(field);
    }

    std::set<std::string> file_chromosomes;
    for(auto chromosome: export_config.annotation_source(i).file_chromosomes()) {
      file_chromosomes.insert(chromosome);
    }

    m_annotation_sources.emplace_back(filename,
                                      export_config.annotation_source(i).data_source(),
                                      fields,
                                      file_chromosomes);
  }

  m_annotation_buffer_size = export_config.annotation_buffer_size();
  m_annotation_buffer.resize(m_annotation_buffer_size);
}

std::vector<annotation_source_t>& AnnotationService::get_annotation_sources() {
  return m_annotation_sources;
}

std::unordered_map<int, int> bcf_ht_type_to_nbytes = {
  {BCF_HT_INT, sizeof(int)},
  {BCF_HT_INT64, sizeof(long)},
  {BCF_HT_REAL, sizeof(float)},
  {BCF_HT_FLAG, 0},
  {BCF_HT_STR, 1}
};

/**
  Create a new genomic_field_t
 */
genomic_field_t AnnotationService::get_genomic_field(const std::string &data_source,
                                                     const std::string &info_attribute,
                                                     const char *value,
                                                     const int32_t value_length,
                                                     const int bcf_ht_type) {
  std::string genomic_field_name = data_source + "_" + info_attribute;
  size_t bytes = bcf_ht_type_to_nbytes.at(bcf_ht_type)*value_length;
  if (m_annotation_buffer_remaining >=  bytes) {
    size_t start = m_annotation_buffer_size - m_annotation_buffer_remaining;
    memcpy(&m_annotation_buffer[start], value, bytes);
    m_annotation_buffer_remaining -= bytes;
    return genomic_field_t(genomic_field_name, &m_annotation_buffer[start], value_length);
  } else {
    throw GenomicsDBException(logger.format("Annotation Buffer size={} specified to store genomic field info is too small", m_annotation_buffer_size));
  }
}

/**
  Annotate a genomic location using all of the configured data sources. Annotations
  are stored in the genomic_fields vector with the name of the field set to the result of
  concatonating the dataSource (ie ClinVar) with a separator (see AnnotationService.DATA_SOURCE_FIELD_SEPARATOR),
  and the INFO field label.
 */
void AnnotationService::annotate(genomic_interval_t& genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields) {
  assert(m_annotation_buffer_size); // Should have been initialized in constructor
  m_annotation_buffer_remaining = m_annotation_buffer_size;

  for(auto annotation_source: m_annotation_sources) {
    htsFile *htsfile_ptr = hts_open(annotation_source.filename.c_str(), "r");
    VERIFY2(htsfile_ptr!=NULL, logger.format("Could not hts_open {} file in read mode", annotation_source.filename));

    // Only supporting vcf files now
    enum htsExactFormat format = hts_get_format(htsfile_ptr)->format;
    VERIFY2(format==vcf, logger.format("File {} is not VCF (a)", annotation_source.filename));

    bcf_hdr_t *hdr = bcf_hdr_read(htsfile_ptr);
    VERIFY2(hdr!=NULL, logger.format("Could not read header in file {}", annotation_source.filename));

    // Read the tbi index file to help with reading the annotation vcfs
    tbx_t *tbx = tbx_index_load3(annotation_source.filename.c_str(), NULL, 0);
    VERIFY2(tbx, logger.format("Could not load the index file for {}", annotation_source.filename));

    // Query using chromosome and position range
    std::string query_range = genomic_interval.contig_name + ":" + std::to_string(genomic_interval.interval.first) + "-" + std::to_string(genomic_interval.interval.second);
    hts_itr_t *itr = tbx_itr_querys(tbx, query_range.c_str());

    VERIFY2(itr, "Could not obtain tbx query iterator. Possibly caused by the vcf not having any variants on the requested chromosome.");

    regidx_t *reg_idx = NULL;
    // reg_idx = regidx_init(annotation_source.filename().c_str(), NULL, NULL, 0, NULL);
    // VERIFY2(reg_idx, "Unable to read file");

    const char **seq = NULL;
    int nseq;
    if (reg_idx) {
      seq = tbx_seqnames(tbx, &nseq);
    }

    kstring_t str = {0,0,0};
    bcf1_t *rec = bcf_init1();
    int idx = 0;

    // Iterate over each matching position in the VCF
    while (tbx_itr_next(htsfile_ptr, tbx, itr, &str) >= 0) {
      if (reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) {
        continue;
      }

      VERIFY2(vcf_parse1(&str, hdr, rec) == 0, "Problem parsing current line of VCF");

      bcf_unpack((bcf1_t*)rec, BCF_UN_ALL); // Using BCF_UN_INFO is probably a little faster
      // bcf_unpack((bcf1_t*)rec, BCF_UN_INFO); // Using BCF_UN_INFO is probably a little faster

      if(ref.compare(rec->d.allele[0]) != 0) {
        // REF doesn't match
        continue;
      } else if(alt.compare(rec->d.allele[1])) {
        // ALT doesn't match
        // NOTE: This only looks at one allele. If there are multiple alternate alleles you need to update the code.
        continue;
      }

      // iteration over the list of INFO fields we are interested in
      for(auto info_attribute: annotation_source.fields) {
        if(info_attribute == "ID") {
          genomic_fields.push_back(get_genomic_field(annotation_source.datasource, info_attribute, rec->d.id, strlen(rec->d.id)));
        } else {
          void* info_value = NULL;
          int32_t info_value_length = 0;
          int info_index = bcf_hdr_id2int(hdr, BCF_DT_ID, info_attribute.c_str());
          VERIFY2(info_index != -1, logger.format("Info field {} not found in annotation source {}", info_attribute,
                                                 annotation_source.filename));
          int info_type = bcf_hdr_id2type(hdr, BCF_HL_INFO, info_index);
          int num_values = bcf_get_info_values(hdr, rec, info_attribute.c_str(), &info_value, &info_value_length, info_type);
          VERIFY2((num_values >= 0 || num_values == -3),
                 logger.format("bcf_get_info_values returned {}. See vcf.h for error code", num_values));
          if (num_values > 0) {
            if (info_type == BCF_HT_STR && ((char *)info_value)[info_value_length-1] == 0) {
              info_value_length--;
            }
            genomic_fields.push_back(get_genomic_field(annotation_source.datasource, info_attribute,
                                                       (char *)info_value, info_value_length, info_type));
            free(info_value);
          }
        }
      }
    }

   regidx_destroy(reg_idx);
   bcf_itr_destroy(itr);
   bcf_hdr_destroy(hdr);
   hts_close(htsfile_ptr);
  }
}
