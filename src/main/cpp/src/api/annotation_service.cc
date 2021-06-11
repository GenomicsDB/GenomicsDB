/**
 * @file annotation_service.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2021 Omics Data Automation, Inc.
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
#include "logger.h"
#include "tiledb_utils.h"

AnnotationService::AnnotationService(const std::string& export_configuration) {
  genomicsdb_pb::ExportConfiguration export_config;
  export_config.ParseFromString(export_configuration);
  for(auto i=0; i<export_config.annotation_source_size(); ++i) {
    auto filename  = export_config.annotation_source(i).filename();
    if (!TileDBUtils::is_file(filename)) {
      throw GenomicsDBException(logger.format("Annotation data source {} not found", filename));
    }
    genomicsdb_htslib_plugin_initialize(filename.c_str());
    std::vector<std::string> fields;
    for(auto field : export_config.annotation_source(i).attributes()) {
      fields.push_back(field);
    }
    annotation_source_t source(filename, export_config.annotation_source(i).data_source(), fields);
    m_annotation_sources.push_back(std::move(source));
  }
  m_annotation_buffer.reserve(16);
}

std::vector<annotation_source_t>& AnnotationService::get_annotation_sources() {
  return m_annotation_sources;
}

/**
  Create a new genomic_field_t
 */
genomic_field_t AnnotationService::get_genomic_field(const std::string &data_source,
                                                     const std::string &info_attribute,
                                                     const char *value,
                                                     const int32_t value_length) {
  std::string genomic_field_name = data_source + "_" + info_attribute;
  m_annotation_buffer.push_back(std::string(value, value_length));
  return genomic_field_t(genomic_field_name, m_annotation_buffer.back().data(),
                         m_annotation_buffer.back().length());
}

/**
  Annotate a genomic location using all of the configured data sources. Annotations
  are stored in the genomic_fields vector with the name of the field set to the result of
  concatonating the dataSource (ie ClinVar) with a separator (see AnnotationService.DATA_SOURCE_FIELD_SEPARATOR),
  and the INFO field label.
 */
void AnnotationService::annotate(genomic_interval_t& genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields) {
  for(auto annotation_source: m_annotation_sources) {
    // Open the VCF file
    htsFile *vcfInFile = hts_open(annotation_source.filename.c_str(), "r");
    VERIFY(vcfInFile != NULL && "Unable to open VCF file");

    // Check file type
    enum htsExactFormat format = hts_get_format(vcfInFile)->format;
    VERIFY(format == vcf && "File is not VCF");

    // Read the header
    bcf_hdr_t *hdr = bcf_hdr_read(vcfInFile);
    VERIFY(hdr && "Could not read VCF header");

    // I'm not sure what this does. Whatever it is, it takes a really long time.
    // Need to look in to how to remove this variable.
    regidx_t *reg_idx = NULL;
    // reg_idx = regidx_init(annotation_source.filename().c_str(), NULL, NULL, 0, NULL);
    // VERIFY(reg_idx && "Unable to read file");

    // Read the tbi index file
    tbx_t *tbx = tbx_index_load3(annotation_source.filename.c_str(), NULL, 0);
    VERIFY(tbx && "Could not load .tbi index");

    // Query using chromosome and position range
    std::string variantQuery = genomic_interval.contig_name + ":" + std::to_string(genomic_interval.interval.first) + "-" + std::to_string(genomic_interval.interval.second);

    hts_itr_t *itr = tbx_itr_querys(tbx, variantQuery.c_str());
    VERIFY(itr && "Problem opening iterator");

    const char **seq = NULL;
    int nseq;
    if (reg_idx) {
      seq = tbx_seqnames(tbx, &nseq);
    }

    kstring_t str = {0,0,0};
    bcf1_t *rec    = bcf_init1();
    int idx = 0;

   
    // Iterate over each matching position in the VCF
    while (tbx_itr_next(vcfInFile, tbx, itr, &str) >= 0)
    {
      if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) {
        continue;
      }

      int readResult = vcf_parse1(&str, hdr, rec);
      VERIFY(readResult==0 && "Problem parsing current line of VCF");

      // jDebug: update this comment. It says "need to add code here", and i think you already have.
      // The query is constrained by chromosome and position but not by allele.  Need to add code here which
      // compares the matches alt and ref alleles to see if they match what we are looking for.
      bcf_unpack((bcf1_t*)rec, BCF_UN_ALL); // Using BCF_UN_INFO is probably a little faster

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
        char* info_value = NULL;
        int32_t info_value_length = 0;

        // bcf_info_t *info = rec->d.info;

          // If the request is for "ID" then give the VCF row ID
          if(info_attribute == "ID") {
            genomic_fields.push_back(get_genomic_field(annotation_source.datasource, info_attribute, rec->d.id, strlen(rec->d.id)));
          } else {
            // Look for the info field.
            // List of return codes:
            // * -1 .. no such INFO tag defined in the header
            // * -2 .. clash between types defined in the header and encountered in the VCF record
            //            - Seems like there should be a way to query to determine what the INFO type is.
            // * -3 .. tag is not present in the VCF record
            int res = bcf_get_info_string(hdr, rec, info_attribute.c_str(), &info_value, &info_value_length);

            if(res > 0) {            
              genomic_field_t jDebugGenomicField = get_genomic_field(annotation_source.datasource,
                                                                     info_attribute, info_value, info_value_length);
              // genomic_fields.push_back(get_genomic_field(annotation_source.data_source(), info_attribute, infoValue, infoValueLength));
              genomic_fields.push_back(std::move(jDebugGenomicField)); // jDebug: delete this line and restore the one above
              // printf("jDebug.cc.a: %s %s=%s\n", annotation_source.data_source().c_str(), info_attribute.c_str(), infoValue);
              // printf("jDebuc.cc.b:    %s=%s\n", jDebugGenomicField.name.c_str(), jDebugGenomicField.str_value().c_str());
            } else if (res == -2) {
              // "clash between types defined in the header and encountered in the VCF record"
              // We need to modify this section so it determines the info field type and then uses the correct
              // bcf_get_info_* method (eg bcf_get_info_int32, bcf_get_info_float, bcf_get_info_flag).

              printf("JDEBUG: Need to deal with error -2, info_attribute=%s\n", info_attribute.c_str());
            } else {
              // jDebug: throw an exception? No, i think it is ok for field not to be found.
              // info field not found
              printf("JDEBUG: Need to deal with ELSE error, info_attribute=%s\n", info_attribute.c_str());
            }
          }
          free(info_value);
      }
   }

   regidx_destroy(reg_idx);

   // jDebug: this function wasn't found so there's probably a leak: regitr_destroy(itr);
   // jDebug: need to find out if the iterator needs to be destroyed.
   bcf_itr_destroy(itr);

   int ret = hts_close(vcfInFile);
   VERIFY(ret == 0 && "Non-zero status when trying to close VCF file");
  }
}
