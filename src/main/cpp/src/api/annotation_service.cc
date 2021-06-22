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
  m_annotation_buffer.reserve(32);
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
 * Read an INFO value for a record. The value is stored in string pointer info_value, and the return
 * value of the htslib bcf_get_info_* function is returned.
 * List of return codes:
 * >0 .. success, value found
 * =0 .. success, but row does not have the requested field
 * -1 .. no such INFO tag defined in the header
 * -2 .. clash between types defined in the header and encountered in the VCF record
 * -3 .. tag is not present in the VCF record
 * -4 .. the operation could not be completed (e.g. out of memory)
**/
int AnnotationService::get_info_value(bcf_hdr_t *hdr, bcf1_t *rec, std::string *info_attribute, std::string *info_value) {
  int32_t info_value_length = 0;
  int32_t *info_value_ptr = NULL;

  // Determine INFO field type
  bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", info_attribute->c_str(), NULL);
  if(hrec == NULL) {
    return -1;
  }

  int type_index = bcf_hrec_find_key(hrec, "Type");
  std::string field_type = hrec->vals[type_index];

  // Retrieve the value using the bcf_get_info_* function corresponding to INFO field type
  //  - Possible Types for INFO fields are: Integer, Float, Flag, Character, and String.
  //  - The Number entry is an Integer that describes the number of values that can be included with the INFO eld
  int res;
  if (field_type.compare("String") == 0) {
    // jDebug: it seems odd that this pointer needs to be set to NULL. But if i don't then i get an exception at run time:
    // api_tests(50033,0x11c527dc0) malloc: *** error for object 0x7ffee2c7dbf8: pointer being realloc'd was not allocated
    // api_tests(50033,0x11c527dc0) malloc: *** set a breakpoint in malloc_error_break to debug
    // SIGABRT - Abort (abnormal termination) signal
    res = bcf_get_info_string(hdr, rec, info_attribute->c_str(), &info_value_ptr, &info_value_length);
    if(res > 0) {
      info_value->assign((char*)info_value_ptr, info_value_length);
    }
  } else if (field_type.compare("Integer") == 0) {
    res = bcf_get_info_int32(hdr, rec, info_attribute->c_str(), &info_value_ptr, &info_value_length);
    if(res > 0) {
      int integer = *info_value_ptr;
      std::string intStr = std::to_string(integer);
      info_value->assign(intStr);
    }
  } else if (field_type.compare("Float") == 0) {
    res = bcf_get_info_float(hdr, rec, info_attribute->c_str(), &info_value_ptr, &info_value_length);
    if(res > 0) {
      float f = *((float*)info_value_ptr);
      std::string floatStr = std::to_string(f);
      info_value->assign(floatStr);
    }
  } else if (field_type.compare("Flag") == 0) {
    res = bcf_get_info_flag(hdr, rec, info_attribute->c_str(), &info_value_ptr, &info_value_length);
    if(res > 0) {
      // Flags don't have values. They are present or not.
      info_value->assign("");
    }
  } else {
    throw GenomicsDBException(logger.format("INFO field {} is of unknown type {}", *info_attribute, field_type));
  }

  return res;
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
      // jDebug: what is this condition?
      if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) {
        continue;
      }

      int readResult = vcf_parse1(&str, hdr, rec);
      VERIFY(readResult==0 && "Problem parsing current line of VCF");

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
          // If the request is for "ID" then give the VCF row ID
          if(info_attribute == "ID") {
            genomic_fields.push_back(get_genomic_field(annotation_source.datasource, info_attribute, rec->d.id, strlen(rec->d.id)));
          } else {
            std::string info_value;
            int res = get_info_value(hdr, rec, &info_attribute, &info_value);

            // Look for the info field. See get_genomic_field for return code meanings
            if(res > 0) {
              genomic_field_t genomic_field_annotation = get_genomic_field(annotation_source.datasource,
                                                                     info_attribute, info_value.c_str(), info_value.length());
              genomic_fields.push_back(std::move(genomic_field_annotation));
            } else if (res == 0 || res == -3) {
              // Do nothing - valid query, but the INFO attribute is not present in this row.
            } else {
              throw GenomicsDBException(logger.format("Recieved error code {} while fetching INFO field {}", res, info_attribute));
            }
          }
      }
   }

   regidx_destroy(reg_idx);
   bcf_itr_destroy(itr);

   int ret = hts_close(vcfInFile);
   VERIFY(ret == 0 && "Non-zero status when trying to close VCF file");
  }
}
