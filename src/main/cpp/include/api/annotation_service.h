/**
 * @file annotation_service.h
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
 * Specification of the AnnotationService to GenomicsDB
 *
 **/

#include "genomicsdb.h"
#include "logger.h"

#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "htslib/regidx.h"
#include "htslib/tbx.h"

#include <map>
#include <string>
#include <vector>

#include "genomicsdb_export_config.pb.h"

#define VERIFY2(X, Y) if(!(X)) throw GenomicsDBException(Y);

typedef struct annotation_source_t {
  std::string filename;
  std::string datasource;
  std::set<std::string> fields;

  bool field_types_processed = false;

  // Constructor
  annotation_source_t(const std::string& filename, const std::string& datasource,
                      const std::set<std::string>& fields) {
    this->filename = std::move(filename);
    this->datasource = std::move(datasource);
    this->fields = std::move(fields);
  }

  std::map<std::string, genomic_field_type_t> field_types() {
    htsFile *htsfile_ptr = hts_open(this->filename.c_str(), "r");
    VERIFY2(htsfile_ptr!=NULL, logger.format("Could not hts_open {} file in read mode", filename));

    // Only supporting vcf files now
    enum htsExactFormat format = hts_get_format(htsfile_ptr)->format;
    VERIFY2(format==vcf, logger.format("File {} is not VCF", filename));

    bcf_hdr_t *hdr = bcf_hdr_read(htsfile_ptr);
    VERIFY2(hdr!=NULL, logger.format("Could not read header in file {}", filename));

    std::map<std::string, genomic_field_type_t> genomic_field_types;
    std::set<std::string> processed_fields;

    for (int i=0; i<hdr->nhrec; i++) {
      bcf_hrec_t* hrec = hdr->hrec[i];
      VERIFY2(hrec!=NULL, logger.format("Could not get header records from {}", filename));
      construct_field_type(hrec, genomic_field_types, processed_fields);
    }

    bcf_hdr_destroy(hdr);
    hts_close(htsfile_ptr);

    if (fields.size() != processed_fields.size()) {
      for (auto field: fields) {
        if (processed_fields.find(field) == processed_fields.end()) {
          VERIFY2(field.compare("ID") == 0,
                  logger.format("No type information could be deduced for annotation info field {}", field));
          genomic_field_types.insert(std::make_pair(datasource + "_" + field,
                                                    genomic_field_type_t(typeid(char),
                                                                         false,/*is_fixed_num_elements*/
                                                                         1,/*num_elements*/
                                                                         1,/*num_dimensions*/
                                                                         false/*contains_phase_info*/)));
        }
      }
    }

    return genomic_field_types;
  }

 private:
  void construct_field_type(bcf_hrec_t* hrec, std::map<std::string, genomic_field_type_t>& types,
                   std::set<std::string>& processed_fields) {
    int index;
    std::string field;
    if (hrec->type == BCF_HL_INFO) {
      index = bcf_hrec_find_key(hrec, "ID");
      if (index >= 0) {
        field = hrec->vals[index];
        if (fields.find(field) != fields.end()) {
          field = hrec->vals[index];
          processed_fields.insert(field);
          types.insert(std::make_pair(datasource + "_" + field, construct_field_type(hrec)));
        }
      }
    }
  }

  genomic_field_type_t construct_field_type(bcf_hrec_t* hrec);
} annotation_source_t;

/**
  Use this service to add annotations from VCF datasources to variant's genomic fields
*/
class AnnotationService {
 public:
  /**
   * Read in annotation sources from ExportConfiguration.
   **/
  AnnotationService(const std::string& export_configuration);

  std::vector<annotation_source_t>& get_annotation_sources();

  void annotate(genomic_interval_t &genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields);

  genomic_field_t get_genomic_field(const std::string &data_source, const std::string &info_attribute, const char *value, const int32_t value_length, int bcf_ht_type=BCF_HT_STR);

 private:
  // List of configured annotation data sources
  std::vector<annotation_source_t> m_annotation_sources;

  // Buffer for annotation values
  // TODO: This could be per attribute/field
  std::vector<std::shared_ptr<std::string>> m_annotation_buffer;
};
