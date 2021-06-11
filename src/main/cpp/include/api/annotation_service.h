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

#include "htslib/regidx.h"
#include "htslib/tbx.h"

#include "genomicsdb.h"

#include <map>
#include <string>
#include <vector>

#include "genomicsdb_export_config.pb.h"

#define VERIFY(X) if(!(X)) throw GenomicsDBException(#X);

typedef struct annotation_source_t {
  std::string filename;
  std::string datasource;
  std::vector<std::string> fields;
  // Constructor
  annotation_source_t(std::string filename, std::string datasource,
                      std::vector<std::string> fields) {
    this->filename = std::move(filename);
    this->datasource = std::move(datasource);
    this->fields = std::move(fields);
  }
  std::map<std::string, genomic_field_type_t> field_types() {
    std::map<std::string, genomic_field_type_t> field_types;
    for(std::string field: this->fields) {
      const std::string attribute_name = this->datasource + "_" + field;
      // TODO : support types other than string too - see src/resources/genomicsdb_vid_mapping.proto
      const std::type_index type_index = typeid(char);
      field_types.insert(std::make_pair(attribute_name,
                                        genomic_field_type_t(typeid(char),
                                                             false,/*is_fixed_num_elements*/
                                                             1,/*num_elements*/
                                                             1,/*num_dimensions*/
                                                             false/*contains_phase_info*/)));
    }
    return field_types;
  }
} annotation_source_t;

/**
  Use this service to add annotations from VCF datasources to variant's genomic ields
*/
class AnnotationService {
 public:
  /**
   * Read in annotation sources from ExportConfiguration.
   **/
  AnnotationService(const std::string& export_configuration);

  std::vector<annotation_source_t>& get_annotation_sources();

  void annotate(genomic_interval_t &genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields);

  genomic_field_t get_genomic_field(const std::string &data_source, const std::string &info_attribute, const char *value, const int32_t value_length);

 private:
  // List of configured annotation data sources
  std::vector<annotation_source_t> m_annotation_sources;

  // Buffer for annotation values
  // TODO: This could be per attribute/field
  std::vector<std::string> m_annotation_buffer;
};
