/**
 * @file genomicsdb_field.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2020,2022 Omics Data Automation, Inc.
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
 * Implementation of GenomicsDB field api methods
 *
 **/

#include "genomicsdb.h"
#include "gt_common.h"

#include <iostream>
#include <sstream>
#include <string>

#define NON_REF "<NON_REF>"

std::string genomic_field_t::recombine_ALT_value(const genomic_field_type_t& field_type, std::string separator) const {
  std::string output;
  std::string item;
  if(field_type.is_cppstring()) {
    for(unsigned int i = 0; i < num_elements; i++) {
      item = cpp_str_value_at(i);
      if (IS_NON_REF_ALLELE(item)) {
        output.empty()?output=NON_REF:output=output+separator+NON_REF;
      } else {
        output.empty()?output=item:output=output+separator+item;
      }
    }
  } else {
    // otherwise treat as char*
    std::stringstream ss(str_value());
    while (std::getline(ss, item, PHASED_ALLELE_SEPARATOR)){
      if (IS_NON_REF_ALLELE(item)) {
        output.empty()?output=NON_REF:output=output+separator+NON_REF;
      } else {
        output.empty()?output=item:output=output+separator+item;
      }
    }
  }

  return "[" + output + "]";
}

std::string genomic_field_t::combine_GT_vector(const genomic_field_type_t& field_type) const {
  assert(field_type.is_int());
  
  std::string output;
  if (field_type.contains_phase_information()) {
    for (auto i=0ul; i<num_elements; i++) {
      output = output + (int_value_at(i)==get_bcf_gt_no_call_allele_index<int>()?".":to_string(i, field_type));
      if (i+1<num_elements) {
        if (int_value_at(i+1) == 0) {
          output = output + UNPHASED_ALLELE_SEPARATOR;
        } else {
          output = output + PHASED_ALLELE_SEPARATOR;
        }
        i++;
      }
    }
  } else {
    for (auto i=0ul; i<num_elements; i++) {
      std::string gt_element = to_string(i, field_type);
      output = output + (int_value_at(i)==get_bcf_gt_no_call_allele_index<int>()?".":to_string(i, field_type));
      if (i+1<num_elements) {
        output = output + UNPHASED_ALLELE_SEPARATOR;
      }
    }
  }
  return output;
}
