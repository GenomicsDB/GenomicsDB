/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#ifndef GT_COMMON_H
#define GT_COMMON_H

#include "headers.h"
#include "profiling.h"

#include "vcf.h"

#define CHECK_MISSING_SAMPLE_GIVEN_REF(REF) (((REF).size() == 0) || ((REF)[0] == '$'))
#define CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(REF) ((REF).size() == 0 || ((REF) == ""))
#define CHECK_IN_THE_MIDDLE_REF(REF) ((REF)[0] == 'N')

extern std::string g_vcf_NON_REF;
extern std::string g_vcf_SPANNING_DELETION;

inline bool IS_NON_REF_ALLELE(const std::string& allele) {
  return allele.length() > 0 && ((allele)[0] == '&');
}

inline bool IS_NON_REF_ALLELE(const char allele_char) {
  return allele_char == '&';
}

#define TILEDB_NON_REF_VARIANT_REPRESENTATION "&"
#define TILEDB_ALT_ALLELE_SEPARATOR "|"

#define PHASED_ALLELE_SEPARATOR '|'
#define UNPHASED_ALLELE_SEPARATOR '/'

#define MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED 50u
#define MAX_GENOTYPE_COUNT ((unsigned)INT32_MAX) //this is likely larger than what can be handled by the VCF format
#define DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE 1048576u

#define UNDEFINED_ATTRIBUTE_IDX_VALUE 0xFFFFFFFFu
#define UNDEFINED_NUM_ROWS_VALUE 0xFFFFFFFFFFFFFFFFull
#define UNDEFINED_UINT64_T_VALUE 0xFFFFFFFFFFFFFFFFull

extern std::unordered_map<std::type_index, int> g_variant_field_type_index_to_tiledb_type;
extern std::unordered_map<std::type_index, int> g_variant_field_type_index_to_bcf_ht_type;
extern std::vector<std::type_index> g_tiledb_type_to_variant_field_type_index;

extern std::string g_tmp_scratch_dir;

extern const char* g_json_indent_unit;

#endif
