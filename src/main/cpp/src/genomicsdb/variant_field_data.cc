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

#include "variant_field_data.h"

std::string g_vcf_NON_REF="<NON_REF>";

fi_union g_tiledb_null_float = { .f = TILEDB_EMPTY_FLOAT32 };
di_union g_tiledb_null_double = { .d = TILEDB_EMPTY_FLOAT64 };

std::unordered_map<std::type_index, int> g_variant_field_type_index_to_tiledb_type =
std::unordered_map<std::type_index, int> {
  { std::type_index(typeid(int)), TILEDB_INT32 },
  { std::type_index(typeid(int64_t)), TILEDB_INT64 },
  { std::type_index(typeid(unsigned)), TILEDB_INT32 },
  { std::type_index(typeid(uint64_t)), TILEDB_INT64 },
  { std::type_index(typeid(float)), TILEDB_FLOAT32 },
  { std::type_index(typeid(double)), TILEDB_FLOAT64 },
  { std::type_index(typeid(char)), TILEDB_CHAR },
  { std::type_index(typeid(bool)), TILEDB_CHAR }
};

std::vector<std::type_index> g_tiledb_type_to_variant_field_type_index = {
  std::type_index(typeid(int)),
  std::type_index(typeid(int64_t)),
  std::type_index(typeid(float)),
  std::type_index(typeid(double)),
  std::type_index(typeid(char))
};

std::unordered_map<std::type_index, int> g_variant_field_type_index_to_bcf_ht_type =
std::unordered_map<std::type_index, int> {
  { std::type_index(typeid(void)), BCF_HT_VOID },
  { std::type_index(typeid(bool)), BCF_HT_FLAG },
  { std::type_index(typeid(int)), BCF_HT_INT },
  { std::type_index(typeid(int64_t)), BCF_HT_INT64 },
  { std::type_index(typeid(unsigned)), BCF_HT_UINT },
  { std::type_index(typeid(uint64_t)), BCF_HT_UINT64 },
  { std::type_index(typeid(float)), BCF_HT_REAL },
  { std::type_index(typeid(double)), BCF_HT_DOUBLE },
  { std::type_index(typeid(std::string)), BCF_HT_STR },
  { std::type_index(typeid(char)), BCF_HT_CHAR }
};


size_t VariantFieldTypeUtil::size(const int bcf_ht_type) {
  switch(bcf_ht_type) {
  case BCF_HT_INT:
    return sizeof(int);
  case BCF_HT_INT64:
    return sizeof(int64_t);
  case BCF_HT_UINT:
    return sizeof(unsigned);
  case BCF_HT_UINT64:
    return sizeof(uint64_t);
  case BCF_HT_REAL:
    return sizeof(float);
  case BCF_HT_DOUBLE:
    return sizeof(double);
  case BCF_HT_STR:
  case BCF_HT_CHAR:
    return sizeof(char);
  default:
    return 0;
  }
  return 0;
}
