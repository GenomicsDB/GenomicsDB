/**
 * @file genomicsdb_nanoarrow_processor.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2024 dātma, inc™
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
 * GenomicsDB pyarrow output
 *
 **/
#include "genomicsdb.h"
#include "genomicsdb_logger.h"

#include <nanoarrow/nanoarrow.hpp>

#include <mutex>

#define TO_ARROW_SCHEMA(X) (reinterpret_cast<ArrowSchema *>(X))
#define TO_ARROW_ARRAY(X) (reinterpret_cast<ArrowArray *>(X))

ArrowVariantCallProcessor::ArrowVariantCallProcessor(bool non_blocking) : m_non_blocking(non_blocking) {
  m_arrow_schema = new ArrowSchema();
}

void ArrowVariantCallProcessor::cleanup_schema(void* schema) {
  if (schema) {
    ArrowSchema* arrow_schema = TO_ARROW_SCHEMA(schema);
    if (arrow_schema->release) {
      ArrowSchemaRelease(arrow_schema);
    }
    delete arrow_schema;
  }
}

void ArrowVariantCallProcessor::cleanup_array(void* array) {
  if (array) {
    ArrowArray* arrow_array = TO_ARROW_ARRAY(array);
    if (arrow_array->n_children) {
      if (arrow_array->release) {
        ArrowArrayRelease(arrow_array);
      }
      delete arrow_array;
    }
  }
}

int set_schema(ArrowSchema* schema, const std::string& name, ArrowType type) {
  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetType(schema, type));
  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetName(schema, name.c_str()));
  return NANOARROW_OK;
}

int create_arrow_schema(ArrowSchema* schema,
                        std::shared_ptr<std::map<std::string, genomic_field_type_t>> genomic_field_types,
                        std::map<std::string, int64_t>& genomic_fields_map) {
  ArrowSchemaInit(schema);
  ArrowSchemaSetName(schema, "Main");

  int64_t num_children = 3; // sample_name, chr, pos
  for (auto& field_type_pair : *genomic_field_types) {
    std::string field_name = field_type_pair.first;
    if (field_name.compare("END") && field_name.compare("ALT")) {
      num_children++;
    }
  }
  
  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetTypeStruct(schema, num_children));

  int64_t i = 0;
  int rc = 0;
  rc |= set_schema(schema->children[i++], "SAMPLE_NAME", NANOARROW_TYPE_STRING);
  rc |= set_schema(schema->children[i++], "CHR", NANOARROW_TYPE_STRING);
  rc |= set_schema(schema->children[i++], "POS", NANOARROW_TYPE_UINT64);

  for (auto& field_type_pair : *genomic_field_types) {
    if (rc) break;
    std::string field_name = field_type_pair.first;
    auto field_type = field_type_pair.second;
    if (field_name.compare("END") && field_name.compare("ALT")) {
      genomic_fields_map.insert(std::pair<std::string,int64_t>(field_name, i));
      if (STRING_FIELD(field_name, field_type)) {
        rc |= set_schema(schema->children[i++], field_name, NANOARROW_TYPE_STRING);
      } else if (INT_FIELD(field_type)) {
        rc |= set_schema(schema->children[i++], field_name, NANOARROW_TYPE_INT32);
      } else if (FLOAT_FIELD(field_type)) {
        rc |= set_schema(schema->children[i++], field_name, NANOARROW_TYPE_FLOAT);
      }
    }
  }
  if (rc) {
    ArrowVariantCallProcessor::cleanup_schema(schema);
    logger.fatal(GenomicsDBException(), "Could not set arrow schema rc={} errno={} {}", rc, errno, errno?strerror(errno):"");
  }
  
  return NANOARROW_OK;
}

ArrowArray* build_arrow_array(void* array, ArrowSchema* schema) {
  //  ArrowVariantCallProcessor::cleanup_array(array);
  ArrowArray* arrow_array = new ArrowArray();
  int rc = ArrowArrayInitFromSchema(arrow_array, schema, NULL);
  for (auto i=0l; i<arrow_array->n_children; i++) {
    rc |= ArrowArrayStartAppending(arrow_array->children[i]);
  }
  if (rc) {
    ArrowVariantCallProcessor::cleanup_schema(schema);
    logger.fatal(GenomicsDBException(), "Could not build arrow array for given schema rc={} errno={} {}", rc, errno, errno?strerror(errno):"");
  }
  return arrow_array;
}

void ArrowVariantCallProcessor::process(const interval_t& interval) {
  if (!m_is_initialized) {
    m_is_initialized = true;

    ArrowSchema* arrow_schema = TO_ARROW_SCHEMA(m_arrow_schema);
    int rc = create_arrow_schema(arrow_schema, get_genomic_field_types(), m_genomic_fields_map);
    if (rc) {
      ArrowVariantCallProcessor::cleanup_schema(m_arrow_schema);
      logger.fatal(GenomicsDBException(), "Could not create arrow schema rc={} errno={} {}", rc, errno, errno?strerror(errno):"");
    }
    m_arrow_array = build_arrow_array(TO_ARROW_ARRAY(m_arrow_array), arrow_schema);
  }
}

int ArrowVariantCallProcessor::append_arrow_array(const std::string& sample_name,
                       const int64_t* coordinates,
                       const genomic_interval_t& genomic_interval,
                       const std::vector<genomic_field_t>& genomic_fields) {
  ArrowArray* array = TO_ARROW_ARRAY(m_arrow_array);

  int64_t i = 0;
  NANOARROW_RETURN_NOT_OK(ArrowArrayAppendString(array->children[i++], ArrowCharView(sample_name.c_str())));
  NANOARROW_RETURN_NOT_OK(ArrowArrayAppendString(array->children[i++], ArrowCharView(genomic_interval.contig_name.c_str())));
  NANOARROW_RETURN_NOT_OK(ArrowArrayAppendInt(array->children[i++], genomic_interval.interval.first));
  auto num_constant_fields = i;

  std::vector<bool> found(m_genomic_fields_map.size()+num_constant_fields, false);
  for (auto& genomic_field : genomic_fields) {
    std::string field_name = genomic_field.name;
    if (field_name.compare("END") && field_name.compare("ALT")) {
      auto field_type = get_genomic_field_types()->at(field_name);
      i = m_genomic_fields_map[field_name];
      found[i] = true;
      if (STRING_FIELD(field_name, field_type)) {
        if (field_name == "GT") {
          NANOARROW_RETURN_NOT_OK(ArrowArrayAppendString(array->children[i],
                                                         ArrowCharView(resolve_gt(genomic_fields).c_str())));
        } else {
          NANOARROW_RETURN_NOT_OK(ArrowArrayAppendString(array->children[i], ArrowCharView( genomic_field.to_string(field_type).c_str())));
        }
      } else if (INT_FIELD(field_type)) {
        NANOARROW_RETURN_NOT_OK(ArrowArrayAppendInt(array->children[i], genomic_field.int_value_at(0)));
      } else if (FLOAT_FIELD(field_type)) {
        NANOARROW_RETURN_NOT_OK(ArrowArrayAppendInt(array->children[i], genomic_field.float_value_at(0)));
      }
    }
  }

  for (auto j=num_constant_fields; j<found.size() && !found[j]; j++) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayAppendNull(array->children[j], 1));
  }
    
  return NANOARROW_OK;
}

void ArrowVariantCallProcessor::process(const std::string& sample_name,
                                        const int64_t* coordinates,
                                        const genomic_interval_t& genomic_interval,
                                        const std::vector<genomic_field_t>& genomic_fields) {
  if (m_non_blocking && m_last_column != -1) {
    if (m_last_column != coordinates[1]) {
      full.release();
    } else {
      empty.release();
    }
  }
  if (m_non_blocking) empty.acquire();
  {
    std::lock_guard<std::mutex> g(m_mtx);
    m_last_column = coordinates[1];
    if (append_arrow_array(sample_name, coordinates, genomic_interval, genomic_fields)) {
      logger.error("Could not append sample_name={}", sample_name);
    }
  }
}

void ArrowVariantCallProcessor::finalize() {
  m_is_finalized = true;
  if (m_non_blocking) full.release();
}

void *ArrowVariantCallProcessor::arrow_schema() {
  return m_arrow_schema;
}

void *ArrowVariantCallProcessor::arrow_array() {
  if (m_is_finalized && TO_ARROW_ARRAY(m_arrow_array)->children[0]->length == 0) return NULL;
  if (m_non_blocking) full.acquire();
  ArrowArray *arrow_array = new ArrowArray();
  ArrowArrayInitFromSchema(arrow_array, TO_ARROW_SCHEMA(m_arrow_schema), NULL);
  {
    std::lock_guard<std::mutex> g(m_mtx);
    if (m_arrow_array) {
      ArrowArrayFinishBuildingDefault(TO_ARROW_ARRAY(m_arrow_array), NULL);
      for (auto i=0l; i<TO_ARROW_ARRAY(m_arrow_array)->n_children; i++) {
        ArrowArrayMove(TO_ARROW_ARRAY(m_arrow_array)->children[i], arrow_array->children[i]);
      }
      m_arrow_array = nullptr;
      m_arrow_array = build_arrow_array(TO_ARROW_ARRAY(m_arrow_array), TO_ARROW_SCHEMA(m_arrow_schema));
    }
  }
  if (m_non_blocking) empty.release();
  return arrow_array;
}

/** Helper arrow utilities for Python and other bindings */
int ArrowVariantCallProcessor::allocate_schema(void** schema, void *src) {
  int rc = 0;
  ArrowSchema** arrow_schema = reinterpret_cast<ArrowSchema**>(schema);
  *arrow_schema = (ArrowSchema *)(ArrowMalloc(sizeof(ArrowSchema)));
  if (src) {
    rc = ArrowSchemaDeepCopy(TO_ARROW_SCHEMA(src), *arrow_schema);
  }
  return rc;
}


