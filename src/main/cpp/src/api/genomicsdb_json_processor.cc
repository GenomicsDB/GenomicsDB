/**
 * @file genomicsdb_json_processor.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
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
 * GenomicsDB json processing payloads
 *
 **/
#include "genomicsdb.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "rapidjson/prettywriter.h"

#define TO_JSON_DOCUMENT(X) (reinterpret_cast<rapidjson::Document *>(X))

#define STRING_FIELD(NAME, TYPE) (TYPE.is_string() || TYPE.is_char() || TYPE.num_elements > 1 || (NAME.compare("GT") == 0))
#define INT_FIELD(TYPE) (TYPE.is_int())
#define FLOAT_FIELD(TYPE) (TYPE.is_float())

JSONVariantCallProcessor::JSONVariantCallProcessor(payload_t payload_mode)
    : m_payload_mode(payload_mode) {
  m_json_document = new rapidjson::Document();
  TO_JSON_DOCUMENT(m_json_document)->SetArray();
  m_json_buffer = new rapidjson::StringBuffer();
}

JSONVariantCallProcessor::~JSONVariantCallProcessor() {
  delete reinterpret_cast<rapidjson::StringBuffer *>(m_json_buffer);
  delete TO_JSON_DOCUMENT(m_json_document);
}

const char *JSONVariantCallProcessor::construct_json_output() {
  rapidjson::Document *json_doc = TO_JSON_DOCUMENT(m_json_document);
  json_doc->Clear();
  
  if (m_payload_mode == just_ncalls) {
    json_doc->SetObject();
    json_doc->AddMember("num_calls", m_num_calls, json_doc->GetAllocator());
  } else {
    json_doc->SetObject();
    auto& allocator = json_doc->GetAllocator();
    for (auto const& [sample_name, fields] : m_sample_info) {
      auto i = 0u;
      rapidjson::Value sample_value = rapidjson::Value(rapidjson::kObjectType);
      rapidjson::Value sample_name_value(rapidjson::StringRef(sample_name.c_str()), allocator);
      if (m_payload_mode == samples_with_ncalls) {
        auto chrom_val = reinterpret_cast<rapidjson::Value *>(fields[0]);
        json_doc->AddMember(sample_name_value, chrom_val->Size(), allocator);
      } else {
        rapidjson::Value field_values = rapidjson::Value(rapidjson::kObjectType);
        field_values.AddMember("CHROM",
                               *(reinterpret_cast<rapidjson::Value *>(fields[i++])),
                               allocator);
        field_values.AddMember("POS",
                               *(reinterpret_cast<rapidjson::Value *>(fields[i++])),
                               allocator);
      
        for (auto field_name : m_field_names) {
          rapidjson::Value field_name_value(rapidjson::StringRef(field_name.c_str()), allocator);
          field_values.AddMember(field_name_value,
                                 *(reinterpret_cast<rapidjson::Value *>(fields[i++])),
                                 allocator);
        }
        json_doc->AddMember(sample_name_value, field_values, allocator);
      }
    }
  }

  
  rapidjson::StringBuffer *buffer = reinterpret_cast<rapidjson::StringBuffer *>(m_json_buffer);
  buffer->Clear();
  rapidjson::Writer<rapidjson::StringBuffer> writer(*buffer);
  json_doc->Accept(writer);
  
  return buffer->GetString();
}

void JSONVariantCallProcessor::process(const interval_t& interval) {
  if (!m_is_initialized) {
    m_is_initialized = true;
    auto& genomic_field_types = get_genomic_field_types();
    for (auto& field_type_pair : *genomic_field_types) {
      std::string field_name = field_type_pair.first;
      genomic_field_type_t field_type = field_type_pair.second;
      if (field_name.compare("END") && field_name.compare("REF") && field_name.compare("ALT")) {
        m_field_names.push_back(field_name);
      }
    }
  }
}

void JSONVariantCallProcessor::process(const std::string& sample_name,
                                       const int64_t* coordinates,
                                       const genomic_interval_t& genomic_interval,
                                       const std::vector<genomic_field_t>& genomic_fields) {
  m_num_calls++;
  if (m_payload_mode == just_ncalls) return;

  rapidjson::Document *json_doc = TO_JSON_DOCUMENT(m_json_document);
  auto& allocator = json_doc->GetAllocator();
  auto sample_info = m_sample_info.find(sample_name);
  if (sample_info == m_sample_info.end()) {
    // Create Value Arrays and add to map
    std::vector<void *> fields;
    rapidjson::Value *chrom = new rapidjson::Value(rapidjson::kArrayType);
    fields.push_back(chrom);
    if (m_payload_mode == all) {
      rapidjson::Value *pos = new rapidjson::Value(rapidjson::kArrayType);
      fields.push_back(pos);
      for (auto field_name : m_field_names) {
        rapidjson::Value *field_array = new rapidjson::Value(rapidjson::kArrayType);
        fields.push_back(field_array);
      }
    }
    m_sample_info.emplace(std::make_pair(sample_name, std::move(fields)));
    sample_info = m_sample_info.find(sample_name);
  }

  auto i = 0u;
  auto& fields = sample_info->second;

  rapidjson::Value *chrom = reinterpret_cast<rapidjson::Value *>(fields[i++]);
  rapidjson::Value chrom_value = rapidjson::Value(rapidjson::kStringType);
  chrom_value.SetString(genomic_interval.contig_name.c_str(),
                        genomic_interval.contig_name.length(),
                        allocator);
  chrom->PushBack(chrom_value, allocator);

  if (m_payload_mode != all) return;

  rapidjson::Value *pos = reinterpret_cast<rapidjson::Value *>(fields[i++]);
  pos->PushBack(rapidjson::Value(genomic_interval.interval.first).Move(), allocator);

  for (auto field_name : m_field_names) {
    auto field_type = get_genomic_field_types()->at(field_name);
    
    bool found = false;
    for (auto genomic_field: genomic_fields) {
      if (genomic_field.name.compare(field_name) == 0) {
        rapidjson::Value *values = reinterpret_cast<rapidjson::Value *>(fields[i++]);
        if (STRING_FIELD(field_name, field_type)) {
          rapidjson::Value value = rapidjson::Value(rapidjson::kStringType);
          std::string str;
          if (field_name == "GT") {
            str = resolve_gt(genomic_fields);
          } else {
            str = genomic_field.to_string(field_type);
          }
          value.SetString(str.c_str(), str.length(), allocator);
          values->PushBack(value, allocator);
        } else if (INT_FIELD(field_type)) {
          values->PushBack(rapidjson::Value(genomic_field.int_value_at(0)).Move(), allocator);
        } else if (FLOAT_FIELD(field_type)) {
          values->PushBack(rapidjson::Value(genomic_field.float_value_at(0)).Move(), allocator);
        } else {
          std::string msg = "Genomic field type for " + field_name + " not supported";
          throw GenomicsDBException(msg.c_str());
        }
        found = true;
        break;
      }
    }

    if (!found) {
      rapidjson::Value *values = reinterpret_cast<rapidjson::Value *>(fields[i++]);
      rapidjson::Value value = rapidjson::Value(rapidjson::kStringType);
      value.SetString(".", allocator);
      values->PushBack(value, allocator);
    }
  }
}
