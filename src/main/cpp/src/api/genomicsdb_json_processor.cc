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
#include "variant_query_config.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "rapidjson/prettywriter.h"

#define TO_JSON_DOCUMENT(X) (reinterpret_cast<rapidjson::Document *>(X))

// Prototypes to internal methods in genomicsdb.cc declared here instead of header to keep the api opaque
std::map<std::string, genomic_field_type_t> create_genomic_field_types(const VariantQueryConfig &query_config,
                                                   void *annotation_service, bool change_alt_to_string=false);
JSONVariantCallProcessor::JSONVariantCallProcessor(payload_t payload_mode)
    : m_payload_mode(payload_mode) {
  m_json_document = new rapidjson::Document();
  TO_JSON_DOCUMENT(m_json_document)->SetObject();
}

JSONVariantCallProcessor::~JSONVariantCallProcessor() {
  if (m_samples_info)  {
    for (auto const& [sample_name, fields] : *m_samples_info.get()) {
      for (auto const& fields: fields) {
        delete reinterpret_cast<rapidjson::Value *>(fields);
      }
    }
  }
  delete TO_JSON_DOCUMENT(m_json_document);
}


void JSONVariantCallProcessor::initialize(const VariantQueryConfig &query_config,
                                          void *annotation_service) {
  GenomicsDBVariantCallProcessor::initialize(query_config, annotation_service);
  if (!m_is_initialized) {
    m_is_initialized = true;
    switch (m_payload_mode) {
      case just_samples:
        m_samples_set = std::make_unique<std::set<std::string>>();
      case just_ncalls:
        return;
      default:
        ;
    }
    for (auto i=0u; i<query_config.get_num_queried_attributes(); i++) {
      const std::string attribute_name = query_config.get_query_attribute_name(i);
      if (attribute_name.compare("END") && attribute_name.compare("ALT")) {
        m_field_names.push_back(attribute_name);
      }
    }
    switch (m_payload_mode) {
      case samples_with_ncalls:
        m_samples = std::make_unique<std::map<std::string, int64_t>>();
        break;
      case all_by_calls: case all:
        m_samples_info = std::make_unique<std::map<std::string, std::vector<void *>>>();
      default:
        ;
    }
  }
}

std::string JSONVariantCallProcessor::construct_json_output() {
  rapidjson::Document *json_doc = TO_JSON_DOCUMENT(m_json_document);
  auto& allocator = json_doc->GetAllocator();
  switch(m_payload_mode) {
    case just_ncalls: {
      json_doc->AddMember("num_calls", m_num_calls, allocator);
      break;
    }
    case just_samples: {
      json_doc->SetArray();
      if (m_samples_set == nullptr || m_samples_set->size() == 0) break;
      for (auto const& sample : *m_samples_set.get()) {
        rapidjson::Value sample_value(rapidjson::StringRef(sample.c_str()), allocator);
        json_doc->PushBack(sample_value, allocator);
      }
      break;
    }
    case samples_with_ncalls: {
      if (m_samples == nullptr || m_samples->size() == 0) break;
      for (auto const& [sample_name, num_calls] : *m_samples.get()) {
        rapidjson::Value sample_name_value(rapidjson::StringRef(sample_name.c_str()), allocator);
        json_doc->AddMember(sample_name_value, rapidjson::Value(num_calls).Move(), allocator);
      }
      break;
    }
    case all_by_calls: {
      if (m_samples_info == nullptr || m_samples_info->size() == 0) break;
      // Create header with names of the fields in the JSON output
      rapidjson::Value field_names = rapidjson::Value(rapidjson::kArrayType);
      field_names.PushBack("CHR", allocator);
      field_names.PushBack("POS", allocator);
      for (auto field_name : m_field_names) {
        rapidjson::Value value(rapidjson::StringRef(field_name.c_str()), allocator);
        field_names.PushBack(value, allocator);
      }
      json_doc->AddMember("FIELD", field_names, allocator);
      for (auto const& [sample_name, fields_vec] : *m_samples_info.get()) {
        rapidjson::Value fields_array = rapidjson::Value(rapidjson::kArrayType);
        for (auto const& fields: fields_vec) {
          fields_array.PushBack(*(reinterpret_cast<rapidjson::Value *>(fields)), allocator);
        }
        rapidjson::Value sample_value = rapidjson::Value(rapidjson::kObjectType);
        rapidjson::Value sample_name_value(rapidjson::StringRef(sample_name.c_str()), allocator);
        json_doc->AddMember(sample_name_value, fields_array, allocator);
      }
      break;
    }
    case all: {
      if (m_samples_info == nullptr || m_samples_info->size() == 0) break;
      for (auto const& [sample_name, fields] : *m_samples_info.get()) {
        auto i = 0u;
        rapidjson::Value field_values = rapidjson::Value(rapidjson::kObjectType);
        field_values.AddMember("CHR",
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
        rapidjson::Value sample_name_value(rapidjson::StringRef(sample_name.c_str()), allocator);
        json_doc->AddMember(sample_name_value, field_values, allocator);
      }
    }
  }

  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  json_doc->Accept(writer);
  
  return std::string(buffer.GetString(), buffer.GetLength());
}

void JSONVariantCallProcessor::process(const interval_t& interval) {
}

static rapidjson::Value get_chr_value(const genomic_interval_t& genomic_interval,
                                      rapidjson::Document::AllocatorType& allocator) {
  rapidjson::Value value = rapidjson::Value(rapidjson::kStringType);
  value.SetString(genomic_interval.contig_name.c_str(),
                  genomic_interval.contig_name.length(),
                  allocator);
  return value;
}

static rapidjson::Value get_value(const std::string& field_name,
                            genomic_field_type_t field_type,
                            const std::vector<genomic_field_t>& genomic_fields,
                            size_t genomic_idx,
                            rapidjson::Document::AllocatorType& allocator) {
  if (STRING_FIELD(field_name, field_type)) {
    std::string str;
    if (field_name == "GT") {
      str = resolve_gt(genomic_fields);
    } else {
      str = genomic_fields[genomic_idx].to_string(field_type);
    }
    rapidjson::Value value(rapidjson::kStringType);
    value.SetString(str.c_str(), str.length(), allocator);
    return value;
  } else if (INT_FIELD(field_type)) {
    return rapidjson::Value(genomic_fields[genomic_idx].int_value_at(0));
  } else if (FLOAT_FIELD(field_type)) {
    return rapidjson::Value(genomic_fields[genomic_idx].float_value_at(0));
  }
  return rapidjson::Value(rapidjson::kNullType);
}

void JSONVariantCallProcessor::process(const std::string& sample_name,
                                       const int64_t* coordinates,
                                       const genomic_interval_t& genomic_interval,
                                       const std::vector<genomic_field_t>& genomic_fields) {
  auto& allocator = TO_JSON_DOCUMENT(m_json_document)->GetAllocator();
  switch (m_payload_mode) {
    case just_ncalls:
      m_num_calls++;
      break;
    case just_samples:
      m_samples_set->insert(sample_name);
      break;
    case samples_with_ncalls: {
      auto sample = m_samples->find(sample_name);
      if (sample == m_samples->end()) {
        m_samples->emplace(sample_name, 0ul);
        sample = m_samples->find(sample_name);
      }
      m_samples->emplace(sample_name, sample->second++);
      break;
    }
    case all_by_calls: {
      auto sample_info = m_samples_info->find(sample_name);
      if (sample_info == m_samples_info->end()) {
        // Create Value Arrays and add to map
        m_samples_info->emplace(sample_name, std::vector<void *>());
        sample_info = m_samples_info->find(sample_name);
      }
      rapidjson::Value *fields = new rapidjson::Value(rapidjson::kArrayType);
      // CHR value
      fields->PushBack(get_chr_value(genomic_interval, allocator).Move(), allocator);
      // POS value
      fields->PushBack(rapidjson::Value(genomic_interval.interval.first).Move(), allocator);
      auto idx = 0u;
      for (auto field_name : m_field_names) {
        if (idx < genomic_fields.size() && genomic_fields[idx].name.compare("ALT") == 0) {
          idx++;
        }
        if (idx < genomic_fields.size() && genomic_fields[idx].name.compare(field_name) == 0) {
          auto field_type = get_genomic_field_types()->at(field_name);
          fields->PushBack(get_value(field_name, field_type, genomic_fields, idx, allocator).Move(), allocator);
          idx++;
        } else {
          fields->PushBack(rapidjson::Value(rapidjson::kNullType), allocator);
        }
      }
      sample_info->second.push_back(fields);
      break;
    }
    case all: {
      auto sample_info = m_samples_info->find(sample_name);
      if (sample_info == m_samples_info->end()) {
        // Create Value Arrays and add to map
        std::vector<void *> fields;
        rapidjson::Value *chr = new rapidjson::Value(rapidjson::kArrayType);
        fields.push_back(chr);
        rapidjson::Value *pos = new rapidjson::Value(rapidjson::kArrayType);
        fields.push_back(pos);
        for (auto field_name : m_field_names) {
          rapidjson::Value *field_array = new rapidjson::Value(rapidjson::kArrayType);
          fields.push_back(field_array);
        }
        m_samples_info->emplace(sample_name, std::move(fields));
        sample_info = m_samples_info->find(sample_name);
      }
      auto i = 0u;
      auto& fields = sample_info->second;
      rapidjson::Value *chr = reinterpret_cast<rapidjson::Value *>(fields[i++]);
      chr->PushBack(get_chr_value(genomic_interval, allocator).Move(), allocator);
      rapidjson::Value *pos = reinterpret_cast<rapidjson::Value *>(fields[i++]);
      pos->PushBack(rapidjson::Value(genomic_interval.interval.first).Move(), allocator);
      auto idx = 0u;
      for (auto field_name : m_field_names) {
        rapidjson::Value *values = reinterpret_cast<rapidjson::Value *>(fields[i++]);
        if (idx < genomic_fields.size() && genomic_fields[idx].name.compare("ALT") == 0) {
          idx++;
        }
        if (idx < genomic_fields.size() && genomic_fields[idx].name.compare(field_name) == 0) {
          auto field_type = get_genomic_field_types()->at(field_name);
          values->PushBack(get_value(field_name, field_type, genomic_fields, idx, allocator).Move(), allocator);
          idx++;
        } else {
          values->PushBack(rapidjson::Value(rapidjson::kNullType), allocator);
        }
      }
    }
  }
}
