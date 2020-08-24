/**
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation, Inc.
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
 **/

#include <cstring>
#include <iostream>
#include <string>

#include <google/protobuf/util/json_util.h>

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"

#include "genomicsdb_import_config.pb.h"
#include "tiledb_utils.h"
#include "vid_mapper_pb.h"

using namespace genomicsdb_pb;

static GenomicsDBFieldInfo* get_field_info(VidMappingPB* vid_map_protobuf, const std::string& field_name) {
  GenomicsDBFieldInfo* field_info = NULL;
  for (int i=0; i<vid_map_protobuf->fields_size(); i++) {
     field_info = vid_map_protobuf->mutable_fields(i);
     if (!field_info->name().compare(field_name)) {
       return field_info;
     }
  }
  field_info = vid_map_protobuf->add_fields();
  field_info->set_name(field_name);
  return field_info;
}

static int add_filter_fields(VidMappingPB* vid_map_protobuf, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    GenomicsDBFieldInfo* field_info =  get_field_info(vid_map_protobuf, (hrec->vals[index]));
    field_info->set_name(hrec->vals[index]);
    field_info->add_vcf_field_class()->assign("FILTER");
    field_info->add_type()->assign("Integer");
    field_info->add_length()->set_variable_length_descriptor("1");
  }
  return 0;
}

static int add_fmt_fields(VidMappingPB* vid_map_protobuf, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    GenomicsDBFieldInfo* field_info = get_field_info(vid_map_protobuf, hrec->vals[index]);
    field_info->add_vcf_field_class()->assign("FORMAT");
    if (!strcmp(hrec->vals[0], "GT")) {
      field_info->add_type()->assign("int");
      field_info->add_length()->set_variable_length_descriptor("PP");
    } else {
      index = bcf_hrec_find_key(hrec, "Type");
      if (index >= 0) { 
        field_info->add_type()->assign(hrec->vals[index]);
      }
      index = bcf_hrec_find_key(hrec, "Number");
      if (index >= 0) {
        field_info->add_length()->set_variable_length_descriptor(hrec->vals[index]);
      } else {        
        field_info->add_length()->set_variable_length_descriptor("var");
      }
    }
  }
  return 0;
}

static int add_info_fields(VidMappingPB* vid_map_protobuf, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    GenomicsDBFieldInfo* field_info = get_field_info(vid_map_protobuf, hrec->vals[index]);
    field_info->add_vcf_field_class()->assign("INFO");
    field_info->add_type()->assign("int");
    field_info->add_length()->set_variable_length_descriptor("1");
  }
  return 0;
}

static int add_contigs(VidMappingPB* vid_map_protobuf, ImportConfiguration* import_config_protobuf, bcf_hrec_t* hrec, int64_t* tiledb_column_offset) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    Chromosome* chromosome = vid_map_protobuf->add_contigs();
    chromosome->set_name(hrec->vals[index]);
    index = bcf_hrec_find_key(hrec, "length");
    int64_t length = std::stoll(hrec->vals[index]);
    chromosome->set_length(length);
    chromosome->set_tiledb_column_offset(*tiledb_column_offset);

    Partition* partition = import_config_protobuf->add_column_partitions();
    GenomicsDBColumn* genomicsdb_column = new GenomicsDBColumn();
    genomicsdb_column->set_tiledb_column(*tiledb_column_offset);
    partition->set_allocated_begin(genomicsdb_column);
    partition->set_array_name(chromosome->name());
    
    *tiledb_column_offset = *tiledb_column_offset + length;
  }
  return 0;
}

static int generate_json(const std::string& merged_header, const std::string& vidmap_output) {
  htsFile* fptr = hts_open(merged_header.c_str(), "r");
  if (!fptr) {
    return -1;
  }

  bcf_hdr_t* hdr = bcf_hdr_read(fptr);
  if (!hdr) {
    return -1;
  }
  if (bcf_hdr_set_samples(hdr, NULL, 1)) {
    return -1;
  }

  ImportConfiguration* import_config_protobuf = new ImportConfiguration();
  import_config_protobuf->set_vid_mapping_file(vidmap_output);
  import_config_protobuf->set_callset_mapping_file("callset.json");

  VidMappingPB* vid_map_protobuf = new VidMappingPB();
  int64_t tiledb_column_offset = 0;
  
  for (int i=0; i<hdr->nhrec; i++) {
    bcf_hrec_t* hrec = hdr->hrec[i];
    switch (hrec->type) {
      case BCF_HL_FLT:
        add_filter_fields(vid_map_protobuf, hrec);
        break;
      case BCF_HL_FMT:
        add_fmt_fields(vid_map_protobuf, hrec);
        break;
      case BCF_HL_INFO:
        add_info_fields(vid_map_protobuf, hrec);
        break;
      case BCF_HL_CTG:
        add_contigs(vid_map_protobuf, import_config_protobuf, hrec, &tiledb_column_offset);
        break;
      default:
        break;
    }
  }

  bcf_hdr_destroy(hdr);
  hts_close(fptr);

  // Write out vidmap.json
  std::string json;
  google::protobuf::util::JsonPrintOptions json_options;
  json_options.add_whitespace = false;
  json_options.always_print_primitive_fields = false;
  google::protobuf::util::MessageToJsonString(*vid_map_protobuf, &json, json_options);
  TileDBUtils::delete_file(vidmap_output);
  TileDBUtils::write_file(vidmap_output, json.data(), json.length());
  delete vid_map_protobuf;

  // Write out loader.json
  json.clear();
  google::protobuf::util::MessageToJsonString(*import_config_protobuf, &json, json_options);
  TileDBUtils::delete_file("loader.json");
  TileDBUtils::write_file("loader.json", json.data(), json.length());
  delete import_config_protobuf;

  return 0;
}

static int merge_headers(const std::string& sample_filenames, const char* delim, const std::string& merged_header, const std::string& callset_output) {
  htsFile* merged_header_fptr = hts_open(merged_header.c_str(), "w");
  if (!merged_header_fptr) {
    return -1;
  }

  bcf_hdr_t* merged_hdr = bcf_hdr_init("w");

  int64_t row_index = 0;
  CallsetMappingPB* callset_protobuf = new CallsetMappingPB();
  char *sample_filename = std::strtok(const_cast<char *>(sample_filenames.c_str()), delim);
  while (sample_filename) {
    htsFile* fptr = hts_open(sample_filename, "r");
    if (!fptr) {
      return -1;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fptr);
    if (!hdr) {
      return -1;
    }
    for (auto i=0; i<bcf_hdr_nsamples(hdr); i++) {
      auto* callset = callset_protobuf->add_callsets();
      callset->set_sample_name(hdr->samples[i]);
      callset->set_row_idx(row_index);
      callset->set_idx_in_file(i);
      callset->set_filename(sample_filename);
    }

    // Don't care about samples from here on for creating vcf header file
    if (bcf_hdr_set_samples(hdr, NULL, 1)) {
      return -1;
    }

    if (bcf_hdr_write(merged_header_fptr, hdr)) {
      std::cout << "Error merging headers\n";
    }

    bcf_hdr_destroy(hdr);
    hts_close(fptr);

    sample_filename = std::strtok(NULL, delim);
  }

  if (merged_hdr) {
    bcf_hdr_destroy(merged_hdr);
  }
  
  hts_close(merged_header_fptr);

  // Write out callset.json
  std::string json;
  google::protobuf::util::JsonPrintOptions json_options;
  json_options.add_whitespace = false;
  json_options.always_print_primitive_fields = false;
  google::protobuf::util::MessageToJsonString(*callset_protobuf, &json, json_options);

  TileDBUtils::delete_file(callset_output);
  TileDBUtils::write_file(callset_output, json.data(), json.length());

  delete callset_protobuf;  

  return 0;
}

int main(int argc, char** argv) {
  // TODO: Evolve CLI for these input args
  // TODO: Hardcoded for testing
  std::string vidmap_output = "vidmap.json";
  std::string callset_output = "callset.json";
  std::string loader_json = "loader.json";
  // List of Filenames delimited with ";"
  std::string sample_filenames = "sample1.vcf;sample1-copy.vcf";
  const char* delim = " ;";
  std::string merged_header = "vcf_header.vcf";

  if (merge_headers(sample_filenames, delim, merged_header, callset_output)) {
    return -1;
  }

  if (generate_json(merged_header, vidmap_output)) {
    return -1;
  }

  google::protobuf::ShutdownProtobufLibrary();
  std::cout << "Success!\n";
  
  return 0;
}
