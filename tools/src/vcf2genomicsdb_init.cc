/**
 * The MIT License (MIT)
 * Copyright (c) 2020-2021 Omics Data Automation, Inc.
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
#include <errno.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <regex>
#include <sstream>
#include <string>

#include <google/protobuf/util/json_util.h>

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"

#include "common.h"
#include "genomicsdb_config_base.h"
#include "genomicsdb_import_config.pb.h"
#include "hfile_genomicsdb.h"
#include "vid_mapper_pb.h"

#include <spdlog/sinks/stdout_color_sinks.h>

using namespace genomicsdb_pb;

static Logger g_logger(Logger::get_logger("vcf2genomicsdb_init"));

struct partition_config_t {
  std::string interval_list;
  int64_t number_of_partitions = 0;
  int64_t size_of_partitions = 0;
  bool merge_small_contigs = false;
};
partition_config_t g_partition_config;

static GenomicsDBFieldInfo* get_field_info(VidMappingPB* vidmap_pb, const std::string& field_name) {
  GenomicsDBFieldInfo* field_info = NULL;
  for (int i=0; i<vidmap_pb->fields_size(); i++) {
     field_info = vidmap_pb->mutable_fields(i);
     if (!field_info->name().compare(field_name)) {
       return field_info;
     }
  }
  field_info = vidmap_pb->add_fields();
  field_info->set_name(field_name);
  return field_info;
}

static int add_filter_fields(VidMappingPB* vidmap_pb, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    GenomicsDBFieldInfo* field_info =  get_field_info(vidmap_pb, (hrec->vals[index]));
    field_info->set_name(hrec->vals[index]);
    field_info->add_vcf_field_class()->assign("FILTER");
    field_info->add_type()->assign("Integer");
    field_info->add_length()->set_variable_length_descriptor("1");
  }
  return 0;
}

static int add_fields(GenomicsDBFieldInfo* field_info, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "Type");
  std::string field_type = hrec->vals[index];
  field_info->add_type()->assign(field_type);
  if (!field_type.compare("String")) {
    field_info->add_length()->set_variable_length_descriptor("var");
  } else if (!field_type.compare("Flag")) {
    field_info->add_length()->set_variable_length_descriptor("1");
  } else {
    index = bcf_hrec_find_key(hrec, "Number");
    std::string field_number = hrec->vals[index];
    if (index >= 0 && field_number.compare(".")) {
      field_info->add_length()->set_variable_length_descriptor(hrec->vals[index]);
    } else {
      field_info->add_length()->set_variable_length_descriptor("var");
    } 
  }
  return 0;
}

static int add_fmt_fields(VidMappingPB* vidmap_pb, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    GenomicsDBFieldInfo* field_info = get_field_info(vidmap_pb, hrec->vals[index]);
    field_info->add_vcf_field_class()->assign("FORMAT");
    if (!strcmp(hrec->vals[index], "GT")) {
      field_info->add_type()->assign("Integer");
      field_info->add_length()->set_variable_length_descriptor("PP");
    } else {
      add_fields(field_info, hrec);
    }
  }
  return 0;
}

static int add_info_fields(VidMappingPB* vidmap_pb, bcf_hrec_t* hrec) {
  int index = bcf_hrec_find_key(hrec, "ID");
  if (index >= 0) {
    GenomicsDBFieldInfo* field_info = get_field_info(vidmap_pb, hrec->vals[index]);
    field_info->add_vcf_field_class()->assign("INFO");
    add_fields(field_info, hrec);
    // Specify custom combinations - hardcoded for now
    if (field_info->name() == "END" || field_info->name() == "DP") {
      field_info->set_vcf_field_combine_operation("sum");
    }
  }
  return 0;
}

static int64_t total_number_of_partitions = 0;
static void create_partition(const std::string& contig, int64_t start, int64_t end, int64_t tiledb_column_offset,
                            const std::string& workspace, ImportConfiguration* import_config_protobuf) {
  GenomicsDBColumn* genomicsdb_column_begin = new GenomicsDBColumn();
  genomicsdb_column_begin->set_tiledb_column(tiledb_column_offset+start);
  GenomicsDBColumn* genomicsdb_column_end = new GenomicsDBColumn();
  genomicsdb_column_end->set_tiledb_column(tiledb_column_offset+end-1);
  
  g_logger.debug("Processed length for tiledb_offset={} chr={} start={} end={}", tiledb_column_offset,
                 contig, genomicsdb_column_begin->tiledb_column(), genomicsdb_column_end->tiledb_column());
  
  Partition* partition = import_config_protobuf->add_column_partitions();
  partition->set_workspace(workspace);
  partition->set_allocated_begin(genomicsdb_column_begin);
  partition->set_allocated_end(genomicsdb_column_end);
  if (contig.find("scaffold") == std::string::npos) {
    partition->set_array_name(g_logger.format("{}${}${}", contig, start, end));
  } else {
    partition->set_array_name(contig);
  }
  total_number_of_partitions++;
}

static std::pair<int64_t, int64_t> get_contig_positions(const std::string& contig,
                                              std::vector<std::pair<std::string, int64_t>> regions) {
  int64_t tiledb_column_offset = 0;
  for (auto region: regions) {
    if (region.first == contig) {
      return {region.second, tiledb_column_offset};
    }
    tiledb_column_offset += region.second;
  }
  return {-1, -1};
}
  
static int process_from_interval_list(const std::string& interval_list,
                                      std::vector<std::pair<std::string, int64_t>> regions,
                                      const std::string& workspace, VidMappingPB* vidmap_pb,
                                      ImportConfiguration* import_config_protobuf) {
  // Intervals is a vector of contigs with start and end positions. The std::string points to the contig and the
  // int64_t pair point to the start and end positions.
  std::vector<std::pair<std::string, std::pair<int64_t, int64_t>>> intervals;
  std::ifstream interval_list_stream(interval_list);
  std::string line;
  int64_t tiledb_column_offset = 0;
  while (std::getline(interval_list_stream, line)) {
    std::istringstream interval_line(line);
    std::string interval;
    if (!(interval_line >> interval)) {
      continue;
    }
    // Parse interval
    auto colon = interval.find(":");
    auto hyphen = interval.find("-");
    auto len = interval.length();
    std::string contig;
    int64_t start = 0;
    int64_t end = -1;
    if (colon == std::string::npos && hyphen == std::string::npos) {
      // Entire chromosome/contig
      contig = interval;
    } else if (hyphen == std::string::npos) {
      // Chromosome with start position
      contig = interval.substr(0, colon);
      start = std::stoll(interval.substr(colon+1));
    } else {
      // Interval with start and end positions
      contig = interval.substr(0, colon);
      start = std::stoll(interval.substr(colon+1, hyphen-colon));
      end = std::stoll(interval.substr(hyphen+1));
    }

    g_logger.debug("Processing interval {} parsed as contig={} start={} end={}", interval, contig, start, end);
    auto length_and_tiledb_column_offset = get_contig_positions(contig, regions);
    auto length = length_and_tiledb_column_offset.first;
    auto tiledb_column_offset = length_and_tiledb_column_offset.second;
    if (length == -1 || tiledb_column_offset == -1) {
      g_logger.error("Discarding interval {} from {}, not locatable in the sample vcfs provided",
                     interval, interval_list);
      continue;
    }

    end = end==-1?length:end;
    if ((end-start) < 0 || (end-start) > length) {
       g_logger.error("Discarding interval {} from {} as the start/end positions are out-of-bounds for contig {}",
                      interval, interval_list, contig);
      continue;
    }
    create_partition(contig, start, end, length_and_tiledb_column_offset.second, workspace, import_config_protobuf);
  }
  return 0;
}

static int add_contigs(std::vector<std::pair<std::string, int64_t>> regions, const std::string& workspace,
                       VidMappingPB* vidmap_pb, ImportConfiguration* import_config_protobuf) {
  int64_t total_length = 0;
  for (auto region: regions) {
    Chromosome* chromosome = vidmap_pb->add_contigs();
    chromosome->set_name(region.first);
    chromosome->set_length(region.second);
    chromosome->set_tiledb_column_offset(total_length);
    total_length += region.second;
  }

  if (!g_partition_config.interval_list.empty()) {
    return process_from_interval_list(g_partition_config.interval_list, regions, workspace, vidmap_pb,
                                      import_config_protobuf);
  }

  int64_t number_of_partitions = g_partition_config.number_of_partitions;
  int64_t max_size_of_partition = g_partition_config.size_of_partitions;
  bool merge_small_contigs = g_partition_config.merge_small_contigs;

  if (number_of_partitions) {
    max_size_of_partition = total_length/number_of_partitions;
  }
  g_logger.debug("max_size_of_partition={}", max_size_of_partition);
  int64_t tolerance = max_size_of_partition*0.25;
  int64_t merge_small_contigs_size = 1024*1024; // 1M
  bool partition_by_contig = !number_of_partitions && !max_size_of_partition;

  int nscaffold = 0;
  bool processing_scaffold = false;
  int64_t scaffold_length = 0;

  int64_t tiledb_column_offset = 0;
  for (auto region: regions) {
    std::string chromosome = region.first;
    int64_t length = region.second;
    int64_t processed_length = 0;
    int64_t start, end = 0;

    g_logger.debug("Region for chromosome={} length={}", chromosome, length);
    if (partition_by_contig) {
      max_size_of_partition = std::max(max_size_of_partition, length);
    }

    while (processed_length < length) {
      start = processed_length;
      end = start + ((length-processed_length)<(max_size_of_partition+tolerance)?
                     length-processed_length:max_size_of_partition);

      // Process a scaffold
      if (merge_small_contigs && length < merge_small_contigs_size) {
        if (!processing_scaffold) {
          processing_scaffold = true;
          nscaffold++;
        }
        scaffold_length += length;
        if (scaffold_length > max_size_of_partition) {
          create_partition("scaffold"+std::to_string(nscaffold), 0, scaffold_length, tiledb_column_offset,
                         workspace, import_config_protobuf);
          processing_scaffold = false;
          length = scaffold_length;
          scaffold_length = 0;
        }
        break;
      } else if (processing_scaffold) { // New region that is larger than a scaffold
        // Write out the scaffold before proccessing the new region
        create_partition("scaffold"+std::to_string(nscaffold), 0, scaffold_length, tiledb_column_offset,
                         workspace, import_config_protobuf);
        processing_scaffold = false;
        tiledb_column_offset += scaffold_length;
      }

      create_partition(chromosome, start, end, tiledb_column_offset, workspace, import_config_protobuf);        
      processed_length += (end - start);
    }
    //  assert(processed_length == length);
    if (!processing_scaffold) {
      tiledb_column_offset += length;
    }
  }
  
  return 0;
}

#define QUOT "\""
static void fix_protobuf_generated_json_for_int64(std::string& json, const std::vector<std::string>& names) {
  std::string str;
  for (auto name: names) {
    str.empty()?(str=QUOT+name+QUOT):(str=str+"|"+QUOT+name+QUOT);
  }
  //        $1          $2             $3
  str = "("+str+")" + "(:\\s)" + QUOT+"(.*)"+QUOT;
  std::regex pattern(str);
  while (std::regex_search(json, pattern)) {
    json = std::regex_replace(json, pattern, "$1$2$3");
  }
}

static int generate_json(const std::string& merged_header, const std::string& workspace,
                         const std::string& vidmap_output, const std::string& loader_output,
                         const std::string& callset_output) {
  htsFile* fptr = hts_open(merged_header.c_str(), "r");
  if (!fptr) {
    g_logger.error("Could not hts_open {} file in read mode {}", merged_header, strerror(errno));
    return ERR;
  }

  bcf_hdr_t* hdr = bcf_hdr_read(fptr);
  if (!hdr) {
    g_logger.error("Could not read header file {}", merged_header);
    return ERR;
  }

  // Only need the header, so can set samples to NULL to bypass optimization
  if (bcf_hdr_set_samples(hdr, NULL, 1)) {
    g_logger.error("Error in bcf_hdr_set_samples");
    return ERR;
  }

  ImportConfiguration* import_config_protobuf = new ImportConfiguration();
  import_config_protobuf->set_vid_mapping_file(vidmap_output);
  import_config_protobuf->set_callset_mapping_file(callset_output);
  import_config_protobuf->set_size_per_column_partition(43581440);

  // These are defaults, listing them out in the loader.json for easier editing
  // Can make these configurable from the command line as needed
  import_config_protobuf->set_segment_size(1048576);
  import_config_protobuf->set_num_cells_per_tile(1000);
  import_config_protobuf->set_consolidate_tiledb_array_after_load(false);
  import_config_protobuf->set_enable_shared_posixfs_optimizations(false);

  VidMappingPB* vidmap_pb = new VidMappingPB();
  // Using vector instead of map here because map does not keep track of order of insertion as it uses some 
  // form of binary tree for fast retrieval
  std::vector<std::pair<std::string, int64_t>> regions;
  int64_t total_length = 0;
  
  for (int i=0; i<hdr->nhrec; i++) {
    bcf_hrec_t* hrec = hdr->hrec[i];
    switch (hrec->type) {
      case BCF_HL_FLT:
        add_filter_fields(vidmap_pb, hrec);
        break;
      case BCF_HL_FMT:
        add_fmt_fields(vidmap_pb, hrec);
        break;
      case BCF_HL_INFO:
        add_info_fields(vidmap_pb, hrec);
        break;
      case BCF_HL_CTG: {
        int index = bcf_hrec_find_key(hrec, "ID");
        assert(index >= 0);
        std::string id = hrec->vals[index];
        index = bcf_hrec_find_key(hrec, "length");
        assert(index >= 0);
        int64_t length = std::stoll(hrec->vals[index]);
        total_length += length;
        regions.push_back({id, length});
        break;
      }
      default:
        break;
    }
  }

  add_contigs(regions, workspace, vidmap_pb, import_config_protobuf);
  g_logger.info("Total number of partitions in loader.json={}", total_number_of_partitions);
  
  bcf_hdr_destroy(hdr);
  hts_close(fptr);

  std::string json;
  google::protobuf::util::JsonPrintOptions json_options;
  json_options.add_whitespace = true;
  json_options.always_print_primitive_fields = false;
  json_options.preserve_proto_field_names = true;

  // Write out vidmap.json  
  google::protobuf::util::MessageToJsonString(*vidmap_pb, &json, json_options);
  g_logger.debug("Writing out vidmap json file to {}", vidmap_output);
  fix_protobuf_generated_json_for_int64(json, {"tiledb_column_offset", "length"});
  TileDBUtils::write_file(vidmap_output, json.data(), json.length());
  delete vidmap_pb;

  // Write out loader.json
  json.clear();
  google::protobuf::util::MessageToJsonString(*import_config_protobuf, &json, json_options);
  g_logger.debug("Writing out loader json file to {}", loader_output);
  fix_protobuf_generated_json_for_int64(json, {"tiledb_column", "size_per_column_partition", "segment_size",
                                               "num_cells_per_tile"});
  TileDBUtils::write_file(loader_output, json.data(), json.length());
  delete import_config_protobuf;

  return 0;
}

static int merge_headers_and_generate_callset(const std::string& sample_list, const std::string& merged_header, const std::string& callset_output) {
  htsFile* merged_header_fptr = hts_open(merged_header.c_str(), "w");
  if (!merged_header_fptr) {
    g_logger.error("Could not hts_open {} file in write mode {}", merged_header, strerror(errno));
    return ERR;
  }

  bcf_hdr_t* merged_hdr = bcf_hdr_init("w");
  if (!merged_hdr) {
    g_logger.error("Could not initialize using bcf_hdr_init in write mode for {}", merged_header);
    return ERR;
  }

  int64_t row_index = 0;
  CallsetMappingPB* callset_protobuf = new CallsetMappingPB();
  std::ifstream sample_list_stream(sample_list);
  std::string line;
  while (std::getline(sample_list_stream, line)) {
    std::istringstream sample_line(line);
    std::string sample_name, sample_uri;
    if (!(sample_line >> sample_uri)) {
      continue;
    }

    htsFile* fptr = hts_open(sample_uri.c_str(), "r");
    if (!fptr) {
      g_logger.error("Could not open sample {} with hts_open {}", sample_uri, strerror(errno));
      return ERR;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fptr);
    if (!hdr) {
      g_logger.error("Could not read sample {} with htslib", sample_uri);
      return ERR;
    }
    for (auto i=0; i<bcf_hdr_nsamples(hdr); i++) {
      auto* callset = callset_protobuf->add_callsets();
      callset->set_sample_name(hdr->samples[i]);
      callset->set_row_idx(row_index++);
      callset->set_idx_in_file(i);
      callset->set_filename(sample_uri);
    }

    // no need to set samples for just creating a vcf header file, so invoke with samples=NULL
    if (bcf_hdr_set_samples(hdr, NULL, 1)) {
      g_logger.error("Could not set samples for {} with htslib", sample_uri);
      return ERR;
    }

    bcf_hdr_t* resultant_hdr = bcf_hdr_merge(merged_hdr, hdr);
    assert(resultant_hdr = merged_hdr);
    
    bcf_hdr_destroy(hdr);
    hts_close(fptr);
  }

  g_logger.debug("Writing out template header file to {}", merged_header);
  if (bcf_hdr_write(merged_header_fptr, merged_hdr)) {
    g_logger.error("Error while writing out merged header file at {}", merged_header);
    return ERR;
  }
  bcf_hdr_destroy(merged_hdr);
  hts_close(merged_header_fptr);

  // Write out callset.json
  std::string json;
  google::protobuf::util::JsonPrintOptions json_options;
  json_options.add_whitespace = true;
  json_options.always_print_primitive_fields = false;
  json_options.preserve_proto_field_names = true;
  google::protobuf::util::MessageToJsonString(*callset_protobuf, &json, json_options);
  fix_protobuf_generated_json_for_int64(json, {"row_idx", "idx_in_file"});
  
  g_logger.debug("Writing out callset json file to {}", callset_output);
  int rc = TileDBUtils::write_file(callset_output, json.data(), json.length());
  delete callset_protobuf;
  
  return rc;
}

enum GenomicsDBArgsEnum {
  VERSION=1000,
};

void print_usage() {
  std::cerr << "Usage: vcf2genomicsdb_init [options]\n"
            << "where options include:\n"
            << "\t \e[1m--help\e[0m, \e[1m-h\e[0m Print a usage message summarizing options available and exit\n"
            << "\t \e[1m--workspace\e[0m=<GenomicsDB workspace URI>, \e[1m-w\e[0m <GenomicsDB workspace URI>\n"
            << "\t\t If workspace does not exist, it is created first\n"
            << "\t\t exits if workspace exists and is invoked without the overwrite-workspace option\n"
            << "\t \e[1m--overwrite-workspace\e[0m, \e[1m-o\e[0m\n"
            << "\t\t Allow for workspace json artifacts to be overwritten\n"
            << "\t \e[1m--sample-list\e[0m=<sample list>, \e[1m-s\e[0m <sample list file>\n"
            << "\t\t Specify sample URIs for import, one line per sample path\n"
            << "\t \e[1m--interval-list\e[0m=<genomic interval list>, \e[1m-i\e[0m <genomic interval list file>\n"
            << "\t\t Optional, create array partitions from intervals in interval list, one line per interval,\n" 
            << "\t\t default is partition by chromosome/contig, overrides --number-of-array-partitions and \n"
            << "\t\t --size-of-array-partitions\n"
            << "\t \e[1m--number-of-array-partitions\e[0m=<number>, \e[1m-n\e[0m <number>\n"
            << "\t\t Optional, suggested number of array partitions\n"
            << "\t \e[1m--size-of-array-partitions\e[0m=<size>, \e[1m-z\e[0m <size>\n"
            << "\t\t Optional, suggest size of arrays partitions, overrides --number-of-array-partitions\n"
            << "\t \e[1m--merge-small-contigs\e[0m, \e[1m-m\e[0m\n"
            << "\t\t Optional, default is false\n"
            << "\t \e[1m--verbose\e[0m, \e[1m-v\e[0m\n"
            << "\t\t Allow verbose messages to be logged\n"
            << "\t \e[1m--version\e[0m Print version and exit\n"
            << std::endl;
}

#define SLASHIFY(path) (path.back()!='/'?path+'/':path)
#define UNSLASHIFY(path) (path.back()=='/'?path.substr(0,path.length()-2):path)

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage();
    return ERR;
  }
  
  static struct option long_options[] = {
    {"workspace",                  1, 0, 'w'},
    {"overwrite-workspace",        0, 0, 'o'},
    {"sample-list",                1, 0, 's'},
    {"interval-list",              1, 0, 'i'},
    {"number-of-array-partitions", 1, 0, 'n'},
    {"size-of-array-partitions",   1, 0, 'z'},
    {"merge-small-contigs",        0, 0, 'm'},
    {"verbose",                    0, 0, 'v'},
    {"version",                    0, 0, VERSION},
    {"help",                       0, 0, 'h'},
    {0,0,0,0}
  };

  std::string workspace;
  bool overwrite_workspace = false;
  std::string sample_list;
  bool verbose = false;

  int c;
  while ((c=getopt_long(argc, argv, "w:os:i:n:z:mvh", long_options, NULL)) >= 0) {
    switch (c) {
      case 'w':
        workspace = std::move(std::string(optarg));
        break;
      case 'o':
        overwrite_workspace = true;
        break;
      case 's':
        sample_list = std::move(std::string(optarg));
        break;
      case 'i':
        g_partition_config.interval_list = std::move(std::string(optarg));
        break;
      case 'n':
        g_partition_config.number_of_partitions = std::stoi(optarg);
        break;
      case 'z':
        g_partition_config.size_of_partitions = std::stoi(optarg);
        break;
      case 'm':
        g_partition_config.merge_small_contigs = true;
        break;
      case 'v':
        verbose = true;
        break;
      case 'h':
        print_usage();
        return 0;
      case VERSION:
        std::cout << GENOMICSDB_VERSION << "\n";
        return OK;
      default:
        std::cerr << "Unknown command line argument\n";
        print_usage();
        return ERR;
    }
  } 

  if (workspace.empty()) {
    g_logger.error("Workspace required to be specified");
    return ERR;
  } else if (sample_list.empty()) {
    g_logger.error("Sample list is required to be specified");
    return ERR;
  } else if (TileDBUtils::is_file(workspace) || TileDBUtils::is_dir(workspace)) {
    if (!overwrite_workspace) {
      g_logger.error("Workspace {} exists, retry with overwrite-workspace option", workspace);
      return ERR;
    }
  }

  try {
    if (TileDBUtils::create_workspace(workspace, overwrite_workspace) != TILEDB_OK) {
      g_logger.error("Could not create workspace {}", workspace);
      return ERR;
    }

    workspace = TileDBUtils::real_dir(workspace);

    std::string vidmap_output = SLASHIFY(workspace)+"vidmap.json";
    std::string callset_output = SLASHIFY(workspace)+"callset.json";
    std::string loader_json = SLASHIFY(workspace)+"loader.json";
    std::string merged_header = SLASHIFY(workspace)+"vcfheader.vcf";

    if (merge_headers_and_generate_callset(sample_list, merged_header, callset_output)) {
      return -1;
    }
    if (generate_json(merged_header, UNSLASHIFY(workspace), vidmap_output, loader_json, callset_output)) {
      return -1;
    }
    
    google::protobuf::ShutdownProtobufLibrary();

  } catch(const GenomicsDBConfigException& genomicsdb_ex) {
    g_logger.error(genomicsdb_ex.what());
    return ERR;
  } catch (const std::exception& ex) {
    g_logger.error(ex.what());
    return ERR;
  } 
  
  g_logger.info("Success!");
  return OK;
}