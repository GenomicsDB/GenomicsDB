/**
 * The MIT License (MIT)
 * Copyright (c) 2026 Tim Fennell
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
 * Tests the combined gVCF produced by BroadCombinedGVCFOperator
 */

#include <catch2/catch.hpp>

#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <zlib.h>

#include <algorithm>
#include <map>

#include "genomicsdb.h"
#include "test_base.h"
#include "tiledb_loader.h"
#include "tiledb_utils.h"

static std::string tests_src_dir(GENOMICSDB_TESTS_SRC_DIR);

static std::string ref_block(const int begin, const int end, const std::string& GT = "0/0") {
  return "1\t" + std::to_string(begin) + "\t.\tA\t<NON_REF>\t.\t.\tEND=" + std::to_string(end) + "\tGT\t" + GT;
}

//A call at pos with REF A, the comma-separated alts followed by <NON_REF>, and the given FORMAT keys and values
static std::string variant(const int pos, const std::string& alts, const std::string& keys, const std::string& values) {
  return "1\t" + std::to_string(pos) + "\t.\tA\t" + alts + ",<NON_REF>\t.\t.\t.\t" + keys + "\t" + values;
}

//as_raw_mq and as_raw_mq_rank_sum hold one entry each for REF, alt and <NON_REF>
static std::string snp(const int pos, const std::string& alt, const std::string& as_raw_mq,
                       const std::string& as_raw_mq_rank_sum) {
  return "1\t" + std::to_string(pos) + "\t.\tA\t" + alt + ",<NON_REF>\t.\t.\tAS_RAW_MQ=" + as_raw_mq
         + ";AS_RAW_MQRankSum=" + as_raw_mq_rank_sum + "\tGT\t0/1";
}

//Writes a bgzipped and tabix-indexed single-sample gVCF, returning its path
static std::string write_gvcf(TempDir& temp_dir, const std::string& sample, const std::vector<std::string>& records) {
  std::string vcf = "##fileformat=VCFv4.2\n"
                    "##ALT=<ID=NON_REF,Description=\"Any other allele\">\n"
                    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n"
                    "##INFO=<ID=AS_RAW_MQ,Number=1,Type=String,Description=\"Allele-specific raw MQ\">\n"
                    "##INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description=\"Allele-specific MQ rank sums\">\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods\">\n"
                    "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype\">\n"
                    "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID\">\n"
                    "##contig=<ID=1,length=249250621>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n";
  for (const auto& record : records)
    vcf += record + "\n";
  auto filename = temp_dir.append(sample + ".g.vcf.gz");
  auto* fp = bgzf_open(filename.c_str(), "w");
  REQUIRE(fp);
  REQUIRE(bgzf_write(fp, vcf.c_str(), vcf.length()) == static_cast<ssize_t>(vcf.length()));
  REQUIRE(bgzf_close(fp) == 0);
  REQUIRE(tbx_index_build(filename.c_str(), 0, &tbx_conf_vcf) == 0);
  return filename;
}

//Loads the samples' gVCFs, in row order, into a new workspace and returns the path of the loader JSON
static std::string load(TempDir& temp_dir, const std::vector<std::pair<std::string, std::vector<std::string>>>& samples) {
  std::string callsets = "{\"callsets\": {";
  for (auto row=0u; row<samples.size(); ++row) {
    auto filename = write_gvcf(temp_dir, samples[row].first, samples[row].second);
    callsets += (row ? ", \"" : "\"") + samples[row].first + "\": {\"row_idx\": " + std::to_string(row)
                + ", \"idx_in_file\": 0, \"filename\": \"" + filename + "\"}";
  }
  callsets += "}}";
  auto callsets_json = temp_dir.append("callsets.json");
  REQUIRE(TileDBUtils::write_file(callsets_json, callsets.c_str(), callsets.length()) == TILEDB_OK);
  auto workspace = temp_dir.append("ws");
  REQUIRE(TileDBUtils::create_workspace(workspace) == TILEDB_OK);
  const std::string loader = R"({"row_based_partitioning": false,
      "column_partitions": [{"begin": 0, "workspace": ")" + workspace + R"(", "array_name": "as_info"}],
      "callset_mapping_file": ")" + callsets_json + R"(",
      "vid_mapping_file": ")" + tests_src_dir + R"(inputs/vid_all_asa.json",
      "reference_genome": ")" + tests_src_dir + R"(inputs/chr1_10MB.fasta.gz",
      "size_per_column_partition": 16384, "treat_deletions_as_intervals": true, "num_parallel_vcf_files": 1,
      "discard_vcf_index": true, "produce_combined_vcf": false, "produce_tiledb_array": true,
      "delete_and_create_tiledb_array": true, "segment_size": 1048576, "num_cells_per_tile": 3})";
  auto loader_json = temp_dir.append("loader.json");
  REQUIRE(TileDBUtils::write_file(loader_json, loader.c_str(), loader.length()) == TILEDB_OK);
  VCF2TileDBLoader loader_obj(loader_json, 0);
  loader_obj.read_all();
  return loader_json;
}

static std::vector<std::string> split(const std::string& str, const char delimiter) {
  std::vector<std::string> tokens;
  for (size_t begin = 0, end = 0; end != std::string::npos; begin = end + 1) {
    end = str.find(delimiter, begin);
    tokens.push_back(str.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
  }
  return tokens;
}

//Queries the combined gVCF, with GT, of all rows for the attributes (a JSON list) and returns each record's columns
//keyed by position
static std::map<int, std::vector<std::string>> query_records(TempDir& temp_dir, const std::string& loader_json,
                                                             const std::string& attributes) {
  const std::string query = R"({"workspace": ")" + temp_dir.append("ws") + R"(", "array_name": "as_info",
      "query_column_ranges": [{"range_list": [{"low": 0, "high": 1000000000}]}],
      "query_row_ranges": [{"range_list": [{"low": 0, "high": 2}]}],
      "reference_genome": ")" + tests_src_dir + R"(inputs/chr1_10MB.fasta.gz",
      "vcf_header_filename": [")" + tests_src_dir + R"(inputs/template_vcf_header.vcf"],
      "produce_GT_field": true, "attributes": )" + attributes + R"(, "segment_size": 1048576})";
  auto vcf_file = temp_dir.append("combined.vcf.gz");
  GenomicsDB gdb(query, GenomicsDB::JSON_STRING, loader_json);
  gdb.generate_vcf(vcf_file, "z", true);
  std::map<int, std::vector<std::string>> records;
  gzFile input = gzopen(vcf_file.c_str(), "r");
  REQUIRE(input);
  char buffer[65536];
  while (gzgets(input, buffer, sizeof(buffer))) {
    std::string line(buffer);
    if (line.empty() || line[0] == '#')
      continue;
    if (line.back() == '\n')
      line.pop_back();
    auto columns = split(line, '\t');
    REQUIRE(columns.size() >= 8u);
    records[std::stoi(columns[1])] = columns;
  }
  gzclose(input);
  return records;
}

//Returns the named field's value at each position whose INFO column has it
static std::map<int, std::string> INFO_field_values(const std::map<int, std::vector<std::string>>& records,
                                                    const std::string& field) {
  std::map<int, std::string> values;
  const auto key = field + "=";
  for (const auto& position_columns_pair : records) {
    for (const auto& key_value : split(position_columns_pair.second[7], ';'))
      if (key_value.compare(0, key.length(), key) == 0)
        values[position_columns_pair.first] = key_value.substr(key.length());
  }
  return values;
}

//Returns each sample's value of the named FORMAT field at each position whose record has the field
static std::map<int, std::vector<std::string>> FORMAT_field_values(
  const std::map<int, std::vector<std::string>>& records, const std::string& field) {
  std::map<int, std::vector<std::string>> values;
  for (const auto& position_columns_pair : records) {
    const auto& columns = position_columns_pair.second;
    REQUIRE(columns.size() > 9u);
    const auto keys = split(columns[8], ':');
    const auto key_idx = static_cast<size_t>(std::find(keys.begin(), keys.end(), field) - keys.begin());
    if (key_idx == keys.size())
      continue;
    auto& sample_values = values[position_columns_pair.first];
    for (auto i=9u; i<columns.size(); ++i) {
      const auto sample_fields = split(columns[i], ':');
      sample_values.push_back(key_idx < sample_fields.size() ? sample_fields[key_idx] : "");
    }
  }
  return values;
}

TEST_CASE_METHOD(TempDir, "combined gVCF allele-specific INFO fields combine only the samples that have values",
                 "[broad_combined_gvcf_allele_specific_INFO]") {
  //Position 100: only S0, whose values must not carry over to 200, where S1 and S2 but not S0, the first row, have
  //values. 300: only S2, the last row. 400: S0 and S1 with different alts, so each sample's missing alt takes
  //its <NON_REF> value
  auto loader_json = load(*this, {
    {"S0", {ref_block(1, 99), snp(100, "C", "1.00|2.00|3.00", "|0.5,1|NaN"), ref_block(101, 399),
            snp(400, "C", "1.00|2.00|3.00", "|0.5,1|NaN"), ref_block(401, 1000)}},
    {"S1", {ref_block(1, 199), snp(200, "C", "10.00|20.00|30.00", "|0.5,2,1.5,1|NaN"), ref_block(201, 399),
            snp(400, "G", "10.00|20.00|30.00", "|1.5,2|NaN"), ref_block(401, 1000)}},
    {"S2", {ref_block(1, 199), snp(200, "C", "100.00|200.00|300.00", "|1.5,4|NaN"), ref_block(201, 299),
            snp(300, "T", "100.00|200.00|300.00", "|2.5,1|NaN"), ref_block(301, 1000)}}
  });

  auto records = query_records(*this, loader_json, R"(["END", "REF", "ALT", "GT", "AS_RAW_MQ", "AS_RAW_MQRankSum"])");

  const std::map<int, std::string> expected_AS_RAW_MQ = {
    {100, "1.000|2.000|3.000"},
    {200, "110.000|220.000|330.000"},
    {300, "100.000|200.000|300.000"},
    {400, "11.000|32.000|23.000|33.000"}
  };
  CHECK(INFO_field_values(records, "AS_RAW_MQ") == expected_AS_RAW_MQ);

  const std::map<int, std::string> expected_AS_RAW_MQRankSum = {
    {100, "|0.500,1|"},
    {200, "|0.500,2,1.500,5|"},
    {300, "|2.500,1|"},
    {400, "|0.500,1|1.500,2|"}
  };
  CHECK(INFO_field_values(records, "AS_RAW_MQRankSum") == expected_AS_RAW_MQRankSum);
}

TEST_CASE_METHOD(TempDir, "combined gVCF AD and PL are missing for samples without them, even after a record where they had them",
                 "[broad_combined_gvcf_FORMAT_missing]") {
  //S0 has AD and PL at 100 only, S1 at 200 only and S2 at both; the others are in reference blocks
  auto loader_json = load(*this, {
    {"S0", {ref_block(1, 99), variant(100, "C", "GT:AD:PL", "0/1:5,3,1:30,0,50,40,60,90"), ref_block(101, 1000)}},
    {"S1", {ref_block(1, 199), variant(200, "C", "GT:AD:PL", "0/1:4,6,2:20,0,70,35,80,95"), ref_block(201, 1000)}},
    {"S2", {ref_block(1, 99), variant(100, "C", "GT:AD:PL", "1/1:0,7,1:90,20,0,85,25,99"), ref_block(101, 199),
            variant(200, "C", "GT:AD:PL", "0/1:3,3,0:25,0,25,40,40,80"), ref_block(201, 1000)}}
  });

  auto records = query_records(*this, loader_json, R"(["END", "REF", "ALT", "GT", "AD", "PL"])");

  auto AD = FORMAT_field_values(records, "AD");
  CHECK(AD.at(100) == std::vector<std::string>({"5,3,1", ".", "0,7,1"}));
  CHECK(AD.at(200) == std::vector<std::string>({".", "4,6,2", "3,3,0"}));
  auto PL = FORMAT_field_values(records, "PL");
  CHECK(PL.at(100) == std::vector<std::string>({"30,0,50,40,60,90", ".", "90,20,0,85,25,99"}));
  CHECK(PL.at(200) == std::vector<std::string>({".", "20,0,70,35,80,95", "25,0,25,40,40,80"}));
}

TEST_CASE_METHOD(TempDir, "combined gVCF string FORMAT fields keep each sample's value and are missing for samples without them",
                 "[broad_combined_gvcf_FORMAT_strings]") {
  auto loader_json = load(*this, {
    {"S0", {ref_block(1, 99), variant(100, "C", "GT:PGT:PID", "0|1:0|1:100_A_C"), ref_block(101, 1000)}},
    {"S1", {ref_block(1, 1000)}},
    {"S2", {ref_block(1, 99), variant(100, "C", "GT:PGT:PID", "1|0:1|0:90_T_GAA"), ref_block(101, 1000)}}
  });

  auto records = query_records(*this, loader_json, R"(["END", "REF", "ALT", "GT", "PGT", "PID"])");

  CHECK(FORMAT_field_values(records, "PGT").at(100) == std::vector<std::string>({"0|1", ".", "1|0"}));
  CHECK(FORMAT_field_values(records, "PID").at(100) == std::vector<std::string>({"100_A_C", ".", "90_T_GAA"}));
}

//At 100 S0 has alt G and S1 alts C and T, so the merged alleles are A,G,C,T,<NON_REF>. At 200 S1 has alt C and S2
//alt G, so they are A,C,G,<NON_REF>. The value of <NON_REF> in each AD differs from its other values
static std::string load_samples_with_different_alts(TempDir& temp_dir) {
  return load(temp_dir, {
    {"S0", {ref_block(1, 99), variant(100, "G", "GT:AD", "0/1:5,3,1"), ref_block(101, 1000)}},
    {"S1", {ref_block(1, 99), variant(100, "C,T", "GT:AD", "1/2:1,4,6,2"), ref_block(101, 199),
            variant(200, "C", "GT:AD", "0/1:7,2,1"), ref_block(201, 1000)}},
    {"S2", {ref_block(1, 199, "./."), variant(200, "G", "GT:AD", "1/1:3,8,1"), ref_block(201, 1000)}}
  });
}

TEST_CASE_METHOD(TempDir, "combined gVCF AD takes each sample's value for <NON_REF> for merged alleles the sample lacks",
                 "[broad_combined_gvcf_AD_remap]") {
  auto loader_json = load_samples_with_different_alts(*this);

  auto records = query_records(*this, loader_json, R"(["END", "REF", "ALT", "GT", "AD"])");

  auto AD = FORMAT_field_values(records, "AD");
  CHECK(AD.at(100) == std::vector<std::string>({"5,3,1,1,1", "1,2,4,6,2", "."}));
  //At 100 S1's C was merged allele 2, which at 200 is G, which S1 lacks
  CHECK(AD.at(200) == std::vector<std::string>({".", "7,2,1,1", "3,1,8,1"}));
}

TEST_CASE_METHOD(TempDir, "combined gVCF GT alleles are renumbered to the merged alleles, keeping REF and no-calls",
                 "[broad_combined_gvcf_GT_remap]") {
  auto loader_json = load_samples_with_different_alts(*this);

  auto records = query_records(*this, loader_json, R"(["END", "REF", "ALT", "GT"])");

  auto GT = FORMAT_field_values(records, "GT");
  CHECK(GT.at(100) == std::vector<std::string>({"0/1", "2/3", "./."}));
  CHECK(GT.at(200) == std::vector<std::string>({"0/0", "0/1", "2/2"}));
}
