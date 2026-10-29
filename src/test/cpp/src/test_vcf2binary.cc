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
 */

#include <catch2/catch.hpp>

#include "genomicsdb.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include "test_base.h"
#include "tiledb_loader.h"
#include "tiledb_utils.h"
#include "vcf2binary.h"
#include "vid_mapper.h"

static std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

//A single-sample VCF on contig 1 of the ctests vid mapping
static std::string write_single_sample_vcf(TempDir& temp_dir) {
  const std::string vcf = "##fileformat=VCFv4.2\n"
                          "##contig=<ID=1,length=249250621>\n"
                          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n"
                          "1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
  auto filename = temp_dir.append("single_sample.vcf");
  REQUIRE(TileDBUtils::write_file(filename, vcf.c_str(), vcf.length()) == TILEDB_OK);
  return filename;
}

//A single-sample VCF whose header IDX values skip 2, which leaves slot 2 of htslib's header dictionary empty
static std::string write_vcf_with_header_idx_gap(TempDir& temp_dir) {
  const std::string vcf = "##fileformat=VCFv4.2\n"
                          "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n"
                          "##contig=<ID=1,length=249250621>\n"
                          "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=1>\n"
                          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\",IDX=3>\n"
                          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n"
                          "1\t100\t.\tA\tC\t.\t.\tDP=5\tGT\t0/1\n";
  auto filename = temp_dir.append("header_idx_gap.vcf");
  REQUIRE(TileDBUtils::write_file(filename, vcf.c_str(), vcf.length()) == TILEDB_OK);
  return filename;
}

//Creates the VCF2Binary for the file, with its sample at idx_in_file in the callset mapping
static void create_vcf2binary(const std::string& vcf_filename, const int idx_in_file) {
  FileBasedVidMapper vid_mapper(ctests_input_dir+"vid.json");
  vid_mapper.parse_callsets_json("{\"callsets\": [{\"sample_name\": \"S0\", \"row_idx\": 0, \"idx_in_file\": "
                                 +std::to_string(idx_in_file)+", \"filename\": \""+vcf_filename+"\"}]}", false);
  std::vector<std::vector<std::string>> vcf_fields;
  vid_mapper.build_vcf_fields_vectors(vcf_fields);
  int64_t file_idx = -1;
  REQUIRE(vid_mapper.get_global_file_idx(vcf_filename, file_idx));
  //close_file so that the file is only opened to read its header, which needs no index
  VCF2Binary vcf2binary(vcf_filename, vcf_fields, file_idx, vid_mapper, {ColumnRange(0, INT64_MAX-1)}, 1024u,
                        true, false, true, true);
}

TEST_CASE_METHOD(TempDir, "VCF2Binary accepts a callset mapping to a sample in the file", "[vcf2binary_sample_index_in_file]") {
  auto vcf_filename = write_single_sample_vcf(*this);
  CHECK_NOTHROW(create_vcf2binary(vcf_filename, 0));
}

TEST_CASE_METHOD(TempDir, "VCF2Binary rejects a callset mapping to a sample past the end of the file", "[vcf2binary_sample_index_past_file]") {
  auto vcf_filename = write_single_sample_vcf(*this);
  CHECK_THROWS_WITH(create_vcf2binary(vcf_filename, 1),
                    Catch::Contains("Callset mapping refers to sample index 1 of "+vcf_filename+", which has 1 sample(s)"));
}

TEST_CASE_METHOD(TempDir, "VCF2Binary accepts a header with an empty dictionary slot", "[vcf2binary_header_idx_gap]") {
  auto vcf_filename = write_vcf_with_header_idx_gap(*this);
  CHECK_NOTHROW(create_vcf2binary(vcf_filename, 0));
}

//Writes a single-sample gVCF on contig 1 with the given records, bgzipped and indexed as the importer needs it
static std::string write_indexed_gvcf(TempDir& temp_dir, const std::string& records) {
  const std::string vcf = "##fileformat=VCFv4.2\n"
                          "##ALT=<ID=NON_REF,Description=\"Any other allele\">\n"
                          "##contig=<ID=1,length=249250621>\n"
                          "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the record\">\n"
                          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n"
                          + records;
  auto filename = temp_dir.append("sample.g.vcf.gz");
  auto bgzf = bgzf_open(filename.c_str(), "w");
  REQUIRE(bgzf != nullptr);
  REQUIRE(bgzf_write(bgzf, vcf.c_str(), vcf.length()) == static_cast<ssize_t>(vcf.length()));
  REQUIRE(bgzf_close(bgzf) == 0);
  REQUIRE(tbx_index_build(filename.c_str(), 0, &tbx_conf_vcf) == 0);
  return filename;
}

//Collects each call as "begin-end ALT", with 1-based positions and ALT as stored, & standing for <NON_REF>
class CallALTCollector : public GenomicsDBVariantCallProcessor {
 public:
  void process(const interval_t& interval) override {}
  void process(const std::string& sample_name, const int64_t* coordinates, const genomic_interval_t& genomic_interval,
               const std::vector<genomic_field_t>& genomic_fields) override {
    for (const auto& field : genomic_fields)
      if (field.name == "ALT")
        m_calls.push_back(std::to_string(genomic_interval.interval.first) + "-"
                          + std::to_string(genomic_interval.interval.second) + " "
                          + field.str_value());
  }
  std::vector<std::string> m_calls;
};

//Imports a gVCF with the given records into a new workspace, treating deletions as intervals, and returns its calls
static std::vector<std::string> import_and_query_calls(TempDir& temp_dir, const std::string& records) {
  auto vcf_filename = write_indexed_gvcf(temp_dir, records);
  auto callsets = temp_dir.append("callsets.json");
  const std::string callsets_json = "{\"callsets\": [{\"sample_name\": \"S0\", \"row_idx\": 0, \"idx_in_file\": 0, "
                                    "\"filename\": \""+vcf_filename+"\"}]}";
  REQUIRE(TileDBUtils::write_file(callsets, callsets_json.c_str(), callsets_json.length()) == TILEDB_OK);
  auto workspace = temp_dir.append("ws");
  REQUIRE(TileDBUtils::create_workspace(workspace, false) == TILEDB_OK);
  auto loader = temp_dir.append("loader.json");
  const std::string loader_json = "{\"column_partitions\": [{\"begin\": 0, \"workspace\": \""+workspace+"\", "
                                  "\"array_name\": \"calls\"}], "
                                  "\"callset_mapping_file\": \""+callsets+"\", "
                                  "\"vid_mapping_file\": \""+ctests_input_dir+"vid.json\", "
                                  "\"size_per_column_partition\": 16384, "
                                  "\"treat_deletions_as_intervals\": true, "
                                  "\"delete_and_create_tiledb_array\": true}";
  REQUIRE(TileDBUtils::write_file(loader, loader_json.c_str(), loader_json.length()) == TILEDB_OK);

  VCF2TileDBLoader(loader, 0).read_all();

  CallALTCollector collector;
  GenomicsDB(workspace, callsets, ctests_input_dir+"vid.json").query_variant_calls(collector, "calls");
  return collector.m_calls;
}

TEST_CASE_METHOD(TempDir, "VCF import stores a <NON_REF> ALT allele as &", "[vcf2binary_alleles]") {
  CHECK(import_and_query_calls(*this, "1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=110\tGT\t0/0\n")
        == std::vector<std::string>{"100-110 &"});
}

TEST_CASE_METHOD(TempDir, "VCF import stores ALT alleles other than <NON_REF> as written", "[vcf2binary_alleles]") {
  CHECK(import_and_query_calls(*this, "1\t100\t.\tA\tG,<*>,<DEL>,<NON_REF>\t.\t.\t.\tGT\t0/1\n")
        == std::vector<std::string>{"100-100 G|<*>|<DEL>|&"});
}

TEST_CASE_METHOD(TempDir, "VCF import ends a deletion without END at the end of its REF", "[vcf2binary_alleles]") {
  CHECK(import_and_query_calls(*this, "1\t100\t.\tACGT\tA,<NON_REF>\t.\t.\t.\tGT\t0/1\n")
        == std::vector<std::string>{"100-103 A|&"});
}

TEST_CASE_METHOD(TempDir, "VCF import ends a spanning deletion without END at the end of its REF", "[vcf2binary_alleles]") {
  CHECK(import_and_query_calls(*this, "1\t100\t.\tAC\t*,<NON_REF>\t.\t.\t.\tGT\t0/1\n")
        == std::vector<std::string>{"100-101 *|&"});
}

TEST_CASE_METHOD(TempDir, "VCF import ends a record with an insertion and a deletion at the end of its REF", "[vcf2binary_alleles]") {
  CHECK(import_and_query_calls(*this, "1\t100\t.\tACGT\tACGTT,A,<NON_REF>\t.\t.\t.\tGT\t1/2\n")
        == std::vector<std::string>{"100-103 ACGTT|A|&"});
}
