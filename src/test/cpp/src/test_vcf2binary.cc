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

#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include "test_base.h"
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

//A bgzipped, indexed single-sample VCF on contig 1 with a record at each of the given positions, long enough to take
//many reads with htslib's default buffer
static std::string write_indexed_vcf(TempDir& temp_dir, const std::vector<int>& positions) {
  std::string vcf = "##fileformat=VCFv4.2\n"
                    "##contig=<ID=1,length=249250621>\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n";
  for (auto pos : positions)
    vcf += "1\t"+std::to_string(pos)+"\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
  auto filename = temp_dir.append("indexed.vcf.gz");
  auto bgzf = bgzf_open(filename.c_str(), "w");
  REQUIRE(bgzf != nullptr);
  REQUIRE(bgzf_write(bgzf, vcf.c_str(), vcf.length()) == static_cast<ssize_t>(vcf.length()));
  REQUIRE(bgzf_close(bgzf) == 0);
  REQUIRE(tbx_index_build(filename.c_str(), 0, &tbx_conf_vcf) == 0);
  return filename;
}

//Reads the file from its start with a VCFReader that has the given read buffer size, and returns the 1-based
//positions of its records
static std::vector<int> read_positions(const std::string& vcf_filename, const int read_buffer_size) {
  FileBasedVidMapper vid_mapper(ctests_input_dir+"vid.json");
  std::vector<std::vector<std::string>> vcf_fields;
  vid_mapper.build_vcf_fields_vectors(vcf_fields);
  VCFReader reader(read_buffer_size);
  reader.initialize(vcf_filename.c_str(), vcf_fields, &vid_mapper, false);
  std::vector<int> positions;
  reader.seek_read_advance("1", 0, true);
  for (auto line = reader.get_line(); line; line = reader.get_line()) {
    positions.push_back(line->pos+1);
    reader.read_and_advance();
  }
  return positions;
}

static std::vector<int> every_tenth_position(const int num_records) {
  std::vector<int> positions;
  for (auto i=1; i<=num_records; ++i)
    positions.push_back(10*i);
  return positions;
}

TEST_CASE_METHOD(TempDir, "VCFReader reads every record with htslib's read buffer", "[vcf_reader_read_buffer_size]") {
  auto positions = every_tenth_position(20000);
  CHECK(read_positions(write_indexed_vcf(*this, positions), 0) == positions);
}

TEST_CASE_METHOD(TempDir, "VCFReader reads every record with a larger read buffer", "[vcf_reader_read_buffer_size]") {
  auto positions = every_tenth_position(20000);
  CHECK(read_positions(write_indexed_vcf(*this, positions), 1048576) == positions);
}

TEST_CASE_METHOD(TempDir, "VCFReader reads every record with a 1-byte read buffer", "[vcf_reader_read_buffer_size]") {
  auto positions = every_tenth_position(20000);
  CHECK(read_positions(write_indexed_vcf(*this, positions), 1) == positions);
}
