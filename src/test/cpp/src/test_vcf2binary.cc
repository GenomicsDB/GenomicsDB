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
