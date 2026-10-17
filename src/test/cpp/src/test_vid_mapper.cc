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

#include "vid_mapper.h"

TEST_CASE("callset indexes stored as JSON strings are read as integers", "[vid_mapper_callset_int64_strings]") {
  //Protobuf JSON, as in the callset mapping GATK writes, stores int64 fields as strings
  FileBasedVidMapper vid_mapper;
  vid_mapper.parse_callsets_json(R"({"callsets": [
      {"sample_name": "S0", "row_idx": "0", "idx_in_file": "0", "filename": "a.vcf.gz"},
      {"sample_name": "S1", "row_idx": "1", "idx_in_file": "1", "filename": "a.vcf.gz"},
      {"sample_name": "S2", "row_idx": "2", "idx_in_file": "0", "filename": "b.vcf.gz"}
    ]})", false);

  int64_t row_idx = -1;
  REQUIRE(vid_mapper.get_tiledb_row_idx(row_idx, "S1"));
  CHECK(row_idx == 1);
  CHECK(vid_mapper.get_callset_info(1).m_idx_in_file == 1);
  REQUIRE(vid_mapper.get_tiledb_row_idx(row_idx, "S2"));
  CHECK(row_idx == 2);
  CHECK(vid_mapper.get_callset_info(2).m_idx_in_file == 0);
}

TEST_CASE("callset index string that is not an integer is rejected", "[vid_mapper_callset_non_integer_string]") {
  FileBasedVidMapper vid_mapper;
  CHECK_THROWS_WITH(vid_mapper.parse_callsets_json(R"({"callsets": [
      {"sample_name": "S0", "row_idx": "0", "idx_in_file": "1x", "filename": "a.vcf.gz"}
    ]})", false), Catch::Contains("\"idx_in_file\" for sample/callset S0 must be an integer"));
}

TEST_CASE("callset index that is not an integer or a string is rejected", "[vid_mapper_callset_non_integer_value]") {
  FileBasedVidMapper vid_mapper;
  CHECK_THROWS_WITH(vid_mapper.parse_callsets_json(R"({"callsets": [
      {"sample_name": "S0", "row_idx": 1.5, "filename": "a.vcf.gz"}
    ]})", false), Catch::Contains("\"row_idx\" for sample/callset S0 must be an integer"));
}

TEST_CASE("contig lengths and offsets stored as JSON strings are read as integers", "[vid_mapper_contig_int64_strings]") {
  //Protobuf JSON, as in the vid mapping GATK writes, stores int64 fields as strings
  rapidjson::Document json_doc;
  json_doc.Parse(R"({
      "fields": [{"name": "PASS", "type": "int"}],
      "contigs": [
        {"name": "chr1", "length": "248956422", "tiledb_column_offset": "0"},
        {"name": "chr14", "length": "107043718", "tiledb_column_offset": "2191407310"}
      ]})");
  FileBasedVidMapper vid_mapper(json_doc);

  ContigInfo contig_info;
  REQUIRE(vid_mapper.get_contig_info("chr14", contig_info));
  CHECK(contig_info.m_length == 107043718);
  CHECK(contig_info.m_tiledb_column_offset == 2191407310);
  int64_t position = -1;
  REQUIRE(vid_mapper.get_tiledb_position(position, "chr14", 100));
  CHECK(position == 2191407410);
  REQUIRE(vid_mapper.get_contig_info("chr1", contig_info));
  CHECK(contig_info.m_length == 248956422);
  CHECK(contig_info.m_tiledb_column_offset == 0);
}

TEST_CASE("contig length string that is not an integer is rejected", "[vid_mapper_contig_non_integer_string]") {
  rapidjson::Document json_doc;
  json_doc.Parse(R"({
      "fields": [{"name": "PASS", "type": "int"}],
      "contigs": [{"name": "chr1", "length": "248956422bp", "tiledb_column_offset": "0"}]})");
  CHECK_THROWS_WITH(FileBasedVidMapper(json_doc), Catch::Contains("\"length\" for contig chr1 must be an integer"));
}

TEST_CASE("contig offset that is not an integer or a string is rejected", "[vid_mapper_contig_non_integer_value]") {
  rapidjson::Document json_doc;
  json_doc.Parse(R"({
      "fields": [{"name": "PASS", "type": "int"}],
      "contigs": [{"name": "chr1", "length": 248956422, "tiledb_column_offset": [0]}]})");
  CHECK_THROWS_WITH(FileBasedVidMapper(json_doc),
                    Catch::Contains("\"tiledb_column_offset\" for contig chr1 must be an integer"));
}
