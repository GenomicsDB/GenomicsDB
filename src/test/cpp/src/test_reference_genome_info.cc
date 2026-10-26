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
#include "vcf_adapter.h"

//Pseudo-random bases, long enough that reading them all takes several reads from the FASTA
static std::string make_bases(size_t length, unsigned seed) {
  std::string bases(length, 'N');
  for (auto& base : bases) {
    seed = seed * 1103515245u + 12345u;
    base = "ACGT"[(seed >> 16) & 3u];
  }
  return bases;
}

//Writes a FASTA with 60 bases per line; the reference builds its .fai index when it loads it
static std::string write_fasta(TempDir& temp_dir, const std::vector<std::pair<std::string, std::string>>& contigs) {
  std::string fasta;
  for (const auto& [name, bases] : contigs) {
    fasta += ">" + name + "\n";
    for (size_t i = 0; i < bases.length(); i += 60)
      fasta += bases.substr(i, 60) + "\n";
  }
  auto filename = temp_dir.append("reference.fa");
  REQUIRE(TileDBUtils::write_file(filename, fasta.c_str(), fasta.length()) == TILEDB_OK);
  return filename;
}

TEST_CASE_METHOD(TempDir, "reference returns the base at each position of a contig", "[reference_genome_info]") {
  auto bases = make_bases(100000, 1);
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", bases}}));

  std::string read_back;
  for (int64_t pos = 0; pos < static_cast<int64_t>(bases.length()); ++pos)
    read_back += reference.get_reference_base_at_position("chr1", pos);
  CHECK(read_back == bases);
}

TEST_CASE_METHOD(TempDir, "reference returns the base at a position before the last one read", "[reference_genome_info]") {
  auto bases = make_bases(100000, 2);
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", bases}}));

  CHECK(reference.get_reference_base_at_position("chr1", 70000) == bases[70000]);
  CHECK(reference.get_reference_base_at_position("chr1", 10) == bases[10]);
  CHECK(reference.get_reference_base_at_position("chr1", 69999) == bases[69999]);
}

TEST_CASE_METHOD(TempDir, "reference returns the base of the contig asked for after switching contigs", "[reference_genome_info]") {
  auto chr1 = make_bases(1000, 3);
  auto chr2 = make_bases(1000, 4);
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", chr1}, {"chr2", chr2}}));

  for (int64_t pos : {0, 500, 999}) {
    CHECK(reference.get_reference_base_at_position("chr1", pos) == chr1[pos]);
    CHECK(reference.get_reference_base_at_position("chr2", pos) == chr2[pos]);
  }
}

TEST_CASE_METHOD(TempDir, "reference keeps the case of the FASTA", "[reference_genome_info]") {
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", "ACGTacgtNn"}}));

  CHECK(reference.get_reference_base_at_position("chr1", 3) == 'T');
  CHECK(reference.get_reference_base_at_position("chr1", 4) == 'a');
  CHECK(reference.get_reference_base_at_position("chr1", 9) == 'n');
}

TEST_CASE_METHOD(TempDir, "reference returns N for a contig it does not have", "[reference_genome_info]") {
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", "ACGTACGTAC"}}));

  CHECK(reference.get_reference_base_at_position("chr1", 1) == 'C');
  CHECK(reference.get_reference_base_at_position("chr2", 1) == 'N');
}

TEST_CASE_METHOD(TempDir, "reference returns N for a position past the end of a contig", "[reference_genome_info]") {
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", "ACGTACGTAG"}, {"chr2", "TTTTTTTTTT"}}));

  CHECK(reference.get_reference_base_at_position("chr1", 9) == 'G');
  CHECK(reference.get_reference_base_at_position("chr1", 10) == 'N');
  CHECK(reference.get_reference_base_at_position("chr1", 5000) == 'N');
  CHECK(reference.get_reference_base_at_position("chr2", 10) == 'N');
}

TEST_CASE_METHOD(TempDir, "reference returns N for a negative position", "[reference_genome_info]") {
  ReferenceGenomeInfo reference;
  reference.initialize(write_fasta(*this, {{"chr1", "ACGTACGTAC"}}));

  CHECK(reference.get_reference_base_at_position("chr1", -1) == 'N');
}
