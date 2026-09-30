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
 * Tests the allele index look-up tables
 */

#include <catch2/catch.hpp>

#include "lut.h"

static const int64_t num_inputs = 4;

// True if no input has a mapping to or from any of the first num_alleles alleles
static bool all_missing(const CombineAllelesLUT& lut, const int64_t num_alleles) {
  for (auto input=0; input<num_inputs; ++input)
    for (auto allele=0; allele<num_alleles; ++allele)
      if (!CombineAllelesLUT::is_missing_value(lut.get_merged_idx_for_input(input, allele))
          || !CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(input, allele)))
        return false;
  return true;
}

TEST_CASE("lut reset_luts clears every mapping added since the last reset", "[lut_reset_clears_mappings]") {
  CombineAllelesLUT lut(num_inputs);
  lut.add_input_merged_idx_pair(0, 0, 0);
  lut.add_input_merged_idx_pair(1, 1, 2);
  lut.add_input_merged_idx_pair(3, 2, 1);
  CHECK(lut.get_merged_idx_for_input(1, 1) == 2);
  CHECK(lut.get_input_idx_for_merged(3, 1) == 2);
  lut.reset_luts();
  CHECK(all_missing(lut, 10));
}

TEST_CASE("lut reset_luts clears mappings in columns added by a resize", "[lut_reset_after_resize]") {
  CombineAllelesLUT lut(num_inputs);
  lut.add_input_merged_idx_pair(2, 1, 1);
  lut.resize_luts_if_needed(30);
  lut.add_input_merged_idx_pair(0, 25, 3);
  lut.add_input_merged_idx_pair(3, 4, 29);
  lut.reset_luts();
  CHECK(all_missing(lut, 30));
}

TEST_CASE("lut mappings added after a reset are kept while older ones stay cleared", "[lut_mappings_after_reset]") {
  CombineAllelesLUT lut(num_inputs);
  lut.resize_luts_if_needed(20);
  lut.add_input_merged_idx_pair(0, 15, 18);
  lut.reset_luts();
  lut.add_input_merged_idx_pair(1, 0, 0);
  lut.add_input_merged_idx_pair(1, 1, 2);
  CHECK(lut.get_merged_idx_for_input(1, 1) == 2);
  CHECK(lut.get_input_idx_for_merged(1, 2) == 1);
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_merged_idx_for_input(0, 15)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(0, 18)));
  lut.reset_luts();
  CHECK(all_missing(lut, 20));
}
