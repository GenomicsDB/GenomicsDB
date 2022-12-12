/**
 *
 * The MIT License (MIT)
 * Copyright (c) 2021 Omics Data Automation, Inc.
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
 * Test the optimized GT remapper 
 *
 */

#include <catch2/catch.hpp>
#include "alleles_combiner_template_definition.h"
#include "gt_remapper_template_definition.h"
#include "test_valid_row_tracker.h"
#include "test_remapped_data_receiver.h"
#include "test_data_provider_for_remapper.h"

class ValidRowTrackerAndGTDataProviderForUnitTest : public ValidRowTrackerForUnitTest, public RemappedDataReceiverForUnitTest,
  public TestDataProviderForRemapper<int>
{
  public:
    ValidRowTrackerAndGTDataProviderForUnitTest(const size_t num_rows)
      : ValidRowTrackerForUnitTest(num_rows),
      RemappedDataReceiverForUnitTest(num_rows),
      TestDataProviderForRemapper<int>(num_rows) {}

    void set_GT_for_row_query_idx(const uint64_t row_query_idx, const std::vector<int>& data) {
      set_data_for_row_query_idx(row_query_idx, data);
    }
};

TEST_CASE("gt_remapper") {
  auto num_rows = 5u;
  ValidRowTrackerAndGTDataProviderForUnitTest tracker(num_rows);
  AllelesCombiner<ValidRowTrackerAndGTDataProviderForUnitTest> combiner(tracker, num_rows);
  GTRemapper<ValidRowTrackerAndGTDataProviderForUnitTest> gt_remapper(0u, tracker, combiner);
  std::string buffer_str;
  //Insert data
  std::vector<std::string> REF_vec(num_rows);
  std::vector<std::string> delimited_ALT_vec(num_rows);
  combiner.reset_before_adding_new_sample_info_at_current_position();
  //simple alleles - we're testing gt_remapper not alleles_combiner
  REF_vec[0u] = REF_vec[1u] = REF_vec[2u] = REF_vec[3u] = "A";
  delimited_ALT_vec[0u] = "T";
  delimited_ALT_vec[1u] = "G";
  delimited_ALT_vec[2u] = "C";
  delimited_ALT_vec[3u] = "CAA";
  insert_allele_info(tracker, combiner, 0u, REF_vec[0u], delimited_ALT_vec[0u]);
  insert_allele_info(tracker, combiner, 1u, REF_vec[1u], delimited_ALT_vec[1u]);
  insert_allele_info(tracker, combiner, 2u, REF_vec[2u], delimited_ALT_vec[2u]);
  insert_allele_info(tracker, combiner, 3u, REF_vec[3u], delimited_ALT_vec[3u]);
  //row_query_idx = 4 is invalid
  combiner.finished_updating_allele_info_for_current_position();
  tracker.clear_remapped_GT();
  tracker.set_GT_for_row_query_idx(0u, { 0, 0, 1}); // 0/1
  //<op, contains_phase, produce_GT_field, do_remap>
  gt_remapper.remap_all_queried_valid_rows<ValidRowTrackerAndGTDataProviderForUnitTest, true, true, true>(tracker);
  CHECK(tracker.get_remapped_GT(0u) == std::vector<int>({ 0, 0, 1}));
  //even without remapping, should return correct values for row 0 since its alleles have the same indexes in the
  //merged list
  tracker.clear_remapped_GT();
  gt_remapper.remap_all_queried_valid_rows<ValidRowTrackerAndGTDataProviderForUnitTest, true, true, false>(tracker);
  CHECK(tracker.get_remapped_GT(0u) == std::vector<int>({ 0, 0, 1}));
  //Add multiple samples
  tracker.clear_remapped_GT();
  tracker.set_GT_for_row_query_idx(0u, { 0, 0, 1}); // 0/1
  tracker.set_GT_for_row_query_idx(1u, { 0, 1, 1, 0, -1, 0, 1}); // 0|1/./1
  tracker.set_GT_for_row_query_idx(2u, { -1, 0, 1, 1, 1, 0, 0}); // ./1|1/0
  //<op, contains_phase, produce_GT_field, do_remap>
  gt_remapper.remap_all_queried_valid_rows<ValidRowTrackerAndGTDataProviderForUnitTest, true, true, true>(tracker);
  //REF=0, T=1, G=2, C=3
  CHECK(tracker.get_remapped_GT(0u) == std::vector<int>({ 0, 0, 1}));
  CHECK(tracker.get_remapped_GT(1u) == std::vector<int>({ 0, 1, 2, 0, -1, 0, 2}));
  CHECK(tracker.get_remapped_GT(2u) == std::vector<int>({ -1, 0, 3, 1, 3, 0, 0}));
  //GT for row_query_idx = 3 is missing
  CHECK(tracker.get_remapped_GT(3u) == std::vector<int>());
  CHECK(tracker.is_GT_missing(3u));
  //row_query_idx = 4 is invalid
  CHECK(tracker.get_remapped_GT(4u) == std::vector<int>());
  CHECK(!tracker.is_GT_missing(4u));
}
