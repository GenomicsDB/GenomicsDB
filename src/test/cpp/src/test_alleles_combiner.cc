/**
 * src/test/cpp/src/test_alleles_combiner.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2020-21 Omics Data Automation, Inc.
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
 * Test the optimized alleles combiner 
 *
 */

#include <catch2/catch.hpp>
#include "alleles_combiner_template_definition.h"
#include "test_valid_row_tracker.h"

TEST_CASE("alleles_combiner") {
  auto num_rows = 3u;
  ValidRowTrackerForUnitTest tracker(num_rows);
  AllelesCombiner<ValidRowTrackerForUnitTest> combiner(tracker, num_rows);
  std::string buffer_str;
  //Convenience vars
  const auto& row_query_idx_with_deletion_MNV_added_in_current_iteration_vec = combiner.get_row_query_idx_with_deletions_or_MNVs_at_current_location_vec();
  const auto& merged_vec = combiner.get_merged_normalized_REF_ALT_vec();
  const auto& lut = combiner.get_merged_alleles_lut();
  //Insert data
  std::vector<std::string> REF_vec(num_rows);
  std::vector<std::string> delimited_ALT_vec(num_rows);
  REF_vec[2u] = "ATGC";
  delimited_ALT_vec[2u] = "A|AGC|TTGC|*|TAGC";
  insert_allele_info(tracker, combiner, 2u, REF_vec[2u], delimited_ALT_vec[2u]);
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 1u);
  CHECK(combiner.get_num_calls_with_NON_REF_allele() == 0u);
  CHECK(combiner.contains_deletion_or_MNV(2u));
  CHECK(!combiner.contains_NON_REF_allele(2u));
  REQUIRE(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec.size() == 1u); //1 sample has deletion/MNV
  CHECK(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec[0u] == 2u);      //row query idx is 2
  auto ptr_length_pair = combiner.get_deletion_MNV_indexes_for_row_query_idx_at_index(0u);
  CHECK(ptr_length_pair.second == 4u);   //row query idx 2 has 4 deletions/MNVs
  auto ptr = ptr_length_pair.first;
  CHECK(ptr[0u] == 1u); //ATGC->A REF is allele idx 0
  CHECK(ptr[1u] == 2u); //AT->A
  CHECK(ptr[2u] == 4u); //A->*
  CHECK(ptr[3u] == 5u); //AT->TA (MNV)
  REQUIRE(merged_vec.size() == 6u); //REF and ALT
  //REF
  CHECK(merged_vec[0u].allele == REF_vec[2u]);
  //ATGC->A
  CHECK(merged_vec[1u].REF_length == REF_vec[2u].length());
  CHECK(merged_vec[1u].allele == "A");
  CHECK(!merged_vec[1u].is_symbolic_allele);
  //AT->A
  CHECK(merged_vec[2u].REF_length == 2u);
  CHECK(merged_vec[2u].allele == "A");
  CHECK(!merged_vec[2u].is_symbolic_allele);
  //A->T
  CHECK(merged_vec[3u].REF_length == 1u);
  CHECK(merged_vec[3u].allele == "T");
  CHECK(!merged_vec[3u].is_symbolic_allele);
  //A->*
  CHECK(merged_vec[4u].REF_length == 1u);
  CHECK(merged_vec[4u].allele == "*");
  CHECK(merged_vec[4u].is_symbolic_allele);
  //AT->TA
  CHECK(merged_vec[5u].REF_length == 2u);
  CHECK(merged_vec[5u].allele == "TA");
  CHECK(!merged_vec[5u].is_symbolic_allele);
  //Throw in another sample with the same set of alleles + NON_REF but in a different order
  REF_vec[1u] = "ATGC";
  delimited_ALT_vec[1u] = "TTGC|*|A|TAGC|&|AGC";
  insert_allele_info(tracker, combiner, 1u, REF_vec[1u], delimited_ALT_vec[1u]);
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 2u);
  CHECK(combiner.get_num_calls_with_NON_REF_allele() == 1u);
  CHECK(combiner.contains_deletion_or_MNV(1u));
  CHECK(combiner.contains_NON_REF_allele(1u));
  CHECK(merged_vec.size() == 6u); //since NON_REF alleles aren't added to merged_vec
  REQUIRE(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec.size() == 2u); //one more sample has deletion/MNV
  CHECK(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec[1u] == 1u); //row query idx 1u
  ptr_length_pair = combiner.get_deletion_MNV_indexes_for_row_query_idx_at_index(1u); //get deletion indexes for row query idx 1u
  CHECK(ptr_length_pair.second == 4u);   //row query idx 1 has 4 deletions/MNVs
  ptr = ptr_length_pair.first;
  CHECK(ptr[0u] == 2u); //A->*, REF is allele idx 0
  CHECK(ptr[1u] == 3u); //ATGC->A
  CHECK(ptr[2u] == 4u); //AT->TA (MNV)
  CHECK(ptr[3u] == 6u); //AT->A
  //throw in one more sample
  REF_vec[0u] = "A";
  delimited_ALT_vec[0u] = "T|<SYM>|G|AGGT";
  insert_allele_info(tracker, combiner, 0u, REF_vec[0u], delimited_ALT_vec[0u]);
  //No changes to NON_REF or deletion counters
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 2u);
  CHECK(combiner.get_num_calls_with_NON_REF_allele() == 1u);
  CHECK(!combiner.contains_deletion_or_MNV(0u));
  CHECK(!combiner.contains_NON_REF_allele(0u));
  REQUIRE(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec.size() == 2u); //no samples with deletion/MNV added
  REQUIRE(merged_vec.size() == 9u); //3 new alt alleles added
  //A-><SYM>
  CHECK(merged_vec[6u].REF_length ==1u);
  CHECK(merged_vec[6u].allele == "<SYM>");
  CHECK(merged_vec[6u].is_symbolic_allele);
  //A->G
  CHECK(merged_vec[7u].REF_length == 1u);
  CHECK(merged_vec[7u].allele == "G");
  CHECK(!merged_vec[7u].is_symbolic_allele);
  //A->AGGT
  CHECK(merged_vec[8u].REF_length == 1u);
  CHECK(merged_vec[8u].allele == "AGGT");
  CHECK(!merged_vec[8u].is_symbolic_allele);
  //Check mapping
  //For row query idx 2 "A|AGC|TTGC|*|TAGC";
  CHECK(lut.get_merged_idx_for_input(2u, 0u) == 0u);
  CHECK(lut.get_merged_idx_for_input(2u, 1u) == 1u);
  CHECK(lut.get_merged_idx_for_input(2u, 2u) == 2u);
  CHECK(lut.get_merged_idx_for_input(2u, 3u) == 3u);
  CHECK(lut.get_merged_idx_for_input(2u, 4u) == 4u);
  CHECK(lut.get_merged_idx_for_input(2u, 5u) == 5u);
  //Must be invalid mapping for other alleles in merged list
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(2u, 6u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(2u, 7u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(2u, 8u)));
  //For row query idx 1 "TTGC|*|A|TAGC|&|AGC"
  CHECK(lut.get_merged_idx_for_input(1u, 0u) == 0u);
  CHECK(lut.get_merged_idx_for_input(1u, 1u) == 3u);
  CHECK(lut.get_merged_idx_for_input(1u, 2u) == 4u);
  CHECK(lut.get_merged_idx_for_input(1u, 3u) == 1u);
  CHECK(lut.get_merged_idx_for_input(1u, 4u) == 5u);
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_merged_idx_for_input(1u, 5u))); //NON_REF isn't mapped to LUT
  CHECK(lut.get_merged_idx_for_input(1u, 6u) == 2u);
  //Must be invalid mapping for other alleles in merged list
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(1u, 6u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(1u, 7u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(1u, 8u)));
  //row query idx 0 "T|<SYM>|G|AGGT"
  CHECK(lut.get_merged_idx_for_input(0u, 0u) == 0u);
  CHECK(lut.get_merged_idx_for_input(0u, 1u) == 3u);
  CHECK(lut.get_merged_idx_for_input(0u, 2u) == 6u);
  CHECK(lut.get_merged_idx_for_input(0u, 3u) == 7u);
  CHECK(lut.get_merged_idx_for_input(0u, 4u) == 8u);
  //Must be invalid mapping for other alleles in merged list
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(0u, 1u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(0u, 2u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(0u, 4u)));
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(0u, 5u)));
  //Done with all allele additions
  combiner.finished_updating_allele_info_for_current_position();
  CHECK(merged_vec.size() == 10u); //NON_REF appended
  {
    //Check if merged ALT alleles (de-normalized as per VCF spec) is correct
    //For row query idx 2 "A|AGC|TTGC|*|TAGC";
    //For row query idx 1 "TTGC|*|A|TAGC|&|AGC"
    //row query idx 0 "T|<SYM>|G|AGGT"
    std::vector<STRING_VIEW> gold_alleles = { "ATGC","A","AGC","TTGC","*","TAGC","<SYM>","GTGC","AGGTTGC","<NON_REF>" };
    std::vector<STRING_VIEW> alleles_vec;
    combiner.get_merged_VCF_spec_alleles_vec(buffer_str, alleles_vec);
    REQUIRE(gold_alleles.size() == alleles_vec.size());
    for(auto i=0u;i<gold_alleles.size();++i)
      CHECK(gold_alleles[i] == alleles_vec[i]);
    CHECK(combiner.get_NON_REF_allele_index_in_merged_allele_list() == merged_vec.size()-1u);
    CHECK(combiner.get_spanning_deletion_allele_index_in_merged_allele_list() == 4u);
  }
  //Move on to next iteration
  remove_allele_info(tracker, combiner, 0u);
  //samples 1 and 2 have spanning deletions and NON_REF
  combiner.reset_before_adding_new_sample_info_at_current_position();
  CHECK(combiner.get_num_calls_with_NON_REF_allele() == 1u);
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 2u);
  REQUIRE(merged_vec.size() == 1u); //must have reset, just a position for REF
  //No new samples added in this iteration
  combiner.finished_updating_allele_info_for_current_position();
  CHECK(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec.size() == 0u); //no new deletions/MNVs added in this iteration
  //Merged list should have just REF + 2 alt alleles - spanning_deletion and NON_REF
  REQUIRE(merged_vec.size() == 3u);
  {
    //REF is empty since it's the middle of a variant and we're not getting the base from the reference genome
    std::vector<STRING_VIEW> gold_alleles = { "", "*","<NON_REF>" };
    std::vector<STRING_VIEW> alleles_vec;
    combiner.get_merged_VCF_spec_alleles_vec(buffer_str, alleles_vec);
    REQUIRE(gold_alleles.size() == alleles_vec.size());
    for(auto i=0u;i<gold_alleles.size();++i)
      CHECK(gold_alleles[i] == alleles_vec[i]);
    CHECK(combiner.get_NON_REF_allele_index_in_merged_allele_list() == merged_vec.size()-1u);
    CHECK(combiner.get_spanning_deletion_allele_index_in_merged_allele_list() == 1u);
  }
  //Move on to next iteration
  //Remove row query idx 1u
  remove_allele_info(tracker, combiner, 1u);
  CHECK(combiner.get_num_calls_with_NON_REF_allele() == 0u); //reduced by 1
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 1u); //reduced by 1
  combiner.reset_before_adding_new_sample_info_at_current_position();
  REQUIRE(merged_vec.size() == 1u); //must have reset, just a position for REF
  //No new samples added in this iteration
  combiner.finished_updating_allele_info_for_current_position();
  CHECK(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec.size() == 0u);
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 1u); //remains same - row_query_idx=2 should still be counted
  //Merged list should have just REF + 1 alt alleles - spanning_deletion
  REQUIRE(merged_vec.size() == 2u);
  {
    //REF is empty since it's the middle of a variant and we're not getting the base from the reference genome
    std::vector<STRING_VIEW> gold_alleles = { "", "*" };
    std::vector<STRING_VIEW> alleles_vec;
    combiner.get_merged_VCF_spec_alleles_vec(buffer_str, alleles_vec);
    REQUIRE(gold_alleles.size() == alleles_vec.size());
    for(auto i=0u;i<gold_alleles.size();++i)
      CHECK(gold_alleles[i] == alleles_vec[i]);
    CHECK(!combiner.merged_alleles_list_contains_NON_REF());
    CHECK(combiner.get_spanning_deletion_allele_index_in_merged_allele_list() == 1u);
  }
  //Move on to next iteration
  remove_allele_info(tracker, combiner, 2u);
  combiner.reset_before_adding_new_sample_info_at_current_position();
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 0u);
  //New set of alleles
  REF_vec[0u] = "A";
  delimited_ALT_vec[0u] = "T|AGGT";
  REF_vec[1u] = "ATGC";
  delimited_ALT_vec[1u] = "A|AGGTTGC";
  REF_vec[2u] = "ATGC";
  delimited_ALT_vec[2u] = "TTGC|A";
  insert_allele_info(tracker, combiner, 0u, REF_vec[0u], delimited_ALT_vec[0u]);
  insert_allele_info(tracker, combiner, 1u, REF_vec[1u], delimited_ALT_vec[1u]);
  insert_allele_info(tracker, combiner, 2u, REF_vec[2u], delimited_ALT_vec[2u]);
  combiner.finished_updating_allele_info_for_current_position();
  CHECK(combiner.get_num_calls_with_deletions_or_MNVs() == 2u);
  CHECK(combiner.get_num_calls_with_NON_REF_allele() == 0u);
  for(auto i=0u;i<3u;++i)
    CHECK(!combiner.contains_NON_REF_allele(i));
  CHECK(!combiner.contains_deletion_or_MNV(0u));
  CHECK(combiner.contains_deletion_or_MNV(1u));
  CHECK(combiner.contains_deletion_or_MNV(2u));
  REQUIRE(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec.size() == 2u); //2 samples have deletion/MNV
  CHECK(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec[0u] == 1u); //row query idxs 1 and 2 have deletions
  CHECK(row_query_idx_with_deletion_MNV_added_in_current_iteration_vec[1u] == 2u);
  ptr_length_pair = combiner.get_deletion_MNV_indexes_for_row_query_idx_at_index(0u);
  CHECK(ptr_length_pair.second == 1u);   //row query idx 1 has 1 deletions/MNVs
  CHECK(ptr_length_pair.first[0u] == 1u); //deletion is at index 1
  ptr_length_pair = combiner.get_deletion_MNV_indexes_for_row_query_idx_at_index(1u);
  CHECK(ptr_length_pair.second == 1u);   //row query idx 2 has 1 deletions/MNVs
  CHECK(ptr_length_pair.first[0u] == 2u); //deletion is at index 2
  REQUIRE(merged_vec.size() == 4u); //REF and ALT
  //REF
  CHECK(merged_vec[0u].allele == "ATGC");
  //A->T
  CHECK(merged_vec[1u].REF_length == 1u);
  CHECK(merged_vec[1u].allele == "T");
  CHECK(!merged_vec[1u].is_symbolic_allele);
  //A->AGGT
  CHECK(merged_vec[2u].REF_length == 1u);
  CHECK(merged_vec[2u].allele == "AGGT");
  CHECK(!merged_vec[2u].is_symbolic_allele);
  //ATGC->A
  CHECK(merged_vec[3u].REF_length == 4u);
  CHECK(merged_vec[3u].allele == "A");
  CHECK(!merged_vec[3u].is_symbolic_allele);
  //Check mapping
  CHECK(lut.get_merged_idx_for_input(0u, 0u) == 0u); //REF
  CHECK(lut.get_merged_idx_for_input(1u, 0u) == 0u); //REF
  CHECK(lut.get_merged_idx_for_input(2u, 0u) == 0u); //REF
  //For row query idx 0 "T|AGGT";
  CHECK(lut.get_merged_idx_for_input(0u, 1u) == 1u);
  CHECK(lut.get_merged_idx_for_input(0u, 2u) == 2u);
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(0u, 3u)));
  //For row query idx 1, "A|AGGTTGC"
  CHECK(lut.get_merged_idx_for_input(1u, 1u) == 3u);
  CHECK(lut.get_merged_idx_for_input(1u, 2u) == 2u);
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(1u, 1u)));
  //For row query idx 2, "TTGC|A"
  CHECK(lut.get_merged_idx_for_input(2u, 1u) == 1u);
  CHECK(lut.get_merged_idx_for_input(2u, 2u) == 3u);
  CHECK(CombineAllelesLUT::is_missing_value(lut.get_input_idx_for_merged(2u, 2u)));
  {
    //Check if merged ALT alleles (de-normalized as per VCF spec) is correct
    //For row query idx 0 "T|AGGT";
    //For row query idx 1, "A|AGGTTGC"
    //For row query idx 2, "TTGC|A"
    std::vector<STRING_VIEW> gold_alleles = { "ATGC","TTGC","AGGTTGC","A" };
    std::vector<STRING_VIEW> alleles_vec;
    combiner.get_merged_VCF_spec_alleles_vec(buffer_str, alleles_vec);
    REQUIRE(gold_alleles.size() == alleles_vec.size());
    for(auto i=0u;i<gold_alleles.size();++i)
      CHECK(gold_alleles[i] == alleles_vec[i]);
    CHECK(!combiner.merged_alleles_list_contains_NON_REF());
    CHECK(!combiner.merged_alleles_list_contains_spanning_deletion());
  }
}
