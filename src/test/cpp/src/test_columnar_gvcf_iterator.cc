/**
 * src/test/cpp/src/test_columnar_gvcf_iterator.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation, Inc.
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
 * Test the GenomicsDB Columnar gvcf iterator
 *
 */

#include <catch2/catch.hpp>

#include "query_variants.h"
#include "variant_operations.h"
#include "variant_query_config.h"
#include "variant_storage_manager.h"
#include "test_common.h"
#include <htslib/vcf.h>

#include <iostream>

#include "test_valid_row_tracker.h"
#include "test_remapped_data_receiver.h"
#include "test_data_provider_for_remapper.h"
#include "alleles_combiner_template_definition.hpp"
#include "gt_remapper_template_definition.hpp"

static std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

static std::string query_json(ctests_input_dir+"query.json");
static std::string loader_json(ctests_input_dir+"loader.json");

static unsigned segment_size = 40;

TEST_CASE("gvcf iterator", "[gvcf_iterator]") {
  VariantQueryConfig query_config;
  GenomicsDBImportConfig loader_config;
  loader_config.read_from_file(loader_json);
  query_config.update_from_loader(loader_config);
  query_config.read_from_file(query_json);

  CHECK(query_config.get_num_column_intervals() == 1);
  
  VariantStorageManager storage_manager(query_config.get_workspace(0), query_config.get_segment_size());
  
  VariantQueryProcessor query_processor(&storage_manager, query_config.get_array_name(0), query_config.get_vid_mapper());
  query_processor.do_query_bookkeeping(query_processor.get_array_schema(), query_config, query_config.get_vid_mapper(), true);
  
  ProfilerOperator counter;
  query_processor.iterate_over_gvcf_entries(query_processor.get_array_descriptor(), query_config, counter, true);
  CHECK(counter.get_value() == 4);

  ProfilerOperator *variant_counter = new ProfilerOperator();
  query_processor.scan_and_operate(query_processor.get_array_descriptor(), query_config, *variant_counter, 0, true);
  CHECK(variant_counter->get_value() == 4);
  delete variant_counter;
}

//Holds a bunch of objects that can be passed as argument to functions
struct BagOfObjects
{
  bcf_hdr_t* hdr;
  bcf1_t* rec;
  int* vcf_buffer_ptr;
  int vcf_buffer_length;
  GenomicsDBGVCFIterator* columnar_gvcf_iter;
  std::vector<size_t> vcf_sample_idx_to_row_query_idx;
  unsigned contains_phase;
  unsigned produce_GT_field;
};

enum GoldTestRowQueryIdxEnum {
  GOLD=0u,
  TEST
};

class IntDataProviderAndReceiver : public RemappedDataReceiverForUnitTest, public TestDataProviderForRemapper<int>,
  public ValidRowTrackerForUnitTest
{
  public:
    IntDataProviderAndReceiver(const size_t num_rows)
      : RemappedDataReceiverForUnitTest(num_rows),
      TestDataProviderForRemapper(num_rows),
      ValidRowTrackerForUnitTest(num_rows) {}
};

struct BagOfRemappers
{
  public:
    BagOfRemappers()
      : m_remap_receiver(),
      m_int_handler(2u), //2 - gold and test
      m_gold_test_alleles_combiner(m_int_handler, 2u), //2 - gold and test
      m_test_gold_gt_remapper(0u, m_int_handler, m_gold_test_alleles_combiner) {
        m_int_handler.set_valid_row_query_idx(0u, true);
        m_int_handler.set_valid_row_query_idx(1u, true);
      }
    SingleRowRemappedDataReceiverForUnitTest m_remap_receiver; //receives remapped data for GVCF iterator
    //Provides/receives data to/from remapper which maps data from test to gold
    IntDataProviderAndReceiver m_int_handler;
    //Used to remap alleles from test to gold
    AllelesCombiner<IntDataProviderAndReceiver> m_gold_test_alleles_combiner;
    //GT remapper from test to gold
    GTRemapper<IntDataProviderAndReceiver> m_test_gold_gt_remapper;
};

//<op, contains_phase, produce_GT_field, do_remap>
template<bool contains_phase, bool produce_GT_field>
void do_GT_remap_to_gold(BagOfRemappers& bag_of_remappers, BagOfObjects& bag, const uint64_t row_query_idx)
{
  //First compute GT in test (from iterator) and receive the remapped GT in m_remap_receiver
  bag_of_remappers.m_remap_receiver.clear_remapped_GT();
  bag.columnar_gvcf_iter->get_GT_remapper().remap_for_row_query_idx
    <SingleRowRemappedDataReceiverForUnitTest, contains_phase, produce_GT_field, true>
    (bag_of_remappers.m_remap_receiver, row_query_idx);
  //Now remap from test to gold - m_int_handler is both a data provider and remapped data receiver
  bag_of_remappers.m_int_handler.set_data_for_row_query_idx(GoldTestRowQueryIdxEnum::TEST,
      bag_of_remappers.m_remap_receiver.get_remapped_GT(0u));
  bag_of_remappers.m_int_handler.clear_remapped_GT();
  bag_of_remappers.m_test_gold_gt_remapper.remap_for_row_query_idx
    <IntDataProviderAndReceiver, contains_phase, produce_GT_field, true>
    (bag_of_remappers.m_int_handler, GoldTestRowQueryIdxEnum::TEST);
  //Now m_int_handler.get_remapped_GT(TEST) should have GT values remapped to gold alleles order
}

template<typename DataType>
void compare_field(BagOfRemappers& bag_of_remappers, BagOfObjects& bag, const std::string& vcf_field_name,
    const unsigned genomicsdb_query_field_idx)
{
  auto num_output = bcf_get_format_int32(bag.hdr, bag.rec, vcf_field_name.c_str(), &(bag.vcf_buffer_ptr), &(bag.vcf_buffer_length));
  auto columnar_gvcf_iter = bag.columnar_gvcf_iter;
  REQUIRE(((num_output == -3) || (num_output >= 0)));
  const auto is_GT_field = (vcf_field_name == "GT");
  for(auto i=0;i<bcf_hdr_nsamples(bag.hdr);++i) {
    auto row_query_idx = bag.vcf_sample_idx_to_row_query_idx[i];
    if(num_output == -3) {
      //No data for field in gold VCF, so either the row is invalid or the field is invalid
      REQUIRE((!columnar_gvcf_iter->is_valid(row_query_idx)
            || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx)));
      continue;
    }
    auto num_per_sample = num_output/bcf_hdr_nsamples(bag.hdr);
    const auto base_ptr = bag.vcf_buffer_ptr + i*num_per_sample;
    auto one_valid = false;
    auto vcf_GT_is_dot = false;
    if(is_GT_field) {
      //You can't just drop the GT field in a VCF - you'll at least have a '.' for a sample that doesn't
      //have any data for a given position. The '.' gets converted to a bcf_no_call,bcf_vector_end_int32,..
      //by htslib
      //Yes, the VCF format is terrible
      for(auto j=0;j<num_per_sample;++j) {
        auto val = base_ptr[j];
        if(j == 0)
          vcf_GT_is_dot = (val == bcf_gt_missing);
        else
          vcf_GT_is_dot = vcf_GT_is_dot && (val == get_bcf_vector_end_value<int>());
      }
      one_valid = true;
    }
    else {
      //Check if field has one valid data item
      one_valid = !((base_ptr[0] == get_bcf_missing_value<DataType>())
        && (num_per_sample == 1 || (base_ptr[1] == get_bcf_vector_end_value<DataType>())));
      //for(auto j=0;j<num_per_sample;++j)
      //one_valid = one_valid || is_bcf_valid_value<DataType>(base_ptr[j]);
    }
    if(one_valid) {
      if(!is_GT_field || !vcf_GT_is_dot) {
        REQUIRE(columnar_gvcf_iter->is_valid(row_query_idx));
        REQUIRE(columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx));
      }
      else {
        //GT and gold VCF GT is dot - it's ok for GenomicsDB to have invalid data
        if(!columnar_gvcf_iter->is_valid(row_query_idx)
            || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx))
          return;
        //However, if data is imported from a multi-sample VCF into GenomicsDB, GenomicsDB would have imported a
        //'.' for the GT, ie, a single element with bcf_gt_no_call. This will be checked below
        //Yes, the VCF format is terrible
      }
      if(is_GT_field) {
        //<op, contains_phase, produce_GT_field, do_remap>
        switch((bag.contains_phase << 1u) | bag.produce_GT_field) {
          case 0u:
            do_GT_remap_to_gold<false, false>(bag_of_remappers, bag, row_query_idx);
            break;
          case 1u:
            do_GT_remap_to_gold<false, true>(bag_of_remappers, bag, row_query_idx);
            break;
          case 2u:
            do_GT_remap_to_gold<true, false>(bag_of_remappers, bag, row_query_idx);
            break;
          case 3u:
            do_GT_remap_to_gold<true, true>(bag_of_remappers, bag, row_query_idx);
            break;
        }
      }
      auto ptr_length_pair = columnar_gvcf_iter->get_raw_pointer_and_length_for_query_idx(row_query_idx,
          genomicsdb_query_field_idx);
      auto& test_remapped_GT_vec = bag_of_remappers.m_int_handler.get_remapped_GT(GoldTestRowQueryIdxEnum::TEST);
      for(auto j=0,test_j=0;j<num_per_sample;++j) {
        auto val = base_ptr[j];
        if(is_bcf_vector_end_value<DataType>(val))
          break;
        if(is_GT_field) {
          //Phase information starts from 2nd element in VCF GT
          if(j > 0) {
            REQUIRE(static_cast<size_t>(test_j) < test_remapped_GT_vec.size());
            //If the schema, doesn't store phase with GT - then the query operation
            //always returns no phase
            if(bag.contains_phase)
              CHECK(bcf_gt_is_phased(val) == test_remapped_GT_vec[test_j]);
            else
              CHECK(test_remapped_GT_vec[test_j] == 0);
            ++test_j;
          }
          REQUIRE(static_cast<size_t>(test_j) < test_remapped_GT_vec.size());
          if(g_skip_GT_matching && bcf_gt_allele(val) != test_remapped_GT_vec[test_j]) {
            WARN("GT mismatch for sample " << bcf_hdr_int2id(bag.hdr, BCF_DT_SAMPLE, i)
                << " at index " << j << " gold " << bcf_gt_allele(val) << " test "
                << test_remapped_GT_vec[test_j] << " - PL remapping is WIP");
          }
          else {
            CHECK(bcf_gt_allele(val) == test_remapped_GT_vec[test_j]);
          }
          ++test_j;
        }
        else {
          REQUIRE(static_cast<unsigned>(j) < ptr_length_pair.second);
          CHECK(reinterpret_cast<const DataType*>(ptr_length_pair.first)[j] == val);
        }
      }
    }
    else
      REQUIRE((!columnar_gvcf_iter->is_valid(row_query_idx)
            || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx)));
  }
}

TEST_CASE("columnar_gvcf_iterator_test", "[gvcf_iterator]") {
  if(g_query_json_file.empty() || g_golden_output_file.empty())
    return;
  VariantQueryConfig query_config;
  GenomicsDBImportConfig loader_config;
  if(!g_loader_json_file.empty()) {
    loader_config.read_from_file(g_loader_json_file);
    query_config.update_from_loader(loader_config);
  }
  query_config.read_from_file(g_query_json_file);

  VariantStorageManager storage_manager(query_config.get_workspace(0), query_config.get_segment_size());

  VariantQueryProcessor query_processor(&storage_manager, query_config.get_array_name(0), query_config.get_vid_mapper());
  query_processor.do_query_bookkeeping(query_processor.get_array_schema(), query_config, query_config.get_vid_mapper(), true);
  auto& vid_mapper = query_config.get_vid_mapper();
  BagOfRemappers bag_of_remappers;
  BagOfObjects bag;
  bag.columnar_gvcf_iter = storage_manager.begin_gvcf_iterator(
      query_processor.get_array_descriptor(), query_config,
      true);
  REQUIRE(bag.columnar_gvcf_iter != 0);
  bag.contains_phase = query_config.is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX)
    && query_config.get_length_descriptor_for_query_attribute_idx(
        query_config.get_query_idx_for_known_field_enum(GVCF_GT_IDX)).contains_phase_information();
  bag.produce_GT_field = query_config.produce_GT_field();

  //Golden VCF
  auto fptr = bcf_open(g_golden_output_file.c_str(), "rb");
  REQUIRE(fptr != 0);
  bag.hdr = bcf_hdr_read(fptr);
  REQUIRE(bag.hdr != 0);
  if(!query_config.sites_only_query())
  {
    REQUIRE(bcf_hdr_nsamples(bag.hdr) == query_config.get_num_rows_to_query());
    //Map between samples in gold hdr and test output
    bag.vcf_sample_idx_to_row_query_idx.resize(bcf_hdr_nsamples(bag.hdr));
    for(auto i=0;i<bcf_hdr_nsamples(bag.hdr);++i)
    {
      int64_t row_idx = 0;
      REQUIRE(vid_mapper.get_tiledb_row_idx(row_idx, bcf_hdr_int2id(bag.hdr, BCF_DT_SAMPLE, i)));
      bag.vcf_sample_idx_to_row_query_idx[i] = query_config.get_query_row_idx_for_array_row_idx(row_idx);
    }
  }
  bag.rec = bcf_init();
  //Check SB and GQ fields since they have no allele dependence and are unmodified from input to output
  unsigned SB_query_idx = 0u;
  auto is_SB_queried = query_config.get_query_idx_for_name("SB", SB_query_idx);
  unsigned GQ_query_idx = 0u;
  auto is_GQ_queried = query_config.get_query_idx_for_name("GQ", GQ_query_idx);
  unsigned GT_query_idx = 0u;
  auto is_GT_queried = query_config.get_query_idx_for_name("GT", GT_query_idx);
  //VCF info for SB, GQ, GT
  auto do_check_SB = is_SB_queried && (bcf_hdr_id2int(bag.hdr, BCF_DT_ID, "SB") >= 0);
  auto do_check_GQ = is_GQ_queried && (bcf_hdr_id2int(bag.hdr, BCF_DT_ID, "GQ") >= 0);
  auto do_check_GT = is_GT_queried && (bcf_hdr_id2int(bag.hdr, BCF_DT_ID, "GT") >= 0);
  //Avoid some memory reallocations
  bag.vcf_buffer_length = 4*bcf_hdr_nsamples(bag.hdr); //since SB has 4 elements per sample - will get resized if necessary
  bag.vcf_buffer_ptr = new int[bag.vcf_buffer_length]; //use single buffer
  std::string contig_name;
  //REF-ALT comparison
  std::vector<STRING_VIEW> gold_alleles;
  std::vector<STRING_VIEW> test_alleles;
  std::string alleles_buffer;
  std::string tmp_buffer;
  std::vector<STRING_VIEW> tmp_alleles;
  //Using AllelesCombiner to keep track of mapping between gold alleles and test allelesi
  //Will be used for GT remapping etc
  //Mark both gold and test as valid samples for remapping
  bag_of_remappers.m_int_handler.set_valid_row_query_idx(GoldTestRowQueryIdxEnum::GOLD, true);
  bag_of_remappers.m_int_handler.set_valid_row_query_idx(GoldTestRowQueryIdxEnum::TEST, true);
  for (; !(bag.columnar_gvcf_iter->end()); ++(*(bag.columnar_gvcf_iter))) {
    auto column_interval = bag.columnar_gvcf_iter->get_current_variant_interval();
    REQUIRE(bcf_read(fptr, bag.hdr, bag.rec) == 0);
    int64_t contig_position = -1;
    REQUIRE(vid_mapper.get_contig_location(column_interval.first, contig_name, contig_position));
    CHECK(contig_name == std::string(bcf_seqname(bag.hdr, bag.rec)));
    CHECK(contig_position == bag.rec->pos);
    //std::cerr << rec->pos <<"\n";
    bcf_unpack(bag.rec, BCF_UN_STR);
    //REF-ALT comparison
    gold_alleles.clear();
    for(auto i=0u;i<bag.rec->n_allele;++i)
      gold_alleles.emplace_back(STRING_VIEW(bag.rec->d.allele[i], strlen(bag.rec->d.allele[i])));
    alleles_buffer.clear();
    test_alleles.clear();
    bag.columnar_gvcf_iter->get_alleles_combiner().get_merged_VCF_spec_alleles_vec(alleles_buffer, test_alleles);
    if(test_alleles[0].length() == 0u) //leftover REF block remains
      test_alleles[0] = STRING_VIEW(bag.rec->d.allele[0], strlen(bag.rec->d.allele[0]));
    //Must have same alleles - possibly in a different order
    REQUIRE(gold_alleles[0u] == test_alleles[0u]);
    REQUIRE(std::unordered_set<STRING_VIEW>(gold_alleles.begin(), gold_alleles.end()) ==
        std::unordered_set<STRING_VIEW>(test_alleles.begin(), test_alleles.end()));
    //Find mapping from test alleles to gold alleles using AllelesCombiner
    bag_of_remappers.m_gold_test_alleles_combiner.remove_allele_info(GoldTestRowQueryIdxEnum::GOLD);
    bag_of_remappers.m_gold_test_alleles_combiner.remove_allele_info(GoldTestRowQueryIdxEnum::TEST);
    bag_of_remappers.m_gold_test_alleles_combiner.reset_before_adding_new_sample_info_at_current_position();
    bag_of_remappers.m_gold_test_alleles_combiner.insert_allele_info(GoldTestRowQueryIdxEnum::GOLD, gold_alleles[0],
        std::vector<STRING_VIEW>(gold_alleles.begin()+1u, gold_alleles.end()), false);
    bag_of_remappers.m_gold_test_alleles_combiner.insert_allele_info(GoldTestRowQueryIdxEnum::TEST, test_alleles[0],
        std::vector<STRING_VIEW>(test_alleles.begin()+1u, test_alleles.end()), false);
    bag_of_remappers.m_gold_test_alleles_combiner.finished_updating_allele_info_for_current_position();
    tmp_buffer.clear();
    tmp_alleles.clear();
    bag_of_remappers.m_gold_test_alleles_combiner.get_merged_VCF_spec_alleles_vec(tmp_buffer, tmp_alleles);
    //Merged allele order in bag_of_remappers.m_gold_test_alleles_combiner should be same as gold vcf
    REQUIRE(gold_alleles == tmp_alleles);
    //Check if mapping from test to gold is correct
    for(auto i=0u;i<test_alleles.size();++i) {
      auto merged_idx = bag_of_remappers.m_gold_test_alleles_combiner.get_merged_allele_idx(GoldTestRowQueryIdxEnum::TEST, i);
      REQUIRE(test_alleles[i] == gold_alleles[merged_idx]);
    }
    if(query_config.sites_only_query())
      continue;
    //Check SB and GQ fields since they have no allele dependence and are unmodified from input to output
    if(do_check_SB)
      compare_field<int>(bag_of_remappers, bag, "SB", SB_query_idx);
    if(do_check_GQ)
      compare_field<int>(bag_of_remappers, bag, "GQ", GQ_query_idx);
    if(do_check_GT)
      compare_field<int>(bag_of_remappers, bag, "GT", GT_query_idx);
  }
  CHECK(bcf_read(fptr, bag.hdr, bag.rec) == -1); //no more records
  delete[] bag.vcf_buffer_ptr;
  bcf_destroy(bag.rec);
  bcf_hdr_destroy(bag.hdr);
  bcf_close(fptr);
  delete bag.columnar_gvcf_iter;
}
