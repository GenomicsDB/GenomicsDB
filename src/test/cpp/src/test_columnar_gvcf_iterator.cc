/**
 * src/test/cpp/src/test_columnar_gvcf_iterator.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation, Inc.
 * Copyright (c) 2023 dātma, inc™
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions
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
#include "alleles_combiner_template_definition.h"
#include "gt_remapper_template_definition.h"
#include "broad_combined_gvcf.h"

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

enum GoldTestEnum {
  GOLD=0u,
  TEST
};

//Holds a bunch of objects that can be passed as argument to functions
struct BagOfObjects
{
  bcf_hdr_t* hdr[2u]; //gold and test
  bcf1_t* rec[2u];
  int* vcf_buffer_ptr[2u];
  int vcf_buffer_length[2u];
  GenomicsDBGVCFIterator* columnar_gvcf_iter;
  std::vector<size_t> vcf_sample_idx_to_row_query_idx;
  unsigned contains_phase;
  unsigned produce_GT_field;
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
  bag_of_remappers.m_int_handler.set_data_for_row_query_idx(GoldTestEnum::TEST,
      bag_of_remappers.m_remap_receiver.get_remapped_GT(0u));
  bag_of_remappers.m_int_handler.clear_remapped_GT();
  bag_of_remappers.m_test_gold_gt_remapper.remap_for_row_query_idx
    <IntDataProviderAndReceiver, contains_phase, produce_GT_field, true>
    (bag_of_remappers.m_int_handler, GoldTestEnum::TEST);
  //Now m_int_handler.get_remapped_GT(TEST) should have GT values remapped to gold alleles order
}

template<typename DataType>
bool check_if_field_has_one_valid_data_item(const DataType* sample_ptr, const int num_per_sample) {
  //Check if field has one valid data item
  //This is a problem for multi-sample VCF files. Every sample must have the same number of elements
  //in the FORMAT structure of htslib (similar to bcf). So, for samples which don't have any valid
  //value for a field, htslib inserts bcf_missing followed by bcf_vector_end (if needed). So, the field
  //is valid if the above condition is not true
  //Yes, the VCF format is terrible
  return !((num_per_sample <= 0) || ((is_bcf_missing_value<DataType>(sample_ptr[0]))
      && (num_per_sample == 1 || (is_bcf_vector_end_value<DataType>(sample_ptr[1])))));
}

//To keep the compile happy when DataType = float
template<typename DataType>
bool is_vcf_GT_dot(const DataType* sample_ptr, const int num_per_sample) {
  return false;
}

//template specialization for int - actual GT check
template<>
bool is_vcf_GT_dot(const int* sample_ptr, const int num_per_sample) {
  //You can't just drop the GT field in a VCF - you'll at least have a '.' for a sample that doesn't
  //have any data for a given position. The '.' gets converted to a bcf_no_call,bcf_vector_end_int32,..
  //by htslib
  //Yes, the VCF format is terrible
  auto vcf_GT_is_dot = false;
  for(auto j=0;j<num_per_sample;++j) {
    auto val = sample_ptr[j];
    if(j == 0)
      vcf_GT_is_dot = (val == bcf_gt_missing);
    else
      vcf_GT_is_dot = vcf_GT_is_dot && is_bcf_vector_end_value<int>(val);
  }
  return vcf_GT_is_dot;
}

template<typename DataType>
void compare_field(BagOfRemappers& bag_of_remappers, BagOfObjects& bag, const std::string& vcf_field_name,
    const unsigned genomicsdb_query_field_idx)
{
  auto columnar_gvcf_iter = bag.columnar_gvcf_iter;
  int num_output[2]; //gold and test
  num_output[GoldTestEnum::GOLD] = bcf_get_format_int32(bag.hdr[GoldTestEnum::GOLD], bag.rec[GoldTestEnum::GOLD],
      vcf_field_name.c_str(),
      &(bag.vcf_buffer_ptr[GoldTestEnum::GOLD]), &(bag.vcf_buffer_length[GoldTestEnum::GOLD]));
  num_output[GoldTestEnum::TEST] = bcf_get_format_int32(bag.hdr[GoldTestEnum::TEST], bag.rec[GoldTestEnum::TEST],
      vcf_field_name.c_str(),
      &(bag.vcf_buffer_ptr[GoldTestEnum::TEST]), &(bag.vcf_buffer_length[GoldTestEnum::TEST]));
  REQUIRE(((num_output[GoldTestEnum::GOLD] == -3) || (num_output[GoldTestEnum::GOLD] >= 0)));
  REQUIRE(((num_output[GoldTestEnum::TEST] == -3) || (num_output[GoldTestEnum::TEST] >= 0)));
  REQUIRE(num_output[GoldTestEnum::GOLD] == num_output[GoldTestEnum::TEST]);
  const auto is_GT_field = (vcf_field_name == "GT");
  const auto num_per_sample = bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]) ?
    num_output[GoldTestEnum::GOLD]/bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]) : 0;
  for(auto i=0;i<bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]);++i) {
    auto row_query_idx = bag.vcf_sample_idx_to_row_query_idx[i];
    if(num_output[GoldTestEnum::GOLD] == -3) {
      //No data for field in gold VCF, so either the row is invalid or the field is invalid
      REQUIRE((!columnar_gvcf_iter->is_valid(row_query_idx)
            || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx)));
      REQUIRE(num_output[GoldTestEnum::TEST] == -3);
      continue;
    }
    const DataType* base_ptr[2];
    base_ptr[GoldTestEnum::GOLD] = bag.vcf_buffer_ptr[GoldTestEnum::GOLD] + i*num_per_sample;
    base_ptr[GoldTestEnum::TEST] = bag.vcf_buffer_ptr[GoldTestEnum::TEST] + row_query_idx*num_per_sample;
    auto one_valid = false;
    auto vcf_GT_is_dot = false;
    if(is_GT_field) {
      vcf_GT_is_dot = is_vcf_GT_dot<DataType>(base_ptr[GoldTestEnum::GOLD], num_per_sample);
      one_valid = true;
    }
    else {
      one_valid = check_if_field_has_one_valid_data_item<DataType>(base_ptr[GoldTestEnum::GOLD], num_per_sample);
      //for(auto j=0;j<num_per_sample;++j)
      //one_valid = one_valid || is_bcf_valid_value<DataType>(base_ptr[GoldTestEnum::GOLD][j]);
    }
    if(!one_valid) {
      //Invalid data in GenomicsDB - that's good
      if((!columnar_gvcf_iter->is_valid(row_query_idx)
            || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx))) {
        //Same condition as if stmt - guaranteed to succeed (kept here because I like to see the number of assertions
        //that succeeded in the unit test
        REQUIRE((!columnar_gvcf_iter->is_valid(row_query_idx)
            || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx)));
        //Test vcf should have no valid data either
        REQUIRE((num_output[GoldTestEnum::TEST] == -3
            || !check_if_field_has_one_valid_data_item<DataType>(base_ptr[GoldTestEnum::TEST], num_per_sample)));
        continue;
      }
      //See problem with multi-sample VCF files - here the data is imported from a multi-sample file. So,
      //GenomicsDB thinks the data is valid even though it's filled with bcf_missing values. We let the
      //check play out
      //Yes, the VCF format is terrible
    }
    if(!is_GT_field || !vcf_GT_is_dot) {
      REQUIRE(columnar_gvcf_iter->is_valid(row_query_idx));
      REQUIRE(columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx));
    }
    else {
      //GT and gold VCF GT is dot - it's ok for GenomicsDB to have invalid data
      if(!columnar_gvcf_iter->is_valid(row_query_idx)
          || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(genomicsdb_query_field_idx, row_query_idx)) {
        //Test vcf line must have no data or GT = . (same as gold)
        REQUIRE((num_output[GoldTestEnum::TEST] <= 0 || is_vcf_GT_dot<DataType>(base_ptr[GoldTestEnum::TEST], num_per_sample)));
        return;
      }
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
    //This vector contains the GT values produced by GenomicsDBGVCFIterator and AllelesCombiner/GTRemapper for
    //the queried samples
    auto& test_GT_vec = bag_of_remappers.m_remap_receiver.get_remapped_GT(GoldTestEnum::TEST);
    //This vector contains the GT values produced by running GTRemapper on test_GT_vec to map allele
    //indexes from test order to gold order
    auto& test_GT_vec_remapped_to_gold = bag_of_remappers.m_int_handler.get_remapped_GT(GoldTestEnum::TEST);
    auto hit_vector_end = false;
    for(auto j=0,test_j=0;j<num_per_sample;++j) {
      DataType val[2];
      val[GoldTestEnum::GOLD] = base_ptr[GoldTestEnum::GOLD][j];
      val[GoldTestEnum::TEST] = base_ptr[GoldTestEnum::TEST][j];
      hit_vector_end = hit_vector_end || is_bcf_vector_end_value<DataType>(val[GoldTestEnum::GOLD]);
      if(is_GT_field) {
        if(hit_vector_end) {
          //GenomicsDB must have processed all GT elements
          REQUIRE(static_cast<size_t>(test_j) == test_GT_vec_remapped_to_gold.size());
          //Test vcf line must have hit a vector end value also
          REQUIRE(is_bcf_vector_end_value<DataType>(val[GoldTestEnum::TEST]));
          break;
        }
        //Phase information starts from 2nd element in VCF GT
        if(j > 0) {
          REQUIRE(static_cast<size_t>(test_j) < test_GT_vec_remapped_to_gold.size());
          //If the schema, doesn't store phase with GT - then the query operation
          //always returns no phase
          if(bag.contains_phase)
            CHECK(bcf_gt_is_phased(val[GoldTestEnum::GOLD]) == test_GT_vec_remapped_to_gold[test_j]);
          else
            CHECK(test_GT_vec_remapped_to_gold[test_j] == 0);
          //Just check if the test vcf line matches the test data produced by GenomicsDBGVCFIterator
          //since the test data is already validated against gold
          CHECK(test_GT_vec[test_j] == bcf_gt_is_phased(val[GoldTestEnum::TEST]));
          ++test_j;
        }
        REQUIRE(static_cast<size_t>(test_j) < test_GT_vec_remapped_to_gold.size());
        if(g_skip_GT_matching && bcf_gt_allele(val[GoldTestEnum::GOLD]) != test_GT_vec_remapped_to_gold[test_j]) {
          WARN("GT mismatch for sample " << bcf_hdr_int2id(bag.hdr[GoldTestEnum::GOLD], BCF_DT_SAMPLE, i)
              << " at index " << j << " gold " << bcf_gt_allele(val[GoldTestEnum::GOLD]) << " test "
              << test_GT_vec_remapped_to_gold[test_j] << " - PL remapping is WIP");
        }
        else {
          CHECK(bcf_gt_allele(val[GoldTestEnum::GOLD]) == test_GT_vec_remapped_to_gold[test_j]);
          //Just check if the test vcf line matches the test data produced by GenomicsDBGVCFIterator
          //since the test data is already validated against gold
          CHECK(test_GT_vec[test_j] == bcf_gt_allele(val[GoldTestEnum::TEST]));
        }
        ++test_j;
        continue;
      }
      //Non GT fields
      //GenomicsDB field - all elements checked
      if(hit_vector_end && static_cast<unsigned>(j) == ptr_length_pair.second) {
        //Test vcf line must terminate the field with a vector end
        REQUIRE(is_bcf_vector_end_value<DataType>(val[GoldTestEnum::TEST]));
        break;
      }
      REQUIRE(static_cast<unsigned>(j) < ptr_length_pair.second);
      const auto test_val = (reinterpret_cast<const DataType*>(ptr_length_pair.first))[j];
      //Either both are missing/vector_end or are equal
      CHECK(((is_bcf_missing_value<DataType>(val[GoldTestEnum::GOLD]) && is_bcf_missing_value(test_val))
          || (is_bcf_vector_end_value<DataType>(val[GoldTestEnum::GOLD]) && is_bcf_vector_end_value<DataType>(test_val))
          || (val[GoldTestEnum::GOLD] == test_val)));
      CHECK(((is_bcf_missing_value<DataType>(val[GoldTestEnum::GOLD]) && is_bcf_missing_value(val[GoldTestEnum::TEST]))
          || (is_bcf_vector_end_value<DataType>(val[GoldTestEnum::GOLD]) && is_bcf_vector_end_value<DataType>(val[GoldTestEnum::TEST]))
          || (val[GoldTestEnum::GOLD] == val[GoldTestEnum::TEST])));
      if(hit_vector_end) {
        //Beyond vector end, either missing or vector_end values only
        REQUIRE((is_bcf_missing_value<DataType>(test_val) || is_bcf_vector_end_value<DataType>(test_val)));
        REQUIRE((is_bcf_missing_value<DataType>(val[GoldTestEnum::TEST])
              || is_bcf_vector_end_value<DataType>(val[GoldTestEnum::TEST])));
      }
    }
  }
}

TEST_CASE("columnar_gvcf_iterator_test", "[gvcf_iterator]") {
  if (g_query_json_file.empty() || g_golden_output_file.empty()) {
    // Sanity testing with json files from ctests input dir when --query-json-file is not specified
    g_loader_json_file = loader_json;
    g_query_json_file = ctests_input_dir + "query_all_attributes.json";
    g_golden_output_file = std::string(GENOMICSDB_TESTS_SRC_DIR) + "golden_outputs/t0_1_2_vcf_at_0";
  }
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
  //We'll be running two levels of checks:
  //(a) testing the output of GVCF iterator and remappers
  //(b) testing the output of VCFWriter (which internally invokes the iterator and remapper functions)

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
  bag.hdr[GoldTestEnum::GOLD] = bcf_hdr_read(fptr);
  REQUIRE(bag.hdr[GoldTestEnum::GOLD] != 0);
  if(!query_config.sites_only_query())
  {
    REQUIRE((uint64_t)bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]) == query_config.get_num_rows_to_query());
    //Map between samples in gold hdr[GoldTestEnum::GOLD] and test output
    bag.vcf_sample_idx_to_row_query_idx.resize(bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]));
    for(auto i=0;i<bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]);++i)
    {
      int64_t row_idx = 0;
      REQUIRE(vid_mapper.get_tiledb_row_idx(row_idx, bcf_hdr_int2id(bag.hdr[GoldTestEnum::GOLD], BCF_DT_SAMPLE, i)));
      bag.vcf_sample_idx_to_row_query_idx[i] = query_config.get_query_row_idx_for_array_row_idx(row_idx);
    }
  }
  bag.rec[GoldTestEnum::GOLD] = bcf_init();
  bag.rec[GoldTestEnum::TEST] = bcf_init();
  //Avoid repeated memory reallocations by reusing a buffer
  for(auto i=0u;i<2u;++i) {
    bag.vcf_buffer_length[i] = 4*bcf_hdr_nsamples(bag.hdr[GoldTestEnum::GOLD]); //since SB has 4 elements per sample - will get resized if necessary
    bag.vcf_buffer_ptr[i] = new int[bag.vcf_buffer_length[i]]; //use single buffer
  }
  //For test VCF
  VCFAdapter vcf_adapter(false);
  vcf_adapter.initialize(query_config);
  //Initialize operator so that you write to a string without any size limits
  BroadCombinedGVCFOperator broad_op(vcf_adapter, vid_mapper, query_config, false, true, true);
  //reference to string containing the VCF lines that broad_op will provide
  //initially this will contain the VCF header
  auto& test_vcf_line = broad_op.get_vcf_writer_to_string().get_string_buffer();
  bag.hdr[GoldTestEnum::TEST] = bcf_hdr_init("r");
  REQUIRE(
      bcf_hdr_deserialize(bag.hdr[GoldTestEnum::TEST], reinterpret_cast<const uint8_t*>(test_vcf_line.data()), 0u, test_vcf_line.size(), 0)
       == test_vcf_line.size());
  //Header is generated using htslib - no point in unit testing it
  test_vcf_line.clear();
  //Check SB and GQ fields since they have no allele dependence and are unmodified from input to output
  unsigned SB_query_idx = 0u;
  auto is_SB_queried = query_config.get_query_idx_for_name("SB", SB_query_idx);
  unsigned GQ_query_idx = 0u;
  auto is_GQ_queried = query_config.get_query_idx_for_name("GQ", GQ_query_idx);
  unsigned GT_query_idx = 0u;
  auto is_GT_queried = query_config.get_query_idx_for_name("GT", GT_query_idx);
  //VCF info for SB, GQ, GT
  auto do_check_SB = is_SB_queried && (bcf_hdr_id2int(bag.hdr[GoldTestEnum::GOLD], BCF_DT_ID, "SB") >= 0);
  auto do_check_GQ = is_GQ_queried && (bcf_hdr_id2int(bag.hdr[GoldTestEnum::GOLD], BCF_DT_ID, "GQ") >= 0);
  auto do_check_GT = is_GT_queried && (bcf_hdr_id2int(bag.hdr[GoldTestEnum::GOLD], BCF_DT_ID, "GT") >= 0);
  std::string contig_name;
  //REF-ALT comparison
  std::vector<STRING_VIEW> gold_alleles;
  std::vector<STRING_VIEW> test_alleles;
  std::vector<STRING_VIEW> test_vcf_line_alleles;
  std::string alleles_buffer;
  std::string tmp_buffer;
  std::vector<STRING_VIEW> tmp_alleles;
  //Using AllelesCombiner to keep track of mapping between gold alleles and test alleles
  //Will be used for GT remapping etc
  //Mark both gold and test as valid samples for remapping
  BagOfRemappers bag_of_remappers;
  bag_of_remappers.m_int_handler.set_valid_row_query_idx(GoldTestEnum::GOLD, true);
  bag_of_remappers.m_int_handler.set_valid_row_query_idx(GoldTestEnum::TEST, true);
  for (; !(bag.columnar_gvcf_iter->end()); ++(*(bag.columnar_gvcf_iter))) {
    //Gold
    REQUIRE(bcf_read(fptr, bag.hdr[GoldTestEnum::GOLD], bag.rec[GoldTestEnum::GOLD]) == 0);
    //Test vcf line
    test_vcf_line.clear(); //clear out last contents
    broad_op.operate_on_columnar_cell(**(bag.columnar_gvcf_iter));
    //Test line must be parsed fully by htslib without errors
    REQUIRE(bcf_deserialize(bag.rec[GoldTestEnum::TEST], reinterpret_cast<uint8_t*>(&(test_vcf_line.at(0u))),
        0u, test_vcf_line.length(), false, bag.hdr[GoldTestEnum::TEST]) == test_vcf_line.length());
    REQUIRE(bag.rec[GoldTestEnum::TEST]->errcode == 0);
    //Check contig and position
    auto column_interval = bag.columnar_gvcf_iter->get_current_variant_interval();
    int64_t contig_position = -1;
    REQUIRE(vid_mapper.get_contig_location(column_interval.first, contig_name, contig_position));
    REQUIRE(contig_name == std::string(bcf_seqname(bag.hdr[GoldTestEnum::GOLD], bag.rec[GoldTestEnum::GOLD])));
    REQUIRE(contig_position == bag.rec[GoldTestEnum::GOLD]->pos);
    //test vcf line
    REQUIRE(std::string(bcf_seqname(bag.hdr[GoldTestEnum::GOLD], bag.rec[GoldTestEnum::GOLD]))
         == std::string(bcf_seqname(bag.hdr[GoldTestEnum::TEST], bag.rec[GoldTestEnum::TEST])));
    REQUIRE(bag.rec[GoldTestEnum::GOLD]->pos == bag.rec[GoldTestEnum::TEST]->pos);
    //Variant length value
    //Get END values
    int num_END_value[2];
    num_END_value[GoldTestEnum::GOLD] = bcf_get_info_int32(bag.hdr[GoldTestEnum::GOLD], bag.rec[GoldTestEnum::GOLD], "END",
        &(bag.vcf_buffer_ptr[GoldTestEnum::GOLD]), &(bag.vcf_buffer_length[GoldTestEnum::GOLD]));
    num_END_value[GoldTestEnum::TEST] = bcf_get_info_int32(bag.hdr[GoldTestEnum::TEST], bag.rec[GoldTestEnum::TEST], "END",
        &(bag.vcf_buffer_ptr[GoldTestEnum::TEST]), &(bag.vcf_buffer_length[GoldTestEnum::TEST]));
    REQUIRE(num_END_value[GoldTestEnum::GOLD] == num_END_value[GoldTestEnum::TEST]);
    if(num_END_value[GoldTestEnum::GOLD] < 0) { //no END field, single position variant
      REQUIRE(1u == column_interval.second - column_interval.first + 1u);
    }
    else {
      REQUIRE(bag.rec[GoldTestEnum::GOLD]->rlen == column_interval.second - column_interval.first + 1u);
    }
    REQUIRE(bag.rec[GoldTestEnum::GOLD]->rlen == bag.rec[GoldTestEnum::TEST]->rlen);
    //Unpack
    bcf_unpack(bag.rec[GoldTestEnum::GOLD], BCF_UN_STR);
    bcf_unpack(bag.rec[GoldTestEnum::TEST], BCF_UN_STR);
    //REF-ALT comparison
    gold_alleles.clear();
    for(auto i=0u;i<bag.rec[GoldTestEnum::GOLD]->n_allele;++i)
      gold_alleles.emplace_back(STRING_VIEW(bag.rec[GoldTestEnum::GOLD]->d.allele[i], strlen(bag.rec[GoldTestEnum::GOLD]->d.allele[i])));
    alleles_buffer.clear();
    test_alleles.clear();
    bag.columnar_gvcf_iter->get_alleles_combiner().get_merged_VCF_spec_alleles_vec(alleles_buffer, test_alleles);
    if(test_alleles[0].length() == 0u) //leftover REF block remains
      test_alleles[0] = STRING_VIEW(bag.rec[GoldTestEnum::GOLD]->d.allele[0], strlen(bag.rec[GoldTestEnum::GOLD]->d.allele[0]));
    //Must have same alleles - possibly in a different order
    REQUIRE(gold_alleles[0u] == test_alleles[0u]);
    REQUIRE(std::unordered_set<STRING_VIEW>(gold_alleles.begin(), gold_alleles.end()) ==
        std::unordered_set<STRING_VIEW>(test_alleles.begin(), test_alleles.end()));
    //Find mapping from test alleles to gold alleles using AllelesCombiner
    bag_of_remappers.m_gold_test_alleles_combiner.remove_allele_info(GoldTestEnum::GOLD);
    bag_of_remappers.m_gold_test_alleles_combiner.remove_allele_info(GoldTestEnum::TEST);
    bag_of_remappers.m_gold_test_alleles_combiner.reset_before_adding_new_sample_info_at_current_position();
    bag_of_remappers.m_gold_test_alleles_combiner.insert_allele_info(GoldTestEnum::GOLD, gold_alleles[0],
        std::vector<STRING_VIEW>(gold_alleles.begin()+1u, gold_alleles.end()), false);
    bag_of_remappers.m_gold_test_alleles_combiner.insert_allele_info(GoldTestEnum::TEST, test_alleles[0],
        std::vector<STRING_VIEW>(test_alleles.begin()+1u, test_alleles.end()), false);
    bag_of_remappers.m_gold_test_alleles_combiner.finished_updating_allele_info_for_current_position();
    tmp_buffer.clear();
    tmp_alleles.clear();
    bag_of_remappers.m_gold_test_alleles_combiner.get_merged_VCF_spec_alleles_vec(tmp_buffer, tmp_alleles);
    //Merged allele order in bag_of_remappers.m_gold_test_alleles_combiner should be same as gold vcf
    REQUIRE(gold_alleles == tmp_alleles);
    //Check if mapping from test to gold is correct
    for(auto i=0u;i<test_alleles.size();++i) {
      auto merged_idx = bag_of_remappers.m_gold_test_alleles_combiner.get_merged_allele_idx(GoldTestEnum::TEST, i);
      REQUIRE(test_alleles[i] == gold_alleles[merged_idx]);
    }
    //Test vcf line
    test_vcf_line_alleles.clear();
    for(auto i=0u;i<bag.rec[GoldTestEnum::TEST]->n_allele;++i)
      test_vcf_line_alleles.emplace_back(STRING_VIEW(bag.rec[GoldTestEnum::TEST]->d.allele[i],
            strlen(bag.rec[GoldTestEnum::TEST]->d.allele[i])));
    //Test vcf line must have the same order of alleles in test_alleles (since they're produced from the same source)
    REQUIRE(test_alleles == test_vcf_line_alleles);
    if(query_config.sites_only_query())
      continue;
    //Check SB and GQ fields since they have no allele dependence and are unmodified from input to output
    if(do_check_SB)
      compare_field<int>(bag_of_remappers, bag, "SB", SB_query_idx);
    if(do_check_GQ)
      compare_field<int>(bag_of_remappers, bag, "GQ", GQ_query_idx);
    if(do_check_GT) {
      compare_field<int>(bag_of_remappers, bag, "GT", GT_query_idx);
    }
  }
  CHECK(bcf_read(fptr, bag.hdr[GoldTestEnum::GOLD], bag.rec[GoldTestEnum::GOLD]) == -1); //no more records
  delete[] bag.vcf_buffer_ptr[GoldTestEnum::GOLD];
  delete[] bag.vcf_buffer_ptr[GoldTestEnum::TEST];
  bcf_destroy(bag.rec[GoldTestEnum::GOLD]);
  bcf_destroy(bag.rec[GoldTestEnum::TEST]);
  bcf_hdr_destroy(bag.hdr[GoldTestEnum::GOLD]);
  bcf_hdr_destroy(bag.hdr[GoldTestEnum::TEST]);
  bcf_close(fptr);
  delete bag.columnar_gvcf_iter;
}
