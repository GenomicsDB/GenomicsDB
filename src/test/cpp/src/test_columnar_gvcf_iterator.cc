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
  
  VariantCounter counter;
  query_processor.iterate_over_gvcf_entries(query_processor.get_array_descriptor(), query_config, counter, true);
  CHECK(counter.get_value() == 4);

  VariantCounter *variant_counter = new VariantCounter();
  query_processor.scan_and_operate(query_processor.get_array_descriptor(), query_config, *variant_counter, 0, true);
  CHECK(variant_counter->get_value() == 4);
  delete variant_counter;
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
  GenomicsDBGVCFIterator* columnar_gvcf_iter = storage_manager.begin_gvcf_iterator(
      query_processor.get_array_descriptor(), query_config,
      true);
  REQUIRE(columnar_gvcf_iter != 0);
  //Golden VCF
  auto fptr = bcf_open(g_golden_output_file.c_str(), "rb");
  REQUIRE(fptr != 0);
  auto hdr = bcf_hdr_read(fptr);
  REQUIRE(hdr != 0);
  if(!query_config.sites_only_query())
    REQUIRE(bcf_hdr_nsamples(hdr) == query_config.get_num_rows_to_query());
  auto rec = bcf_init();
  //Check SB and GQ fields since they have no allele dependence and are unmodified from input to output
  unsigned SB_query_idx = 0u;
  auto is_SB_queried = query_config.get_query_idx_for_name("SB", SB_query_idx);
  unsigned GQ_query_idx = 0u;
  auto is_GQ_queried = query_config.get_query_idx_for_name("GQ", GQ_query_idx);
  //VCF info for SB and GQ
  auto SB_hdr_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, "SB");
  auto GQ_hdr_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, "GQ");
  //Avoid some memory reallocations
  auto vcf_buffer_length = 4*bcf_hdr_nsamples(hdr); //since SB has 4 elements per sample
  auto vcf_buffer_ptr = new int[vcf_buffer_length]; //use single buffer
  std::string contig_name;
  //REF-ALT comparison
  std::unordered_set<STRING_VIEW> gold_alleles;
  std::unordered_set<STRING_VIEW> test_alleles;
  std::string alleles_buffer;
  std::vector<STRING_VIEW> alleles_vec;
  for (; !(columnar_gvcf_iter->end()); ++(*columnar_gvcf_iter)) {
    auto column_interval = columnar_gvcf_iter->get_current_variant_interval();
    REQUIRE(bcf_read(fptr, hdr, rec) == 0);
    int64_t contig_position = -1;
    REQUIRE(vid_mapper.get_contig_location(column_interval.first, contig_name, contig_position));
    CHECK(contig_name == std::string(bcf_seqname(hdr, rec)));
    CHECK(contig_position == rec->pos);
    //std::cerr << rec->pos <<"\n";
    bcf_unpack(rec, BCF_UN_STR);
    //REF-ALT comparison
    gold_alleles.clear();
    for(auto i=0u;i<rec->n_allele;++i)
      gold_alleles.insert(STRING_VIEW(rec->d.allele[i], strlen(rec->d.allele[i])));
    alleles_buffer.clear();
    alleles_vec.clear();
    test_alleles.clear();
    columnar_gvcf_iter->get_alleles_combiner().get_merged_VCF_spec_alleles_vec(alleles_buffer, alleles_vec);
    if(alleles_vec[0].length() == 0u) //leftover REF block remains
      alleles_vec[0] = STRING_VIEW(rec->d.allele[0], strlen(rec->d.allele[0]));
    for(auto x : alleles_vec)
      test_alleles.insert(x);
    CHECK(gold_alleles == test_alleles);
    //Check SB and GQ fields since they have no allele dependence and are unmodified from input to output
    for(auto field_name : { "SB", "GQ" }) {
      if((strcmp(field_name, "SB") == 0 && is_SB_queried && SB_hdr_idx >= 0)
          || (strcmp(field_name, "GQ") == 0 && is_GQ_queried && GQ_hdr_idx >= 0)) {
        auto query_field_idx = (strcmp(field_name, "SB") == 0) ? SB_query_idx : GQ_query_idx;
        auto num_output = bcf_get_format_int32(hdr, rec, field_name, &vcf_buffer_ptr, &vcf_buffer_length);
        REQUIRE(((num_output == -3) || (num_output >= 0)));
        for(auto i=0;i<bcf_hdr_nsamples(hdr);++i) {
          int64_t row_idx = 0;
          REQUIRE(vid_mapper.get_tiledb_row_idx(row_idx, bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i)));
          auto row_query_idx_idx = query_config.get_query_row_idx_for_array_row_idx(row_idx);
          if(num_output == -3)
            CHECK((!columnar_gvcf_iter->is_valid(row_query_idx_idx)
                || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(query_field_idx, row_query_idx_idx)));
          else {
            auto num_per_sample = num_output/bcf_hdr_nsamples(hdr);
            auto one_valid = false;
            for(auto j=0;j<num_per_sample;++j)
              one_valid = one_valid || is_bcf_valid_value<int>(vcf_buffer_ptr[i*num_per_sample+j]);
            if(one_valid) {
              CHECK(columnar_gvcf_iter->is_field_valid_for_row_query_idx(query_field_idx, row_query_idx_idx));
              auto ptr_length_pair = columnar_gvcf_iter->get_raw_pointer_and_length_for_query_idx(row_query_idx_idx,
                  query_field_idx);
              for(auto j=0;j<num_per_sample;++j) {
                auto val = vcf_buffer_ptr[i*num_per_sample+j];
                if(is_bcf_valid_value<int>(val)) {
                  REQUIRE(static_cast<unsigned>(j) < ptr_length_pair.second);
                  CHECK(reinterpret_cast<const int*>(ptr_length_pair.first)[j] == val);
                }
              }
            }
            else
              CHECK((!columnar_gvcf_iter->is_valid(row_query_idx_idx)
                    || !columnar_gvcf_iter->is_field_valid_for_row_query_idx(query_field_idx, row_query_idx_idx)));
          }
        }
      }
    }
  }
  CHECK(bcf_read(fptr, hdr, rec) == -1); //no more records
  delete[] vcf_buffer_ptr;
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  bcf_close(fptr);
  delete columnar_gvcf_iter;
}
