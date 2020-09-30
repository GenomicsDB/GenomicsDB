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
