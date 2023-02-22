/**
 * src/test/cpp/src/test_query_variants.cc
 *
 * The MIT License (MIT)
 * Copyright (c) 2023 Omics Data Automation, Inc.
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
 * Test query_variants.cc
 *
 */

#include <catch2/catch.hpp>
#include "test_base.h"

#include "query_variants.h"
#include "tiledb_utils.h"

TEST_CASE_METHOD(TempDir, "constructor", "[query_variants]") {
  CHECK_THROWS_AS(new VariantStorageManager(""), std::exception);
  auto ws = append("ws");
  CHECK(TileDBUtils::create_workspace(ws) == TILEDB_OK);
  VariantStorageManager variant_storage_manager(ws);
  VidMapper vid_mapper;
  CHECK_THROWS_AS(new VariantQueryProcessor(&variant_storage_manager, "non_existent_array", vid_mapper), std::exception);
}

static std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

static std::string workspace(ctests_input_dir+"ws");
static std::string array_name("t0_1_2");
static std::string query_json(ctests_input_dir+"query.json");
static std::string loader_json(ctests_input_dir+"loader.json");

TEST_CASE("test basic queries", "[query_variants]") {
  VariantStorageManager storage_manager(workspace);

  VariantQueryConfig query_config;
  GenomicsDBImportConfig loader_config;
  loader_config.read_from_file(loader_json);
  query_config.update_from_loader(loader_config);
  query_config.read_from_file(query_json);
  query_config.set_array_name(array_name);
  query_config.validate();

  VariantQueryProcessor query_processor(&storage_manager, array_name, query_config.get_vid_mapper());
  query_processor.obtain_TileDB_attribute_idxs(query_processor.get_array_schema(), query_config);

  // Add a non existent attribute to query_config and check if exception is thrown
  // by VariantQueryProcessor
  query_config.clear_attributes_to_query();
  query_config.set_attributes_to_query({"non_existent_attribute"});
  CHECK_THROWS_AS(query_processor.obtain_TileDB_attribute_idxs(query_processor.get_array_schema(), query_config), std::exception);
}
