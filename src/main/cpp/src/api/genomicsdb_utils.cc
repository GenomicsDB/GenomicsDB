/**
 * @file genomicsdb_utils.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2022 Omics Data Automation, Inc.
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
 * @section DESCRIPTION
 *
 * genomicsdb utilities
 *
 **/
#include "genomicsdb_utils.h"
#include "tiledb_utils.h"


namespace genomicsdb {

std::string version() {
  return GENOMICSDB_VERSION;
}

bool is_file(const std::string& filename) {
  return TileDBUtils::is_file(filename);
}

ssize_t file_size(const std::string& filename) {
  return TileDBUtils::file_size(filename);
}

int read_entire_file(const std::string& filename, void **buffer, size_t *length) {
  return TileDBUtils::read_entire_file(filename, buffer, length);
}

std::vector<std::string> get_array_names(const std::string& workspace) {
  return TileDBUtils::get_array_names(workspace);
}

int cache_fragment_metadata(const std::string& workspace, const std::string& array_name) {
  return TileDBUtils::cache_fragment_metadata(workspace, array_name);
}

}
