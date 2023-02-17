/**
 * @file src/test/cpp/src/test_htslib_plugin.cc
 *
 * @section LICENSE
 *
 * The MIT License
 * 
 * @copyright Copyright (c) 2020 Omics Data Automation, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @section DESCRIPTION
 *
 * Test htslib plugin support
 */

#include <catch2/catch.hpp>

#include "hfile_genomicsdb.h"
#include "test_base.h"

#include <fcntl.h>

TEST_CASE_METHOD(TempDir, "test htslib plugin", "[test_htslib_plugin]") {
  // No exceptions should be thrown in genomicsdb_htslib_plugin_initialize
  genomicsdb_htslib_plugin_initialize();

  CHECK(genomicsdb_filesystem_init(NULL, 0) == NULL);
  CHECK(genomicsdb_filesystem_init("", 0) == NULL);
  CHECK(genomicsdb_filesystem_init(get_temp_dir().c_str(), 0) == NULL); // Directories not supported
  CHECK(genomicsdb_filesystem_init("dfdfd://ddd", 0) == NULL);

  std::string filename = get_temp_dir()+"/foo";
  CHECK(genomicsdb_filesystem_init(filename.c_str(), 0) == NULL);
  CHECK(genomicsdb_filesystem_init(filename.c_str(), O_RDONLY) == NULL);
  void* ctx = genomicsdb_filesystem_init(filename.c_str(), O_RDWR);
  CHECK(ctx != NULL);

  char buffer[16];

  CHECK(genomicsdb_filesystem_write(ctx, "/dfdfd/dfdfd/non-existent-file", "Hello", 6) == -1);
  CHECK(genomicsdb_filesystem_read(ctx, "non-existent-file", 0, buffer, 15) == -1);

  CHECK(genomicsdb_filesystem_write(ctx, filename.c_str(), "Hello", 6) == 6);
  CHECK(genomicsdb_filesize(ctx, filename.c_str()) == 6);
  CHECK(genomicsdb_filesystem_read(ctx, filename.c_str(), 0, buffer, 6) == 6);
  CHECK(strcmp(buffer, "Hello") == 0);
  CHECK(genomicsdb_filesystem_write(ctx, filename.c_str(), " World", 7) == 7);
  CHECK(genomicsdb_filesystem_read(ctx, filename.c_str(), 6, buffer, 6) == 6);
  CHECK(strcmp(buffer, " World") == 0);
  CHECK(genomicsdb_filesystem_close(ctx, filename.c_str()) == 0);
}
