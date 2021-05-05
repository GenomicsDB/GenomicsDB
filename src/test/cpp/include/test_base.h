/**
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
 */

#ifndef GENOMICSDB_TEST_BASE_H
#define GENOMICSDB_TEST_BASE_H

#include "tiledb_utils.h"

class TempDir {
 public:
  TempDir() {
    create_temp_directory();
  }

  ~TempDir() {
    TileDBUtils::delete_dir(tmp_dirname_);
  }

  const std::string& get() {
    return tmp_dirname_;
  }

  const std::string append(const std::string path) {
    return append_slash(tmp_dirname_) + path;
  }

 private:
  std::string tmp_dirname_;

  const std::string append_slash(const std::string path) {
    if (path[path.length()-1] == '/') {
      return path;
    } else {
      return path + "/";
    }
  }

  void create_temp_directory() {
    std::string dirname_pattern("GenomicsDBXXXXXX");
    const char *tmp_dir = getenv("TMPDIR");
    if (tmp_dir == NULL) {
      tmp_dir = P_tmpdir; // defined in stdio
    }
    assert(tmp_dir != NULL);
    tmp_dirname_ = mkdtemp(const_cast<char *>((append_slash(tmp_dir)+dirname_pattern).c_str()));
  }
};

#endif /* GENOMICSDB_TEST_BASE_H */
