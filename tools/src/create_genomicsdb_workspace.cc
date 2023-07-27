/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#include <iostream>
#include "tiledb_utils.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Needs 1 argument <workspace_directory>\n";
    exit(-1);
  }
  int rc = TileDBUtils::create_workspace(argv[1], false);
  switch (rc) {
    case 0: // OK
      std::cout << "Workspace successfully created at " << argv[1] << std::endl;
      break;
    case -1: // NOT_DIR
      std::cerr << "Could not create a workspace as " << argv[1] << " exists and is not a directory" << std::endl;
      break;
    case -2: // NOT_CREATED - error condition
      if (strnlen(tiledb_errmsg, TILEDB_ERRMSG_MAX_LEN) > 0) {
         std::cerr << "TileDB error message:" << tiledb_errmsg << std::endl;
      }
      std::cerr << "Could not create a workspace at " << argv[1] << std::endl;
      break;
    case 1: //UNCHANGED
      std::cout << "Workspace already exists at " << argv[1] << " and is unchanged" << std::endl;
      break;
    default:
      if (strnlen(tiledb_errmsg, TILEDB_ERRMSG_MAX_LEN) > 0) {
         std::cerr << "TileDB error message:" << tiledb_errmsg << std::endl;
      }
      std::cerr << "Could not create a workspace at " << argv[1] << " return code=" << rc << std::endl;
  }
  return rc;
}
