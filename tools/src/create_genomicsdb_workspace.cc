/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 * Copyright (c) 2023 dātma, inc™
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

#include "common.h"
#include "tiledb_utils.h"

#include <spdlog/sinks/stdout_color_sinks.h>

static Logger g_logger(Logger::get_logger("create_genomicsdb_workspace"));

void print_usage() {
  std::cerr << "Usage: create_genomicsdb_workspace [options] <workspace_directory>\n"
            << "where options include:\n"
            << "\t \e[1m--help\e[0m, \e[1m-h\e[0m Print a usage message summarizing options available and exit\n"
            << "\t \e[1m--version\e[0m Print version and exit\n"
            << "\t <workspace_directory> GenomicsDB workspace URI\n"
            << std::endl;
      }

int main(int argc, char** argv) {
  static struct option long_options[] = {
    {"help",                       0, 0, 'h'},
    {"version",                    0, 0, VERSION},
    {0,0,0,0}
  };

  int c;
  while ((c=getopt_long(argc, argv, "h", long_options, NULL)) >= 0) {
    switch (c) {
      case 'h':
        print_usage();
        return OK;
      case VERSION:
        std::cout << GENOMICSDB_VERSION << "\n";
        return OK;
      default:
        std::cerr << "Unknown command line argument\n";
        print_usage();
        return ERR;
    }
  }

  if (optind >= argc) {
    g_logger.error("Workspace directory required to be specified. See usage");
    print_usage();
    return ERR;
  }

  int rc = TileDBUtils::create_workspace(argv[optind], false);
  if (rc < 0) {
    g_logger.error("Could not create workspace at {}", argv[optind]);
    LOG_TILEDB_ERRMSG;
    return ERR;
  } else if (rc) {
    g_logger.error("Workspace {} already exists and is unchanged", argv[optind]);
    return OK;
  }

  g_logger.info("Successfully created workspace {}", argv[optind]);
  return OK;
}
