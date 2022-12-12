/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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
*/

#include <getopt.h>
#include <iostream>

#include "common.h"
#include "tiledb_loader.h"

#include <spdlog/sinks/stdout_color_sinks.h>

static Logger g_logger(Logger::get_logger("consolidate_genomicsdb_array"));

enum GenomicsDBArgsEnum {
  VERSION=1000,
};

void print_usage() {
  std::cerr << "Usage: consolidate_genomicsdb_array [options]\n"
            << "where options include:\n"
            << "\t \e[1m--help\e[0m, \e[1m-h\e[0m Print a usage message summarizing options available and exit\n"
            << "\t \e[1m--workspace\e[0m=<GenomicsDB workspace URI>, \e[1m-w\e[0m <GenomicsDB workspace URI>\n"
            << "\t\t Specify path to GenomicsDB workspace\n"
            << "\t \e[1m--array-name\e[0m=<Array Name>, \e[1m-a\e[0m <Array Name>\n"
            << "\t\t Specify name of array that requires consolidation\n"
            << "\t \e[1m--buffer-size\e[0m=<Buffer Size>, \e[1m-z\e[0m <Buffer Size>\n"
            << "\t\t Optional, default is 10M. Specify a buffer size for consolidation\n"
            << "\t \e[1m--batch-size\e[0m=<Batch Size>, \e[1m-z\e[0m <Batch Size>\n"
            << "\t\t (Experimental) Optional, default is all fragments considered. Specify a batch size for consolidating with a set of fragments at a time\n"
            << "\t \e[1m--shared-posixfs-optimizations\e[0m, \e[1m-p\e[0m\n"
            << "\t\t Optional, default is false. If specified, the array folder is not locked for read/write and file handles are kept open until a final close for write\n"
            << "\t \e[1m--version\e[0m Print version and exit\n"
            << std::endl;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage();
    return ERR;
  }

  static struct option long_options[] = {
    {"workspace",                    1, 0, 'w'},
    {"array-name",                   1, 0, 'a'},
    {"buffer-size",                  1, 0, 'z'},
    {"batch-size",                   1, 0, 'b'},
    {"shared-posixfs-optimizations", 0, 0, 'p'},
    {"version",                      0, 0, VERSION},
    {"help",                         0, 0, 'h'},
    {0,0,0,0}
  };

  // Set global log level to info as default
#ifdef DEBUG
  spdlog::set_level(spdlog::level::info);
#else
  spdlog::set_level(spdlog::level::debug);
#endif

  std::string workspace;
  std::string array_name;
  size_t buffer_size = TILEDB_CONSOLIDATION_BUFFER_SIZE;
  int batch_size = -1;
  bool enable_shared_posixfs_optimizations = false;

  int c;
  while ((c=getopt_long(argc, argv, "w:a:z:b:ph", long_options, NULL)) >= 0) {
    switch (c) {
      case 'w':
        workspace = std::move(std::string(optarg));
        break;
      case 'a':
        array_name = std::move(std::string(optarg));
        break;
      case 'z':
        buffer_size = std::stoll(optarg);
        break;
      case 'b':
        batch_size = std::stoi(optarg);
      case 'p':
        enable_shared_posixfs_optimizations = true;
        break;
      case 'h':
        print_usage();
        return 0;
      case VERSION:
        std::cout << GENOMICSDB_VERSION << "\n";
        return OK;
      default:
        std::cerr << "Unknown command line argument\n";
        print_usage();
        return ERR;
    }
  }

  // Validate options
  if (workspace.empty()) {
    g_logger.error("Workspace(--workspace/-w) required to be specified");
    return ERR;
  } else if (array_name.empty()) {
    g_logger.error("Array Name(--array_name/-a) required to be specified");
    return ERR;
  }

  g_logger.info("Starting consolidation of {} in {}", array_name, workspace);
  VCF2TileDBLoader::consolidate_tiledb_array(workspace.c_str(), array_name.c_str(),
                                             buffer_size, batch_size, enable_shared_posixfs_optimizations);
  g_logger.info("Consolidation of {} in {} completed", array_name, workspace);
  return 0;
}
