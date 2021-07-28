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

#include "vcf2binary.h"
#include "tiledb_loader.h"
#include "cli.h"
#include <mpi.h>

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

int vcf_histogram_main(int argc, char** argv) {
  if (argc <= 1) {
    std::cerr << "Needs 1 arg <json_config_file>\n";
    //exit(-1);
    return -1;
  }
  //Converter object
  GenomicsDBImportConfig loader_config;
  loader_config.read_from_file(argv[1], 0);
  VCF2TileDBConverter converter(loader_config, 0);
  converter.create_and_print_histogram(argv[1]);
  return 0;
}
