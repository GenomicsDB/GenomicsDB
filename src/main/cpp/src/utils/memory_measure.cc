/**
 * The MIT License (MIT)
 * Copyright (c) 2024 dātma, inc™
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

#include "genomicsdb_logger.h"
#include "memory_measure.h"
#include <string>

void read_off_memory_status(statm_t& result, const size_t page_size) {
#ifdef __linux__
  const char* statm_path = "/proc/self/statm";

  FILE *f = fopen(statm_path,"r");
  if (!f) {
    perror(statm_path);
    abort();
  }
  if (7 != fscanf(f,"%lu %lu %lu %lu %lu %lu %lu",
                  &result.size,&result.resident,&result.share,&result.text,&result.lib,&result.data,&result.dt)) {
    perror(statm_path);
    abort();
  }
  result.size *= page_size;
  result.resident *= page_size;
  result.share *= page_size;
  result.text *= page_size;
  result.lib *= page_size;
  result.data *= page_size;
  result.dt *= page_size;
  fclose(f);
#endif
}

void print_memory_usage(const std::string& msg) {
#ifdef __linux__
  statm_t mem_result;
  read_off_memory_status(mem_result);
  logger.info("Mem usage {} rss={}M", msg, mem_result.resident/1000000);
#endif
}
