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

#ifndef COMMON_HEADERS_H
#define COMMON_HEADERS_H

#include <stdio.h>
#include <typeinfo>
#include <stdlib.h>
#include <cstdint>
#include <string.h>
#include <iostream>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>
#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <set>
#include <map>
#include <climits>
#include <sstream>
#include <iomanip>
#include <exception>
#include <fstream>
#include <functional>
#ifndef DISABLE_OPENMP
#include <omp.h>
#endif
#include <typeindex>

#ifdef SAFESTRINGLIB_FOUND
#include "safe_mem_lib.h"
#else
#define memcpy_s(dst, dst_size, src, src_size) memcpy((dst), (src), (src_size))
#endif

#ifdef STRING_VIEW_FOUND
#include <string_view>
typedef std::string_view STRING_VIEW;
#elif defined(STRING_VIEW_EXPERIMENTAL_FOUND)
#include <experimental/string_view>
typedef std::experimental::string_view STRING_VIEW;
#else
//TODO: define a simple class to mimic string_view behavior
//to use with compilers which don't support string_view (for example, gcc 4.8)
#error "TBD: Use a compiler with string_view support for now"
#endif

typedef std::pair<int64_t, int64_t> ColumnRange;
typedef std::pair<int64_t, int64_t> TileDBRowRange;
bool ColumnRangeCompare(const ColumnRange& x, const ColumnRange& y);

class CircularBufferController {
 public:
  CircularBufferController(unsigned num_entries) {
    m_num_entries = num_entries;
    m_num_entries_with_valid_data = 0u;
    m_num_reserved_entries = 0u;
    m_curr_write_idx = 0u;
    m_curr_read_idx = 0u;
  }
  //Advance write idx - increase #valid entries, decrease #reserved entries if specified
  inline void advance_write_idx(bool unreserve=false) {
    m_curr_write_idx = (m_curr_write_idx+1u)%m_num_entries;
    assert(m_num_entries_with_valid_data < m_num_entries);
    ++m_num_entries_with_valid_data;
    if (unreserve) {
      assert(m_num_reserved_entries > 0u);
      --m_num_reserved_entries;
    }
  }
  //Reserves an entry without marking it as valid
  void reserve_entry()   {
    ++m_num_reserved_entries;
  }
  inline void advance_read_idx() {
    m_curr_read_idx = (m_curr_read_idx+1u)%m_num_entries;
    assert(m_num_entries_with_valid_data > 0u);
    --m_num_entries_with_valid_data;
  }
  inline unsigned get_num_entries_with_valid_data() const {
    return m_num_entries_with_valid_data;
  }
  inline unsigned get_num_empty_entries() const {
    return m_num_entries - m_num_entries_with_valid_data - m_num_reserved_entries;
  }
  //Get idx
  inline unsigned get_write_idx() const {
    return m_curr_write_idx;
  }
  inline unsigned get_read_idx() const {
    return m_curr_read_idx;
  }
 protected:
  //Points to latest entry with valid data
  unsigned m_curr_write_idx;
  //Points to entry being read currently
  unsigned m_curr_read_idx;
  unsigned m_num_entries;
  unsigned m_num_entries_with_valid_data;
  unsigned m_num_reserved_entries;
};

class RWBuffer {
 public:
  RWBuffer(size_t buffer_capacity=1048576u) {
    m_buffer.resize(buffer_capacity);
    m_next_read_idx = 0u;
    m_num_valid_bytes = 0u;
  }
  size_t get_num_remaining_bytes() const {
    return m_num_valid_bytes-m_next_read_idx;
  }
  uint8_t* get_pointer_at_read_position() {
    return &(m_buffer[m_next_read_idx]);
  }
  const uint8_t* get_pointer_at_read_position() const {
    return &(m_buffer[m_next_read_idx]);
  }
  std::vector<uint8_t> m_buffer;
  size_t m_next_read_idx;
  size_t m_num_valid_bytes;
};

#endif
