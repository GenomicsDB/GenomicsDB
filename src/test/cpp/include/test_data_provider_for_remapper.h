/**
 * The MIT License (MIT)
 * Copyright (c) 2021 Omics Data Automation, Inc.
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

#ifndef TEST_DATA_PROVIDER_FOR_REMAPPER_H
#define TEST_DATA_PROVIDER_FOR_REMAPPER_H

template<typename T>
class TestDataProviderForRemapper 
{
  public:
    TestDataProviderForRemapper(const size_t num_rows)
      : m_row_query_idx_to_data(num_rows) { }

    void set_data_for_row_query_idx(const uint64_t row_query_idx, const std::vector<T>& data) {
      assert(row_query_idx < m_row_query_idx_to_data.size());
      m_row_query_idx_to_data[row_query_idx] = data;
    }

    std::pair<const uint8_t*, size_t> get_raw_pointer_and_length_for_query_idx(const uint64_t row_query_idx,
        const int field_query_idx) const {
      assert(row_query_idx < m_row_query_idx_to_data.size());
      return std::pair<const uint8_t*, size_t>(
          reinterpret_cast<const uint8_t*>(&(m_row_query_idx_to_data[row_query_idx].front())),
          m_row_query_idx_to_data[row_query_idx].size()); 
    }

  private:
    std::vector<std::vector<T>> m_row_query_idx_to_data;
};

#endif
