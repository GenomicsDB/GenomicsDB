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

#ifndef TEST_REMAPPED_DATA_RECEIVER_H
#define TEST_REMAPPED_DATA_RECEIVER_H

class RemappedDataReceiverForUnitTest 
{
  public:
    RemappedDataReceiverForUnitTest(const size_t num_rows)
      : m_row_query_idx_to_remapped_GT(num_rows),
      m_row_query_idx_to_is_GT_missing(num_rows, false) {}

    void clear_remapped_GT() {
      for(auto& vec : m_row_query_idx_to_remapped_GT)
        vec.clear();
      m_row_query_idx_to_is_GT_missing.assign(m_row_query_idx_to_is_GT_missing.size(), false);
    }
    bool write_GT_allele_index(const uint64_t row_query_idx, const int v) {
      assert(row_query_idx < m_row_query_idx_to_remapped_GT.size());
      m_row_query_idx_to_remapped_GT[row_query_idx].push_back(v);
      return true;
    }
    bool write_GT_phase(const uint64_t row_query_idx, const int v) {
      assert(row_query_idx < m_row_query_idx_to_remapped_GT.size());
      m_row_query_idx_to_remapped_GT[row_query_idx].push_back(v);
      return true;
    }
    bool write_GT_empty(const uint64_t row_query_idx) {
      assert(row_query_idx < m_row_query_idx_to_remapped_GT.size());
      m_row_query_idx_to_is_GT_missing[row_query_idx] = true;
      return true;
    }
    const std::vector<int>& get_remapped_GT(const uint64_t row_query_idx) const {
      assert(row_query_idx < m_row_query_idx_to_remapped_GT.size());
      return m_row_query_idx_to_remapped_GT[row_query_idx];
    }
    bool is_GT_missing(const uint64_t row_query_idx) const {
      assert(row_query_idx < m_row_query_idx_to_remapped_GT.size());
      return m_row_query_idx_to_is_GT_missing[row_query_idx];
    }
  private:
    std::vector<std::vector<int>> m_row_query_idx_to_remapped_GT;
    std::vector<bool> m_row_query_idx_to_is_GT_missing;
};

//Useful if you're interested in processing 1 sample at a time and don't care about
//row_query_idx
class SingleRowRemappedDataReceiverForUnitTest 
{
  public:
    SingleRowRemappedDataReceiverForUnitTest()
      :  m_is_GT_missing(false) {}

    void clear_remapped_GT() {
      m_remapped_GT.clear();
      m_is_GT_missing = false;
    }
    bool write_GT_allele_index(const uint64_t row_query_idx, const int v) {
      m_remapped_GT.push_back(v);
      return true;
    }
    bool write_GT_phase(const uint64_t row_query_idx, const int v) {
      m_remapped_GT.push_back(v);
      return true;
    }
    bool write_GT_empty(const uint64_t row_query_idx) {
      m_is_GT_missing = true;
      return true;
    }
    const std::vector<int>& get_remapped_GT(const uint64_t row_query_idx) const {
      return m_remapped_GT;
    }
    bool is_GT_missing(const uint64_t row_query_idx) const {
      return m_is_GT_missing;
    }
  private:
    std::vector<int> m_remapped_GT;
    bool m_is_GT_missing;
};

#endif
