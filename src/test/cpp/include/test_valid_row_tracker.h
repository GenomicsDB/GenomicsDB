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

#ifndef TEST_VALID_ROW_TRACKER_H
#define TEST_VALID_ROW_TRACKER_H

class ValidRowTrackerForUnitTest
{
  public:
    ValidRowTrackerForUnitTest(const size_t num_rows)
    {
    }
    bool is_valid_row_query_idx(const uint64_t row_query_idx) const { return (m_valid_row_query_idx_set.count(row_query_idx) > 0u); }
    std::set<uint64_t>::const_iterator begin_valid_row_query_idx() const { return m_valid_row_query_idx_set.begin(); }
    std::set<uint64_t>::const_iterator end_valid_row_query_idx() const { return m_valid_row_query_idx_set.end(); }
    void set_valid_row_query_idx(const uint64_t row_query_idx, const bool v)
    {
      if(v)
        m_valid_row_query_idx_set.insert(row_query_idx);
      else
        m_valid_row_query_idx_set.erase(row_query_idx);
    }
  private:
    //not optimal, purely for unit test
    std::set<uint64_t> m_valid_row_query_idx_set;
};

template<typename ValidRowTrackerTy>
void insert_allele_info(ValidRowTrackerTy& tracker, AllelesCombiner<ValidRowTrackerTy>& combiner,
    const uint64_t row_query_idx, const std::string& REF, const std::string& ALT_delimited_str)
{
  tracker.set_valid_row_query_idx(row_query_idx, true); //set valid before modifying combiner
  combiner.insert_allele_info(row_query_idx, REF, ALT_delimited_str, false);
}

template<typename ValidRowTrackerTy>
void remove_allele_info(ValidRowTrackerTy& tracker, AllelesCombiner<ValidRowTrackerTy>& combiner,
    const uint64_t row_query_idx)
{
  combiner.remove_allele_info(row_query_idx);
  tracker.set_valid_row_query_idx(row_query_idx, false);
}

#endif
