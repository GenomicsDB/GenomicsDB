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

#include <algorithm>
#include "variant_query_config.h"
#include "query_variants.h"
#include "json_config.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBConfigException(#X);

using namespace std;

bool ColumnRangeCompare(const ColumnRange& x, const ColumnRange& y) {
  return (x.first < y.first);
}

void VariantQueryConfig::add_attribute_to_query(const string& name, unsigned schema_idx) {
  if (m_query_attribute_name_to_query_idx.find(name) == m_query_attribute_name_to_query_idx.end()) {
    auto idx = m_query_attributes_info_vec.size();
    m_query_attributes_info_vec.emplace_back(name, schema_idx);
    m_query_attribute_name_to_query_idx[name] = idx;
  }
}

void VariantQueryConfig::set_attributes_to_query(const vector<string>& attributeNames) {
  for (auto i=0u; i<attributeNames.size(); ++i)
    add_attribute_to_query(attributeNames[i], UNDEFINED_ATTRIBUTE_IDX_VALUE);
}

void VariantQueryConfig::clear_attributes_to_query() {
  m_query_attribute_name_to_query_idx.clear();
  m_query_attributes_info_vec.clear();
  m_query_idx_known_variant_field_enum_LUT.reset_luts();
}

void VariantQueryConfig::add_column_interval_to_query(const int64_t colBegin, const int64_t colEnd) {
  m_query_column_intervals.push_back(make_pair(colBegin, colEnd));
  std::sort(m_query_column_intervals.begin(), m_query_column_intervals.end(), ColumnRangeCompare);
}

void VariantQueryConfig::set_column_interval_to_query(const int64_t colBegin, const int64_t colEnd) {
  m_query_column_intervals.resize(1u);
  m_query_column_intervals[0] = make_pair(colBegin, colEnd);
}

void VariantQueryConfig::update_rows_to_query_to_all_rows() {
  if (m_query_all_rows) //already querying all rows
    return;
  assert(is_bookkeeping_done());
  m_query_all_rows = true;
}

void VariantQueryConfig::reorder_query_fields() {
  auto special_field_names = vector<string> { "END", "REF", "ALT" };
  m_first_normal_field_query_idx = 0u;
  for (auto i=0u; i<special_field_names.size(); ++i) {
    unsigned query_idx = 0u;
    auto& curr_field_name = special_field_names[i];
    if (get_query_idx_for_name(curr_field_name, query_idx)) {
      assert(query_idx >= m_first_normal_field_query_idx);
      if (query_idx > m_first_normal_field_query_idx) { // == implies already in right place
        auto& other_field_name = m_query_attributes_info_vec[m_first_normal_field_query_idx].m_name;
        //Before swap, update name mappings
        m_query_attribute_name_to_query_idx[curr_field_name] = m_first_normal_field_query_idx;
        m_query_attribute_name_to_query_idx[other_field_name] = query_idx;
        //Swap positions in schema idx and attribute names vector
        std::swap(m_query_attributes_info_vec[query_idx], m_query_attributes_info_vec[m_first_normal_field_query_idx]);
        //Now other_field_name can no longer be used
      }
      ++m_first_normal_field_query_idx;
    }
  }
}

void VariantQueryConfig::flatten_composite_fields(const VidMapper& vid_mapper) {
  auto has_composite_fields = false;
  for (auto i=0u; i<get_num_queried_attributes(); ++i) {
    auto field_info = vid_mapper.get_field_info(m_query_attributes_info_vec[i].m_name);
    if (field_info == 0)
      throw UnknownQueryAttributeException(std::string("Field ")
                                           +m_query_attributes_info_vec[i].m_name+" not found in vid mapping");
    auto num_elements_in_tuple = field_info->get_genomicsdb_type().get_num_elements_in_tuple();
    //Multi-element tuple - composite field
    if (num_elements_in_tuple > 1u) {
      has_composite_fields = true;
      for (auto j=0u; j<num_elements_in_tuple; ++j) {
        auto flattened_field_info = vid_mapper.get_flattened_field_info(field_info, j);
        add_attribute_to_query(flattened_field_info->m_name, UNDEFINED_ATTRIBUTE_IDX_VALUE);
      }
    }
  }
  if (has_composite_fields) {
    auto num_composite_fields_seen = 0u;
    //Remove composite fields, shift up m_query_attributes_info_vec and re-assign ids
    for (auto i=0u; i<get_num_queried_attributes(); ++i) {
      auto& field_name = m_query_attributes_info_vec[i].m_name;
      //Multi-element tuple - composite field
      if (vid_mapper.get_field_info(field_name)->get_genomicsdb_type().get_num_elements_in_tuple()
          > 1u)
        ++num_composite_fields_seen;
      else if (num_composite_fields_seen > 0u) { //shift-up elements in m_query_attributes_info_vec
        assert(m_query_attribute_name_to_query_idx[field_name] == i);
        m_query_attribute_name_to_query_idx[field_name] -= num_composite_fields_seen;
        m_query_attributes_info_vec[i-num_composite_fields_seen] =
          std::move(m_query_attributes_info_vec[i]);
      }
    }
    m_query_attributes_info_vec.resize(m_query_attributes_info_vec.size()-num_composite_fields_seen,
                                       VariantQueryFieldInfo("", UNDEFINED_ATTRIBUTE_IDX_VALUE));
  }
}

void VariantQueryConfig::read_from_file(const std::string& filename, const int rank) {
  GenomicsDBConfigBase::read_from_file(filename, rank);
  validate(rank);
}

void VariantQueryConfig::read_from_JSON_string(const std::string& str, const int rank) {
  GenomicsDBConfigBase::read_from_JSON_string(str, rank);
  validate(rank);
}

void VariantQueryConfig::read_from_PB(const genomicsdb_pb::ExportConfiguration* export_config, const int rank) {
  GenomicsDBConfigBase::read_from_PB(export_config, rank);
  validate(rank);
}

void VariantQueryConfig::read_from_PB_binary_string(const std::string& str, const int rank) {
  GenomicsDBConfigBase::read_from_PB_binary_string(str, rank);
  validate(rank);
}

void VariantQueryConfig::validate(const int rank) {
  //Workspace
  VERIFY_OR_THROW(m_workspaces.size() && "No workspace specified");
  VERIFY_OR_THROW((m_single_workspace_path || static_cast<size_t>(rank) < m_workspaces.size())
                  && ("Could not find workspace for rank "+std::to_string(rank)).c_str());
  auto& workspace = m_single_workspace_path ? m_workspaces[0] : m_workspaces[rank];
  VERIFY_OR_THROW(!workspace.empty() && "Empty workspace string");

  //Array
  VERIFY_OR_THROW(m_array_names.size() && "No array specified");
  VERIFY_OR_THROW((m_single_array_name || static_cast<size_t>(rank) < m_array_names.size())
                  && ("Could not find array for rank "+std::to_string(rank)).c_str());
  auto& array_name = m_single_array_name ? m_array_names[0] : m_array_names[rank];
  VERIFY_OR_THROW(!array_name.empty() && "Empty array name");

  //Initialize vid_mapper from file if necessary
  if (!m_vid_mapper.is_initialized()) {
    if (m_vid_mapping_file.size() > 0) {
      m_vid_mapper = std::move(FileBasedVidMapper(m_vid_mapping_file));
    }
    if (m_callset_mapping_file.size() > 0) {
      m_vid_mapper.parse_callsets_json(m_callset_mapping_file, true);
    }
  }
  
  //Query columns
  if (!m_scan_whole_array && m_column_ranges.size()) {
    VERIFY_OR_THROW((m_single_query_column_ranges_vector || static_cast<size_t>(rank) < m_column_ranges.size())
                    && "Rank >= query column ranges vector size");
    for (const auto& range : get_query_column_ranges(rank))
      add_column_interval_to_query(range.first, range.second);
  }

  //Query rows
  if (!m_scan_whole_array && m_row_ranges.size()) {
    VERIFY_OR_THROW((m_single_query_row_ranges_vector || static_cast<size_t>(rank) < m_row_ranges.size())
                    && "Rank >= query row ranges vector size");
    set_rows_to_query(get_query_row_ranges(rank));
  }
  //Attributes
  set_attributes_to_query(m_attributes);
}

void interval_expander::update_rows(const std::vector<std::pair<int64_t, int64_t>>& vec) {
  std::vector<interval> temp_intervals;
  std::swap(temp_intervals, intervals);

  for(auto& e : vec) {
    int64_t first = e.first, second = e.second;
    // make sure interval isn't backwards
    if(e.first <= e.second) {
      first = e.first;
      second = e.second;
    } else {
      first = e.second;
      second = e.first;
    }

    // all values allowed until clamp_low and clamp_high are called
    temp_intervals.push_back(interval(first, second));
  }
  // sort by left of intervals
  std::sort(temp_intervals.begin(), temp_intervals.end());

  bool current_empty = true;
  interval current;
  for(auto& ival : temp_intervals) { // merge overlapping intervals
    if(current_empty) {
      current = ival;
      current_empty = false;
    } else {
      if(current.intersects(ival)) {
        current = current + ival;
      } else {
        intervals.push_back(current);
        current = ival;
      }
    }
  }
  if(!current_empty) {
    intervals.push_back(current);
    current_empty = true;
  }

  for(auto& ival : intervals) {
    std::cerr << ival.first << " -- " << ival.second << std::endl;
  }

  index();
}

void interval_expander::set_rows(const std::vector<std::pair<int64_t, int64_t>>& vec) {
  clear();
  update_rows(vec);
}

void interval_expander::clamp_low(int64_t lo) {
  intervals.erase(std::remove_if(intervals.begin(),
                                 intervals.end(),
                                 [lo](interval& ival) -> bool { return ival.second < lo; }),
                  intervals.end());
  if(!intervals.size()) {
    return;
  }
  if(intervals[0].first < lo) { // only need to check first because intervals should be merged
    intervals[0].first = lo;
  }
}

void interval_expander::clamp_high(int64_t hi) {
  intervals.erase(std::remove_if(intervals.begin(),
                                 intervals.end(),
                                 [hi](interval& ival) -> bool { return ival.first > hi; }),
                  intervals.end());
  if(!intervals.size()) {
    return;
  }
  if(intervals.back().second > hi) { // only need to check last because intervals should be merged
    intervals.back().second = hi;
  }
}

void interval_expander::index() {
  int64_t total_size = 0;
  for(auto& ival : intervals) {
    ival.total_size_before = total_size;
    total_size += ival.size();
  }
}

int64_t interval_expander::get_array_row_from_query_row(int64_t row) const {
  auto it = std::upper_bound(intervals.begin(), intervals.end(), row, [] (int64_t l, const interval& r) { return l < r.size_inclusive(); });
  assert(it != intervals.end());
  int64_t offset_in_interval = row - it->total_size_before;
  return it->first + offset_in_interval;
}

// check if array row is queried
bool interval_expander::is_queried_row(int64_t row) const {
  auto it = std::lower_bound(intervals.begin(), intervals.end(), row, [row] (const interval& l, int64_t r) { return l.second < row; });
  return it != intervals.end() && it->first <= row;
}

int64_t interval_expander::operator[](int64_t idx) const {
  return get_array_row_from_query_row(idx);
}

int64_t interval_expander::get_query_row_from_array_row(int64_t row) const {
  auto it = std::lower_bound(intervals.begin(), intervals.end(), row, [row] (const interval& l, int64_t r) { return l.second < row; });
  assert(it != intervals.end() && it->first <= row);
  return it->total_size_before + (row - it->first);
}


void interval_expander::clear() {
  intervals.clear();
}

int64_t interval_expander::size() const {
  if(!intervals.size()) {
    return 0;
  }
  return intervals.back().size_inclusive();
}
