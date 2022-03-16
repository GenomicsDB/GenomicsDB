/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
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

#ifndef VARIANT_QUERY_CONFIG_H
#define VARIANT_QUERY_CONFIG_H

#include "headers.h"
#include "lut.h"
#include "known_field_info.h"
#include "genomicsdb_config_base.h"
#include "genomicsdb_logger.h" // REMOVE

//Out of bounds query exception
class OutOfBoundsQueryException : public std::exception {
 public:
  OutOfBoundsQueryException(const std::string m="Out of bounds query exception") : msg_(m) { ; }
  ~OutOfBoundsQueryException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};

class UnknownQueryAttributeException : public std::exception {
 public:
  UnknownQueryAttributeException(const std::string m="Invalid queried attribute") : msg_(m) { ; }
  ~UnknownQueryAttributeException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};

class VariantQueryConfig : public GenomicsDBConfigBase {
 private:
  class VariantQueryFieldInfo {
   public:
    VariantQueryFieldInfo(const std::string& name, const int schema_idx)
      : m_name(name), m_schema_idx(schema_idx), m_vid_field_info(0) {
    }
    VariantQueryFieldInfo(const FieldInfo& vid_field_info, const int schema_idx)
      : m_schema_idx(schema_idx), m_vid_field_info(&vid_field_info) {
    }
    std::string m_name;
    int m_schema_idx;
    const FieldInfo* m_vid_field_info;
  };
 public:
  VariantQueryConfig()
    : GenomicsDBConfigBase() {
    clear();
    m_query_idx_known_variant_field_enum_LUT.reset_luts();
    m_done_bookkeeping = false;
    m_query_all_rows = true;
    m_num_rows_in_array = UNDEFINED_NUM_ROWS_VALUE;
    m_smallest_row_idx = 0;
    m_first_normal_field_query_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  }
  void clear() {
    m_query_attributes_info_vec.clear();
    m_query_attribute_name_to_query_idx.clear();
    m_query_rows.clear(); // REMOVE updated
    m_query_column_intervals.clear();
    // m_array_row_idx_to_query_row_idx.clear(); // REMOVE old
  }
  /**
   * Function that specifies which attributes to query from each cell
   */
  void set_attributes_to_query(const std::vector<std::string>& attributeNames);
  /**
   * Function used by query processor to add extra attributes to query
   */
  void add_attribute_to_query(const std::string& name, unsigned schema_idx);
  /*
   * Clear attributes to query
   */
  void clear_attributes_to_query();
  /**
   * Check whether attribute is defined
   */
  inline bool is_schema_idx_defined_for_query_idx(unsigned idx) const {
    assert(idx < m_query_attributes_info_vec.size());
    return (m_query_attributes_info_vec[idx].m_schema_idx != static_cast<int>(UNDEFINED_ATTRIBUTE_IDX_VALUE));
  }
  /**
   * Set TileDB array schema attribute idx (schemaIdx) for the queried attribute idx
   */
  void set_schema_idx_for_query_idx(unsigned idx, unsigned schemaIdx) {
    assert(idx < m_query_attributes_info_vec.size());
    m_query_attributes_info_vec[idx].m_schema_idx = schemaIdx;
  }
  /**
   * Get TileDB array schema attribute idx for the queried attribute idx
   */
  inline unsigned get_schema_idx_for_query_idx(unsigned idx) const {
    assert(idx < m_query_attributes_info_vec.size());
    return m_query_attributes_info_vec[idx].m_schema_idx;
  }
  /**
   * Get idx in the query for given attribute name
   */
  inline bool get_query_idx_for_name(const std::string& name, unsigned& idx) const {
    auto iter = m_query_attribute_name_to_query_idx.find(name);
    if (iter != m_query_attribute_name_to_query_idx.end()) {
      idx = (*iter).second;
      return true;
    } else
      return false;
  }
  /**
   * Get name for idx
   */
  inline std::string get_query_attribute_name(unsigned idx) const {
    assert(idx < m_query_attributes_info_vec.size());
    return m_query_attributes_info_vec[idx].m_name;
  }
  /**
   * Get number of attributes in query
   */
  inline unsigned get_num_queried_attributes() const {
    return m_query_attributes_info_vec.size();
  }
  std::vector<int> get_query_attributes_schema_idxs() const {
    auto schema_idx_vec = std::vector<int>(m_query_attributes_info_vec.size(), UNDEFINED_ATTRIBUTE_IDX_VALUE);
    for (auto i=0u; i<m_query_attributes_info_vec.size(); ++i)
      schema_idx_vec[i] = m_query_attributes_info_vec[i].m_schema_idx;
    return schema_idx_vec;
  }
  /*
   * Attributes info parameters
   */
  void set_query_attribute_info(const unsigned query_field_idx,
                                const FieldInfo& vid_field_info) {
    assert(query_field_idx < m_query_attributes_info_vec.size());
    m_query_attributes_info_vec[query_field_idx].m_vid_field_info = &vid_field_info;
  }
  /*
   * Flattens composite fields into multiple fields, rearranges
   * m_query_attributes_info_vec and sets m_query_attribute_name_to_query_idx
   */
  void flatten_composite_fields(const VidMapper& vid_mapper);
  const FieldLengthDescriptor& get_length_descriptor_for_query_attribute_idx(const unsigned query_idx) const {
    assert(query_idx < m_query_attributes_info_vec.size());
    assert(m_query_attributes_info_vec[query_idx].m_vid_field_info);
    return m_query_attributes_info_vec[query_idx].m_vid_field_info->m_length_descriptor;
  }
  int get_num_elements_for_query_attribute_idx(const unsigned query_idx) const {
    assert(query_idx < m_query_attributes_info_vec.size());
    return get_length_descriptor_for_query_attribute_idx(query_idx).get_num_elements();
  }
  int get_VCF_field_combine_operation_for_query_attribute_idx(const unsigned query_idx) const {
    assert(query_idx < m_query_attributes_info_vec.size());
    assert(m_query_attributes_info_vec[query_idx].m_vid_field_info);
    return m_query_attributes_info_vec[query_idx].m_vid_field_info->m_VCF_field_combine_operation;
  }
  std::type_index get_element_type(const unsigned query_idx) const {
    assert(query_idx < m_query_attributes_info_vec.size());
    assert(m_query_attributes_info_vec[query_idx].m_vid_field_info);
    return m_query_attributes_info_vec[query_idx].m_vid_field_info->get_genomicsdb_type().get_tuple_element_type_index(0u);
  }
  const FieldInfo* get_field_info_for_query_attribute_idx(const unsigned query_idx) const {
    assert(query_idx < m_query_attributes_info_vec.size());
    assert(m_query_attributes_info_vec[query_idx].m_vid_field_info);
    return (m_query_attributes_info_vec[query_idx].m_vid_field_info);
  }
  /*
   * Re-order query fields so that special fields like COORDS,END,NULL,OFFSET,ALT are first
   */
  void reorder_query_fields();
  unsigned get_first_normal_field_query_idx() const {
    return m_first_normal_field_query_idx;
  }
  //Query idx <--> known fields mapping
  void resize_LUT(unsigned num_known_fields) {
    m_query_idx_known_variant_field_enum_LUT.resize_luts_if_needed(get_num_queried_attributes(), num_known_fields);
  }
  //Map queryIdx <--> knownEnumIdx
  inline void add_query_idx_known_field_enum_mapping(unsigned queryIdx, unsigned knownEnumIdx) {
    assert(queryIdx < get_num_queried_attributes());
    assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
    m_query_idx_known_variant_field_enum_LUT.add_query_idx_known_field_enum_mapping(queryIdx, knownEnumIdx);
  }
  //Get query idx for given knownEnumIdx
  inline unsigned get_query_idx_for_known_field_enum(unsigned knownEnumIdx) const {
    assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
    return m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx);
  }
  //Get known field enum for given queryIdx
  inline unsigned get_known_field_enum_for_query_idx(unsigned queryIdx) const {
    assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(queryIdx));
    return m_query_idx_known_variant_field_enum_LUT.get_known_field_enum_for_query_idx(queryIdx);
  }
  //Check whether query contains given knownEnumIdx
  inline bool is_defined_query_idx_for_known_field_enum(unsigned knownEnumIdx) const {
    assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
    return m_query_idx_known_variant_field_enum_LUT.is_defined_value(
             m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx));
  }
  //Check whether query idx is known field
  inline bool is_defined_known_field_enum_for_query_idx(unsigned queryIdx) const {
    assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(queryIdx));
    return m_query_idx_known_variant_field_enum_LUT.is_defined_value(
             m_query_idx_known_variant_field_enum_LUT.get_known_field_enum_for_query_idx(queryIdx));
  }
  inline bool is_bookkeeping_done() const {
    return m_done_bookkeeping;
  }
  inline void set_done_bookkeeping(bool value) {
    m_done_bookkeeping = value;
  }
  /*
   * Function that specifies which rows to query
   */
  void set_rows_to_query(const std::vector<std::pair<int64_t, int64_t>>& rowIntervals) {
    logger.error("REMOVE VariantQueryConfig::set_rows_to_query m_query_all_rows {} -> false", m_query_all_rows);
    m_query_rows.set_rows(rowIntervals);
    m_query_all_rows = false;
  }
  
  /*
   * Clamps all row intervals to between lo and hi
   * Might result in no remaining row intervals
   */ 
  void clamp_query_rows(uint64_t lo, uint64_t hi) {
    m_query_rows.clamp_low(lo);
    m_query_rows.clamp_high(hi);
  }
  /*
   * Used by QueryProcessor to set number of rows if all rows need to be queried.
   */
  void set_num_rows_in_array(uint64_t num_rows, const uint64_t smallest_row_idx) {
    logger.error("REMOVE VariantQueryConfig::set_num_rows_in_array num rows {}, smallest {}", num_rows, smallest_row_idx);
    m_num_rows_in_array = num_rows;
    m_smallest_row_idx = smallest_row_idx;
    m_query_rows.clamp_low(m_smallest_row_idx);
    m_query_rows.clamp_high(m_num_rows_in_array - 1); // 0 based
  }
  /*
   * Function that specifies which rows to query and also updates
   * m_array_row_idx_to_query_row_idx // REMOVE comment
   * Pre-requisite: query bookkeeping should be done before calling this function // REMOVE no longer true
   */
  void update_rows_to_query(const std::vector<std::pair<int64_t, int64_t>>& rowIntervals) { set_rows_to_query(rowIntervals); } // NOTE: interval_expander::update_rows retains old rows
  /*
   * Function that specifies all rows should be queried.
   * m_array_row_idx_to_query_row_idx // REMOVE comment
   * Pre-requisite: query bookkeeping should be done before calling this function
   */
  void update_rows_to_query_to_all_rows();
  /**
   * Rows to query
   */
  inline bool query_all_rows() const {
    return m_query_all_rows;
  }
  /**
   * If all rows are queried, return m_num_rows_in_array (set by QueryProcessor)
   * Else return size of m_query_rows vector // REMOVE comment
   */
  inline uint64_t get_num_rows_to_query() const {
    logger.error("REMOVE VariantQueryConfig::get_num_rows_to_query m_query_all_rows {}", m_query_all_rows);
    /*Either query subset of rows (in m_query_rows) or set the variable m_num_rows_in_array correctly*/ // REMOVE comment
    assert(!m_query_all_rows || (m_num_rows_in_array != UNDEFINED_NUM_ROWS_VALUE));
    return m_query_all_rows ? m_num_rows_in_array : m_query_rows.size(); // REMOVE updated
  }
  inline int64_t get_smallest_row_idx_in_array() const {
    return m_smallest_row_idx;
  }
  inline uint64_t get_num_rows_in_array() const {
    assert(m_num_rows_in_array != UNDEFINED_NUM_ROWS_VALUE);
    return m_num_rows_in_array;
  }
  /**
   * If all rows are queried, return idx
   * Else return m_query_rows[idx] // REMOVE comment
   */
  inline int64_t get_array_row_idx_for_query_row_idx(uint64_t idx) const {
    assert(idx < get_num_rows_to_query());
    return m_query_all_rows ? (idx+m_smallest_row_idx) : m_query_rows[idx]; // REMOVE updated
  }
  /*
   * Index in m_query_rows for given array row idx // REMOVE comment
   */
  inline uint64_t get_query_row_idx_for_array_row_idx(int64_t row_idx) const {
    logger.error("REMOVE VariantQueryConfig::get_query_row_idx_for_array_row_idx {} m_query_all_rows is {}", row_idx, m_query_all_rows);
    if(m_query_all_rows) {
      auto retval = row_idx - m_smallest_row_idx;
      assert(retval >= 0 && retval < get_num_rows_to_query());
      return retval;
    }
    return m_query_rows.get_query_row_from_array_row(row_idx);
  }
  /*
   * Check if this row is being queried or no
   */
  inline bool is_queried_array_row_idx(int64_t row_idx) const {
    return m_query_all_rows ? true : m_query_rows.is_queried_row(row_idx);
  }
  /*
   * Function that specifies which column ranges to query
   * Note that the order in which these intervals are queried is NOT the same order in which
   * the intervals are added. Book-keeping sorts the intervals in ascending order so that
   * intervals close by have a chance of re-using cached tiles in memory
   */
  void add_column_interval_to_query(const int64_t colBegin, const int64_t colEnd);
  /*
   * Function that sets the interval as the only interval to be queried
   */
  void set_column_interval_to_query(const int64_t colBegin, const int64_t colEnd);
  /**
   * Returns number of ranges queried
   */
  inline unsigned get_num_column_intervals() const {
    return  m_query_column_intervals.size();
  }
  inline ColumnRange get_column_interval(unsigned idx) const {
    assert(idx < m_query_column_intervals.size());
    return m_query_column_intervals[idx];
  }
  inline uint64_t get_column_begin(unsigned idx) const {
    return get_column_interval(idx).first;
  }
  inline uint64_t get_column_end(unsigned idx) const {
    return get_column_interval(idx).second;
  }
  /*
   * Read configuration from JSON file
   */
  void read_from_file(const std::string& filename, const int rank=0);
  /*
   * Read configuration from JSON string
   */
  void read_from_JSON_string(const std::string& str, const int rank=0);
  /*
   * Read configuration from protobuf based export configuration
   */
  void read_from_PB(const genomicsdb_pb::ExportConfiguration* export_config, const int rank=0);
  /*
   * Read configuration from bytes that can be parsed by protobuf's export configuration
   */
  void read_from_PB_binary_string(const std::string& str, const int rank=0);
  /*
   * Validates and intializes variant query configuration. GenomicsDBConfigException is thrown on failed checks.
   */
  void validate(const int rank=0);
 private:
  /*
   * Function to invalid TileDB array row idx -> query row idx mapping
   * @param all_rows if true, invalidates all mappings, else invalidates mapping for rows in m_query_rows only // REMOVE comment
   */
  std::vector<VariantQueryFieldInfo> m_query_attributes_info_vec;
  //Map from query name to index in m_query_attributes_info_vec
  std::unordered_map<std::string, unsigned> m_query_attribute_name_to_query_idx;
  //Flag that tracks whether book-keeping is done
  bool m_done_bookkeeping;
  //Idx in m_query_attributes_info_vec with the first common attribute - see reorder_query_fields();
  unsigned m_first_normal_field_query_idx;
  //Mapping between queried idx and known fields enum
  QueryIdxToKnownVariantFieldsEnumLUT m_query_idx_known_variant_field_enum_LUT;
  /*Rows to query*/
  bool m_query_all_rows;
  /*vector of specific row idxs to query*/
  //std::vector<int64_t> m_query_rows; // REMOVE old
  /*vector mapping array row_idx to query row idx*/
  //std::vector<uint64_t> m_array_row_idx_to_query_row_idx; // REMOVE old
  /*Set by query processor*/
  uint64_t m_num_rows_in_array;
  int64_t m_smallest_row_idx;
  /*Column ranges to query*/
  std::vector<ColumnRange> m_query_column_intervals;

  class interval_expander {
    public:
      struct interval {
        // left of interval, right of interval, size of interval, size of all preceeding intervals
        int64_t first, second, total_size_before;
        interval(int64_t first = -1, int64_t second = -1, int64_t total_size_before = -1) : first(first), second(second), total_size_before(total_size_before) { }
        bool operator<(const interval& r) const {
          return first < r.first;
        }
        // number of rows in interval
        inline int64_t size() const { return second - first + 1; }
        inline int64_t span() const { return second - first; }
        // total number of rows in interval and all preceeding
        inline int64_t size_inclusive() const { return second - first + 1 + total_size_before; }
        bool intersects(const interval& r) const {
          return second >= r.first && first <= r.second;
        }
        // bounding interval
        interval operator+(const interval& r) const {
          int64_t nfirst, nsecond;
          nfirst = (first <= r.first ? first : r.first);
          nsecond = (second >= r.second ? second : r.second);

          return interval(nfirst, nsecond);
        }
      };

      void update_rows(const std::vector<std::pair<int64_t, int64_t>>& vec) {
        logger.error("REMOVE interval_expander::update_rows {} pairs", vec.size());
        std::vector<interval> temp_intervals;
        std::swap(temp_intervals, intervals);

        for(auto& e : vec) {
          int64_t first = e.first, second = e.second;
          // make sure interval isn't backwards
          if(e.first <= e.second) {
            first = e.first;
            second = e.second;
          }
          else {
            first = e.second;
            second = e.first;
          }

          // make sure interval is in bounds (on the left)
          if(first < 0) {
            first = 0;
          }
          if(second >= 0) {
            temp_intervals.push_back(interval(first, second));
          }
        }
        // sort by left of intervals
        std::sort(temp_intervals.begin(), temp_intervals.end());
        logger.error("REMOVE after sort");
        for(auto& ival : temp_intervals) {
          std::cerr << ival.first << " -- " << ival.second << std::endl;
        }


        bool current_empty = true;
        interval current;
        for(auto& ival : temp_intervals) { // merge overlapping intervals
          if(current_empty) {
            current = ival;
            current_empty = false;
            logger.error("\tREMOVE set current to {} {}", current.first, current.second);
          }
          else {
            if(current.intersects(ival)) {
              current = current + ival;
              logger.error("\tREMOVE {} {} intersects with {} {}", current.first, current.second, ival.first, ival.second);
            }
            else {
              logger.error("\tREMOVE {} {} does not intersect with {} {}, pushing first interval", current.first, current.second, ival.first, ival.second);
              intervals.push_back(current);
              current = ival;
            }
          }
        }
        if(!current_empty) {
          logger.error("REMOVE pushing current {} {}", current.first, current.second);
          intervals.push_back(current);
          current_empty = true;
        }

        logger.error("REMOVE end of update_rows {} intervals", intervals.size());
        logger.error("REMOVE contents of intervals");
        for(auto& ival : intervals) {
          std::cerr << ival.first << " -- " << ival.second << std::endl;
        }

        index();
      }

      void set_rows(const std::vector<std::pair<int64_t, int64_t>>& vec) {
        logger.error("REMOVE interval_expander::set_rows");
        clear();
        update_rows(vec);
      }

      // makes sure all array rows to be queried are greater than or equal to lo
      void clamp_low(int64_t lo) {
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

      // makes sure all array rows to be queried are greater than or equal to hi
      void clamp_high(int64_t hi) {
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

      void index() {
        logger.error("REMOVE interval_expander::index");

        int64_t total_size = 0;
        for(auto& ival : intervals) {
          ival.total_size_before = total_size;
          total_size += ival.size();
        }
      }

      // TODO buffering/binary search
      int64_t get_array_row_from_query_row(int64_t row) const {
        logger.error("REMOVE interval_expander::get_array_row_from_query_row {}", row);
        auto it = std::upper_bound(intervals.begin(), intervals.end(), row, [] (int64_t l, const interval& r) { logger.error("REMOVE ival {} {} size inc {}", r.first, r.second, r.size_inclusive()); return l < r.size_inclusive(); });
        assert(it != intervals.end());
        int64_t offset_in_interval = row - it->total_size_before;
        return it->first + offset_in_interval;
      };

      // check if array row is queried
      bool is_queried_row(int64_t row) const {
        logger.error("REMOVE interval_expander::is_queried_row");
        auto it = std::find_if(intervals.begin(), intervals.end(), [row] (const interval& ival) { return row >= ival.first && row <= ival.second; });
        return it != intervals.end();
      }

      inline int64_t operator[](int64_t idx) const {
        return get_array_row_from_query_row(idx);
      }

      // TODO use upper_bound to use binary search
      int64_t get_query_row_from_array_row(int64_t row) const {
        logger.error("REMOVE interval_expander::get_query_row_from_array_row {}", row);
        auto it = std::find_if(intervals.begin(), intervals.end(), [row] (const interval& ival) { logger.error("REMOVE ival {} {}", ival.first, ival.second); return row >= ival.first && row <= ival.second; });
        assert(it != intervals.end());
        return it->total_size_before + (row - it->first);
      }

      interval_expander(const std::vector<std::pair<int64_t, int64_t>>& vec = {}) {
        logger.error("REMOVE interval_expander ctor with {} pairs", vec.size());
        update_rows(vec);
      }

      void clear() {
        intervals.clear();
      }

      int64_t size() const {
        if(!intervals.size()) {
          return 0;
        }
        return intervals.back().size_inclusive();
      }

    private:
      std::vector<interval> intervals;
  };

  interval_expander m_query_rows; // REMOVE updated
};

#endif
