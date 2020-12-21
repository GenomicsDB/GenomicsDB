/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
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

#ifndef GENOMICSDB_ITERATORS_H
#define GENOMICSDB_ITERATORS_H

#include "headers.h"
#include "variant_array_schema.h"
#include "genomicsdb_columnar_field.h"
#include "vid_mapper.h"
#include "tiledb.h"
#include "timer.h"
#include "variant_query_config.h"

//Exceptions thrown
class GenomicsDBIteratorException : public std::exception {
 public:
  GenomicsDBIteratorException(const std::string m="") : msg_("GenomicsDBIteratorException exception : "+m) { ; }
  ~GenomicsDBIteratorException() { ; }
  // ACCESSORS
  /** Returns the exception message. */
  const char* what() const noexcept {
    return msg_.c_str();
  }
 private:
  std::string msg_;
};

//Whenever a live cell must be pointed to in a buffer, use this structure
//Useful for the left sweep, tracking cells while producing VCF records etc
class GenomicsDBLiveCellMarker {
 public:
  GenomicsDBLiveCellMarker(const unsigned num_fields, const size_t num_queried_rows) {
    m_store_gvcf_specific_info = false;
    m_valid.resize(num_queried_rows);
    m_initialized.resize(num_queried_rows);
    m_row.resize(num_queried_rows);
    m_begin.resize(num_queried_rows);
    m_end.resize(num_queried_rows);
    for (auto i=0ull; i<num_fields; ++i) {
      m_buffer_ptr_vec.emplace_back(num_queried_rows, static_cast<GenomicsDBBuffer*>(0));
      m_indexes.emplace_back(num_queried_rows);
    }
    reset();
  }
  inline void reset() {
    //Set columns to -1
    m_begin.assign(m_begin.size(), -1ll);
    m_initialized.assign(m_initialized.size(), false);
    m_valid.assign(m_valid.size(), false);
    if(m_store_gvcf_specific_info) {
      m_contains_deletion.assign(m_contains_deletion.size(), false);
      m_contains_MNV.assign(m_contains_MNV.size(), false);
    }
  }
  inline void prepare_for_gvcf_iterator() {
    m_store_gvcf_specific_info = true;
    const auto num_queried_rows = m_valid.size();
    m_contains_deletion.resize(num_queried_rows);
    m_contains_MNV.resize(num_queried_rows);
  }
  inline void set_row_idx(const size_t idx, const int64_t row_idx) {
    assert(idx < m_row.size());
    m_row[idx] = row_idx;
  }
  inline int64_t get_row_idx(const size_t idx) const {
    assert(idx < m_row.size());
    return m_row[idx];
  }
  inline void set_column_interval(const size_t idx, const int64_t b, const int64_t e) {
    assert(idx < m_begin.size() && idx < m_end.size());
    m_begin[idx] = b;
    m_end[idx] = e;
  }
  inline void set_initialized(const size_t idx, const bool val) {
    assert(idx < m_initialized.size());
    m_initialized[idx] = val;
  }
  inline bool is_initialized(const size_t idx) const {
    assert(idx < m_initialized.size());
    return m_initialized[idx];
  }
  inline void set_valid(const size_t idx, const bool val) {
    assert(idx < m_valid.size());
    m_valid[idx] = val;
  }
  inline bool is_valid(const size_t idx) const {
    assert(idx < m_valid.size());
    return m_valid[idx];
  }
  inline void set_field_marker(const unsigned field_idx, const size_t idx,
                               GenomicsDBBuffer* buffer_ptr, const size_t index) {
    assert(field_idx < m_buffer_ptr_vec.size() && field_idx < m_indexes.size());
    assert(idx < m_buffer_ptr_vec[field_idx].size() && idx < m_indexes[field_idx].size());
    m_buffer_ptr_vec[field_idx][idx] = buffer_ptr;
    m_indexes[field_idx][idx] = index;
  }
  inline GenomicsDBBuffer* get_buffer(const unsigned field_idx, const size_t idx) {
    assert(m_initialized[idx] && m_valid[idx]);
    assert(field_idx < m_buffer_ptr_vec.size() && field_idx < m_indexes.size());
    assert(idx < m_buffer_ptr_vec[field_idx].size() && idx < m_indexes[field_idx].size());
    return m_buffer_ptr_vec[field_idx][idx];
  }
  inline const GenomicsDBBuffer* get_buffer(const unsigned field_idx, const size_t idx) const {
    assert(m_initialized[idx] && m_valid[idx]);
    assert(field_idx < m_buffer_ptr_vec.size() && field_idx < m_indexes.size());
    assert(idx < m_buffer_ptr_vec[field_idx].size() && idx < m_indexes[field_idx].size());
    return m_buffer_ptr_vec[field_idx][idx];
  }
  inline size_t get_index(const unsigned field_idx, const size_t idx) const {
    assert(m_initialized[idx] && m_valid[idx]);
    assert(field_idx < m_buffer_ptr_vec.size() && field_idx < m_indexes.size());
    assert(idx < m_buffer_ptr_vec[field_idx].size() && idx < m_indexes[field_idx].size());
    return m_indexes[field_idx][idx];
  }
  inline int64_t get_begin(const size_t idx) const {
    assert(idx < m_begin.size());
    return m_begin[idx];
  }
  inline int64_t get_end(const size_t idx) const {
    assert(idx < m_end.size());
    return m_end[idx];
  }
  inline void set_contains_deletion(const size_t idx, const bool v) {
    assert(idx < m_contains_deletion.size());
    m_contains_deletion[idx] = v;
  }
  inline void set_contains_MNV(const size_t idx, const bool v) {
    assert(idx < m_contains_MNV.size());
    m_contains_MNV[idx] = v;
  }
  inline bool contains_deletion_or_MNV(const size_t idx) const {
    assert(idx < m_contains_deletion.size());
    assert(idx < m_contains_MNV.size());
    return m_contains_deletion[idx] || m_contains_MNV[idx];
  }
  std::vector<bool>::const_iterator valid_begin() const { return m_valid.begin(); }
  std::vector<bool>::const_iterator valid_end() const { return m_valid.end(); }
 private:
  //Outermost vector - 1 per cell
  //Keeping the theme of columnar structures
  std::vector<bool> m_initialized;
  std::vector<bool> m_valid;
  std::vector<int64_t> m_row;
  std::vector<int64_t> m_begin;
  std::vector<int64_t> m_end;
  //Used in the GVCF iterator
  bool m_store_gvcf_specific_info;
  std::vector<bool> m_contains_deletion;
  std::vector<bool> m_contains_MNV;
  //Outer vector - per field
  //Inner vector - per sample/CallSet
  std::vector<std::vector<GenomicsDBBuffer*> > m_buffer_ptr_vec;
  std::vector<std::vector<size_t> >m_indexes;
};

//m_PQ_live_cell_markers is a min-heap ordered in column major order
//Elements of the heap are pairs<column, row_query_idx>
//Earlier only the row_query_idx was stored in the heap and the column was obtained by
//accessing the live cell markers. Storing the column in the heap helps improve
//locality of accesses in the heap.
typedef std::pair<int64_t, size_t> ColumnRowQueryIdxPair;

class GenomicsDBLiveCellMarkerColumnMajorComparator {
  public:
    bool operator()(const ColumnRowQueryIdxPair& a, const ColumnRowQueryIdxPair& b) const {
      return !((a.first < b.first) || (a.first == b.first && a.second < b.second));
    }
};

enum GenomicsDBIteratorStatsEnum {
  TOTAL_CELLS_TRAVERSED=0,
  NUM_CELLS_TRAVERSED_IN_FIND_INTERSECTING_INTERVALS_MODE,
  NUM_USELESS_CELLS_IN_FIND_INTERSECTING_INTERVALS_MODE,
  NUM_USELESS_CELLS_IN_SIMPLE_TRAVERSAL_MODE,
  NUM_STATS
};

class GenomicsDBColumnarCell;
class VariantQueryConfig;
/*
 * Iterates over TileDB cells one at a time
 */
class SingleCellTileDBIterator {
 public:
  //Initializes tiledb_array object internally
  SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx,
                           const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
                           const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size);
  //Uses tiledb_array object provided by caller (if non-NULL)
  SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx,
                           const TileDB_Array* tiledb_array,
                           const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
                           const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size);
  ~SingleCellTileDBIterator();
  //Delete copy and move constructors
  SingleCellTileDBIterator(const SingleCellTileDBIterator& other) = delete;
  SingleCellTileDBIterator& operator=(const SingleCellTileDBIterator& other) = delete;
  SingleCellTileDBIterator(SingleCellTileDBIterator&& other) = delete;
  //Iterator functionality
  inline const GenomicsDBColumnarCell& operator*() const {
    return *m_cell;
  }
  /*
   * Should be called by consumers of this iterator - not internally
   * See advance_*() functions below for internal use
   */
  const SingleCellTileDBIterator& operator++();
  //Get pointer, length, size for current cell for field corresponding to query_idx
  inline const uint8_t* get_raw_field_pointer_for_query_idx(const int query_idx) const {
    assert(static_cast<size_t>(query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto x = get_buffer_and_index(query_idx);
    return genomicsdb_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
             x.first, x.second);
  }
  inline size_t get_field_length(const int query_idx) const {
    assert(static_cast<size_t>(query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto x = get_buffer_and_index(query_idx);
    return genomicsdb_columnar_field.get_length_of_data_in_buffer_at_index(
             x.first, x.second);
  }
  inline size_t get_field_size_in_bytes(const int query_idx) const {
    assert(static_cast<size_t>(query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto x = get_buffer_and_index(query_idx);
    return genomicsdb_columnar_field.get_size_of_data_in_buffer_at_index(
             x.first, x.second);
  }
  inline bool is_valid(const int query_idx) const {
    assert(static_cast<size_t>(query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto x = get_buffer_and_index(query_idx);
    return genomicsdb_columnar_field.is_valid_data_in_buffer_at_index(
             x.first, x.second);
  }
  inline const int64_t* get_coordinates() const {
    auto coords_query_idx = m_fields.size()-1u;
    return reinterpret_cast<const int64_t*>(get_raw_field_pointer_for_query_idx(coords_query_idx));
  }
  void print(const int query_idx, std::ostream& fptr=std::cout) const;
  void print_ALT(const int query_idx, std::ostream& fptr=std::cout) const;
  void print_csv(const int query_idx, std::ostream& fptr=std::cout) const;
  inline bool end() const {
    return m_done_reading_from_TileDB && m_PQ_live_cell_markers.empty();
  }
  inline bool at_new_query_column_interval() const {
    return m_at_new_query_column_interval;
  }
  inline uint64_t get_current_query_column_interval_idx() const {
    return m_query_column_interval_idx;
  }
 protected:
  /*
   * These functions must be called internally only. Consumers of iterator objects
   * should use operator++ only.
   */
  /*
   * Function that starts a new query interval
   * (a) For scan, starts the scan in simple traversal mode
   * (b) With intervals
   *     * Determines all intersecting intervals in find intersecting intervals mode
   *     * Moves into simple traversal mode
   *     * Finds first 'useful' cell in the simple traversal mode
   * The first call to this function must pass the 3 arguments correctly
   * Subsequent calls can skip the arguments
   */
  void begin_new_query_column_interval(TileDB_CTX* tiledb_ctx=0, const char* array_path=0,
                                       std::vector<const char*>* attribute_names=0);
  //Helper functions - sets flags and preps some fields
  void move_into_simple_traversal_mode();
  void reset_for_next_query_interval();
  void handle_current_cell_in_find_intersecting_intervals_mode();
  /*
   * Does one read for the attributes in m_query_attribute_idx_vec from TileDB
   * skip_cells - for queried attributes, skip #cells specified in
   *              m_query_attribute_idx_num_cells_to_increment_vec
   */
  void read_from_TileDB(const bool skip_cells);
  //Helper functions
  void decrement_num_live_entries_in_live_cell_marker(const size_t row_query_idx);
  void initialize_live_cell_marker_from_tail(const size_t row_query_idx, const int64_t* coords,
      const int64_t END);
  /*
   * Advance indexes till a useful cell is found
   * Return true iff the current query interval has valid cells after increment
   * min_num_cells_to_increment - can be set to 0, in that case the "current" cell is checked
   * to see if it's useful. If not, the index is incremented.
   * If min_num_cells_to_increment >0,
   * the index is advanced first and then the cells underneath are checked for "usefulness"
   * A cell is "useful" iff
   *   its row is part of the queried rows
   *   in simple traversal mode - END >= begin
   *   in find intersecting intervals mode the row_query_idx is !initialized
   */
  bool advance_to_next_useful_cell(const uint64_t min_num_cells_to_increment);
 private:
  /*
   * num_cells_incremented - returns the actual #cells passed over
   */
  bool advance_coords_and_END_till_useful_cell_found(
    const uint64_t min_num_cells_to_increment,
    uint64_t& num_cells_incremented);
  bool advance_coords_and_END(const uint64_t num_cells_to_advance);
  void advance_fields_other_than_coords_END(const uint64_t num_cells_to_increment);
  bool keep_advancing_in_find_intersecting_intervals_mode() const;
  inline bool is_duplicate_cell_at_end_position(const int64_t coords_column_value, const int64_t END_field_value) const {
    return coords_column_value > END_field_value;
  }
  inline bool is_duplicate_cell_at_end_position_that_begins_before_query_interval(const int64_t coords_column_value,
      const int64_t END_field_value, const uint64_t query_column) const {
    //interval corresponding to this cell begins before the query column
    return is_duplicate_cell_at_end_position(coords_column_value, END_field_value)
           && (static_cast<uint64_t>(END_field_value) < query_column);
  }
  /*
   * Optimization for skipping useless cells
   * Useless are determined using the coords and the END value
   * a) in simple traversal mode, a cell is useless if END < coords[1] (END cell) OR
   * b) in find intersecting intervals mode, a cell is useless if the corresponding row_query_idx is
   * already initialized OR
   * c) the row is not part of the query
   * Once the number of contiguous useless cells are determined, fields other than coords and END
   * can be skipped in one shot
   */
  void increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
#ifdef DO_PROFILING
  void increment_num_cells_traversed_stats(const uint64_t num_cells_incremented);
#endif
  //Given the columnar field object, get the buffer pointer and index in the buffer
  //Response depends on whether this is a cell that begins before the query interval begin
  //or cell >= query interval begin
  inline std::pair<const GenomicsDBBuffer*,size_t> get_buffer_and_index(const int field_query_idx) const {
    assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[field_query_idx];
    //No more markers for intervals intersecting query interval begin
    //get data from live list tail
    if (m_PQ_live_cell_markers.empty())
      return std::pair<const GenomicsDBBuffer*,size_t>(genomicsdb_columnar_field.get_live_buffer_list_tail() ,
          genomicsdb_columnar_field.get_curr_index_in_live_list_tail());
    else {
      const auto row_query_idx = m_PQ_live_cell_markers.top().second;
      return std::pair<const GenomicsDBBuffer*,size_t>(m_live_cell_markers.get_buffer(field_query_idx, row_query_idx),
          m_live_cell_markers.get_index(field_query_idx, row_query_idx));
    }
  }
 protected:
  bool m_done_reading_from_TileDB;
  bool m_in_find_intersecting_intervals_mode;
  bool m_in_simple_traversal_mode;
  bool m_at_new_query_column_interval;
  //Flag that specifies if this is the first time reading from TileDB
  bool m_first_read_from_TileDB;
  unsigned m_END_query_idx;
  const VariantArraySchema* m_variant_array_schema;
  const VariantQueryConfig* m_query_config;
  uint64_t m_query_column_interval_idx;
  GenomicsDBColumnarCell* m_cell;
  //Buffers for fields
  std::vector<GenomicsDBColumnarField> m_fields;
  //Cell markers for handling the sweep operation
  GenomicsDBLiveCellMarker m_live_cell_markers;
  //Cell markers in column major order
  std::priority_queue<ColumnRowQueryIdxPair, std::vector<ColumnRowQueryIdxPair>,
      GenomicsDBLiveCellMarkerColumnMajorComparator> m_PQ_live_cell_markers;
  uint64_t m_num_queried_rows_initialized;
  //Contains query idx for only the fields that must be fetched from TileDB in the next round
  //The first time all fields are queried - in subsequent iterations only those fields whose
  //buffers are consumed completely are queried
  std::vector<int> m_query_attribute_idx_vec;
  //Contains #cells to skip - 1 for every entry in m_query_attribute_idx_vec
  std::vector<size_t> m_query_attribute_idx_num_cells_to_increment_vec;
  //Since variable length fields have buffers for offsets, need a mapping structure
  std::vector<size_t> m_query_attribute_idx_to_tiledb_buffer_idx;
  std::vector<void*> m_buffer_pointers;
  std::vector<size_t> m_buffer_sizes;
  std::vector<size_t> m_skip_counts;
  //The TileDB array object
  //Could point to an object initialized by the caller or could point
  //to an object owned by the iterator
  const TileDB_Array* m_tiledb_array;
  //If the iterator owns a TileDB array object, then tiledb_array_init() is
  //called for this object and m_tiledb_array points to it
  TileDB_Array* m_owned_tiledb_array;
#ifdef DO_PROFILING
  uint64_t m_num_cells_traversed_stats[GenomicsDBIteratorStatsEnum::NUM_STATS];
  //For a given query, tracks the length of consecutive cell segments that are useless
  std::vector<uint64_t> m_useless_cell_interval_lengths_histogram;
  std::vector<uint64_t> m_num_cells_traversed_in_find_intersecting_intervals_mode_histogram;
#ifdef COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW
  std::vector<uint64_t> m_cell_counts_since_last_cell_from_same_row;
  std::vector<uint64_t> m_histogram_cell_counts_since_last_cell_from_same_row;
#endif //COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW
#ifdef PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL
  std::vector<std::pair<int64_t, int64_t>> m_observed_cells_in_curr_window;
  std::vector<uint64_t> m_row_idx_to_num_observed_cells_in_curr_window;
  uint64_t m_num_observed_row_idxs_in_curr_window;
  uint64_t m_num_observed_cells_in_curr_window;
  void update_sliding_window_to_profile_num_cells_to_traverse(
    const GenomicsDBColumnarField& coords_columnar_field);
#endif //PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL
#endif //DO_PROFILING
 public:
  class valid_row_query_idx_iterator {
    public:
      valid_row_query_idx_iterator(const std::vector<bool>::const_iterator& iter, const std::vector<bool>::const_iterator& end,
          const size_t idx)
        : m_iter(iter), m_end(end), m_idx(idx) {
          if(m_iter != m_end && !(*m_iter)) {
            auto valid_iter = std::find(m_iter, m_end, true);
            m_idx += std::distance(m_iter, valid_iter);
            m_iter = valid_iter;
          }
        }
      bool operator !=(const valid_row_query_idx_iterator& other) const {
        return (m_iter != other.m_iter || m_idx != other.m_idx);
      }
      size_t operator *() const {
        return m_idx;
      }
      valid_row_query_idx_iterator& operator ++() {
        ++m_iter;
        ++m_idx;
        if(m_iter != m_end && !(*m_iter)) {
          auto valid_iter = std::find(m_iter, m_end, true);
          m_idx += std::distance(m_iter, valid_iter);
          m_iter = valid_iter;
        }
        return *this;
      }
    private:
      std::vector<bool>::const_iterator m_iter;
      std::vector<bool>::const_iterator m_end;
      size_t m_idx;
  };
  valid_row_query_idx_iterator begin_valid_row_query_idx() const {
    return valid_row_query_idx_iterator(m_live_cell_markers.valid_begin(), m_live_cell_markers.valid_end(), 0u);
  }
  valid_row_query_idx_iterator end_valid_row_query_idx() const {
    return valid_row_query_idx_iterator(m_live_cell_markers.valid_end(), m_live_cell_markers.valid_end(),
        m_query_config->get_num_rows_to_query());
  }
};

/*
 * For GVCFs, you need to maintain a structure which arranges samples
 * in increasing order of END values. Each sample would have at most 1 entry
 * in this set. Ideally, this would be a min heap. However, due to the problem
 * of intersecting intervals for a given sample,
 * (https://github.com/GenomicsDB/GenomicsDB/wiki/Overlapping-variant-calls-in-a-sample)
 * we might have to delete elements from the structure whose END values
 * are not the minimum. Hence, for a heap, we would need to pop till we find
 * the element to delete and re-push all the elements we popped.
 * This is not optimal (see lines 472-480 in query_variants.cc)
 * Instead of a heap, we use a tree (STL set) whose elements contain
 * (a) END position (b) row_query_idx
 * This allows us to perform 'pops' and the above deletions in O(log(n)) time
 * , albeit with a higher constant for pops relative to a heap.
 * See STL doc for find() function to determine how a set can do the deletion
 * in O(log(N)) time
 */
typedef std::pair<int64_t, size_t> GVCFEndSetElementTy;

class GVCFEndSetElementComparator {
 public:
  bool operator()(const GVCFEndSetElementTy& a, const GVCFEndSetElementTy& b) const {
    return (a.first < b.first) || (a.first == b.first && a.second < b.second);
  }
};

class GenomicsDBGVCFCell;
class GenomicsDBGVCFIterator : public SingleCellTileDBIterator {
  public:
    GenomicsDBGVCFIterator(TileDB_CTX* tiledb_ctx,
	const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
	const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size);
    //Uses tiledb_array object provided by caller (if non-NULL)
    GenomicsDBGVCFIterator(TileDB_CTX* tiledb_ctx,
	const TileDB_Array* tiledb_array,
	const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
	const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size);
    ~GenomicsDBGVCFIterator();
    //Delete default copy and move constructors
    GenomicsDBGVCFIterator(const GenomicsDBGVCFIterator& other) = delete;
    GenomicsDBGVCFIterator& operator=(const GenomicsDBGVCFIterator& other) = delete;
    GenomicsDBGVCFIterator(GenomicsDBGVCFIterator&& other) = delete;
    //iterator functions
    const GenomicsDBGVCFIterator& operator++();
    inline const GenomicsDBGVCFCell& operator*() const {
      return *m_cell;
    }
    inline bool end() const {
      return SingleCellTileDBIterator::end() && m_end_set.empty();
    }
    inline std::pair<const uint8_t*, size_t> get_raw_pointer_and_length_for_query_idx(const size_t row_query_idx,
	const int field_query_idx) const {
      assert(static_cast<size_t>(field_query_idx) < m_fields.size());
      auto& genomicsdb_columnar_field = m_fields[field_query_idx];
      const auto buffer_ptr = m_live_cell_markers.get_buffer(field_query_idx, row_query_idx);
      const auto index = m_live_cell_markers.get_index(field_query_idx, row_query_idx);
      return std::pair<const uint8_t*, size_t>(
	  genomicsdb_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(buffer_ptr, index),
	  genomicsdb_columnar_field.get_length_of_data_in_buffer_at_index(buffer_ptr, index));
    }
    //Returns current variant interval
    ColumnRange get_current_variant_interval() const {
      return ColumnRange(m_current_start_position, m_current_end_position);
    }
    //returns true if the row_query_idx is valid for current iteration (ie the sample
    //corresponding to row_query_idx
    bool is_valid_row_query_idx(const size_t row_query_idx) const {
      return m_live_cell_markers.is_valid(row_query_idx);
    }
    //row corresponding to idx must be valid before using this function
    bool is_field_valid_for_valid_row_query_idx(const unsigned field_idx, const size_t idx) const {
      auto gdb_buffer_ptr = m_live_cell_markers.get_buffer(field_idx, idx);
      auto index = m_live_cell_markers.get_index(field_idx, idx);
      return m_fields[field_idx].is_valid_data_in_buffer_at_index(gdb_buffer_ptr, index);
    }
    //Row corresponding to idx need not be valid before calling this function
    bool is_field_valid_for_row_query_idx(const unsigned field_idx, const size_t idx) const {
      return is_valid_row_query_idx(idx) && is_field_valid_for_valid_row_query_idx(field_idx, idx);
    }
  private:
    /*
     * new_cell = true implies a new call is being added to the live cell markers
     * false meaning call is being invalidated from the markers
     */
    template<bool new_cell>
    void update_num_deletions_and_MNVs(const size_t row_query_idx);
    /*
     * Keep filling m_end_set with cells from TileDB as long as the coords[1] == m_current_start_position
     * The moment you hit a cell with coords[1] > m_current_start_position, store its value
     * in m_next_start_position, update m_current_end_position and exit
     */
    void fill_end_set_in_simple_traversal_mode();
    /*
     * Given m_current_start_position, m_next_start_position and m_num_calls_with_deletions_or_MNVs,
     * update m_current_end_position
     */
    void update_current_end_position();
    /*
     * Should be called every time processing of a new query column interval begins.
     * Assumes that SingleCellTileDBIterator::begin_new_query_column_interval is already
     * called.
     */
    void begin_new_query_column_interval();
  private:
    int64_t m_current_start_position;
    int64_t m_current_end_position;
    int64_t m_next_start_position;
    int64_t m_query_interval_limit;
    uint64_t m_num_calls_with_deletions_or_MNVs;
    std::set<GVCFEndSetElementTy, GVCFEndSetElementComparator> m_end_set;
    unsigned m_REF_query_idx;
    unsigned m_ALT_query_idx;
    GenomicsDBGVCFCell* m_cell;
#ifdef DO_PROFILING
    std::vector<uint64_t> m_num_times_initialized;
    std::vector<uint64_t> m_num_times_invalidated;
    uint64_t m_bin_size;
#endif
};

#endif
