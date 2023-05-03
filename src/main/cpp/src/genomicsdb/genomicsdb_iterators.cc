/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
 * Copyright (c) 2020 Omics Data Automation, Inc.
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

#include "genomicsdb_iterators.h"
#include "variant_cell.h"
#include "variant_query_config.h"
#include "gt_remapper.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBIteratorException(#X);

#if defined(DO_PROFILING) && defined(COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW)
#define COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW_HISTOGRAM_BIN_SIZE 1ull
#endif

SingleCellTileDBIterator::SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx,
    const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
    const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
  : SingleCellTileDBIterator(tiledb_ctx,
                             0,
                             vid_mapper, variant_array_schema,
                             array_path, query_config, buffer_size)
{}

SingleCellTileDBIterator::SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx,
    const TileDB_Array* tiledb_array,
    const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
    const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
  : m_variant_array_schema(&variant_array_schema), m_query_config(&query_config),
    m_cell(new GenomicsDBColumnarCell(this)),
    m_tiledb_array(tiledb_array), m_owned_tiledb_array(0),
    m_query_column_interval_idx(0u),
    m_first_read_from_TileDB(true),
    m_done_reading_from_TileDB(false),
    m_in_find_intersecting_intervals_mode(false),
    m_in_simple_traversal_mode(false),
    m_at_new_query_column_interval(true),
    m_live_cell_markers(query_config.get_num_queried_attributes()+1u, query_config.get_num_rows_to_query()), //+1 for coords
    m_num_queried_rows_initialized(0),
    m_PQ_live_cell_markers() {
#ifdef DO_PROFILING
  memset(m_num_cells_traversed_stats, 0, GenomicsDBIteratorStatsEnum::NUM_STATS*sizeof(uint64_t));
  //No point in tracking beyond 10K cell segments
  m_useless_cell_interval_lengths_histogram.resize(10000u, 0ull);
  m_num_cells_traversed_in_find_intersecting_intervals_mode_histogram.resize(query_config.get_num_rows_to_query(), 0ull);
#ifdef COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW
  m_cell_counts_since_last_cell_from_same_row.resize(query_config.get_num_rows_to_query(), 0ull);
#endif //COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW
#ifdef PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL
  m_num_observed_row_idxs_in_curr_window = 0ull;
  m_num_observed_cells_in_curr_window = 0ull;
  m_row_idx_to_num_observed_cells_in_curr_window.resize(query_config.get_num_rows_to_query(), 0ull);
#endif //PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL
#endif //DO_PROFILING
  const auto& attribute_ids = query_config.get_query_attributes_schema_idxs();
  std::vector<const char*> attribute_names(attribute_ids.size()+1u);  //+1 for the COORDS
  m_query_attribute_idx_vec.resize(attribute_ids.size()+1u);//+1 for the COORDS
  m_query_attribute_idx_num_cells_to_increment_vec.resize(attribute_ids.size()+1u, 0u);
  m_query_attribute_idx_to_tiledb_buffer_idx.resize(attribute_ids.size()+1u);//+1 for the COORDS
  for (auto i=0ull; i<attribute_ids.size(); ++i) {
    auto schema_idx = attribute_ids[i];
    attribute_names[i] = m_variant_array_schema->attribute_name(schema_idx).c_str();
    m_query_attribute_idx_to_tiledb_buffer_idx[i] = m_buffer_pointers.size();
    auto is_variable_length_field = m_variant_array_schema->is_variable_length_field(schema_idx);
    //TileDB does not have a way to distinguish between char, string and int8_t fields
    //Hence, a series of possibly messy checks here
    assert(vid_mapper);
    const auto* vid_field_info = vid_mapper->get_field_info(attribute_names[i]);
    auto curr_type_index = (vid_field_info
                            && vid_field_info->get_genomicsdb_type().get_tuple_element_bcf_ht_type(0u) == BCF_HT_FLAG)
                           ? std::type_index(typeid(bool))
                           : m_variant_array_schema->type(schema_idx);
    //GenomicsDBColumnarField
    m_fields.emplace_back(curr_type_index,
                          is_variable_length_field ? BCF_VL_VAR : BCF_VL_FIXED,
                          m_variant_array_schema->val_num(schema_idx), buffer_size);
    //Buffer pointers and size
    m_buffer_pointers.push_back(0);
    m_buffer_sizes.push_back(0);
    m_skip_counts.push_back(0u);
    if (is_variable_length_field) {
      m_buffer_pointers.push_back(0);
      m_buffer_sizes.push_back(0);
      m_skip_counts.push_back(0u);
    }
  }
  //Co-ordinates
  auto coords_idx = attribute_ids.size();
  attribute_names[coords_idx] = TILEDB_COORDS;
  m_query_attribute_idx_to_tiledb_buffer_idx[coords_idx] = m_buffer_pointers.size();
  //GenomicsDBColumnarField
  m_fields.emplace_back(m_variant_array_schema->dim_type(), BCF_VL_FIXED, m_variant_array_schema->dim_length(), buffer_size);
  //Buffer pointers and size for COORDS
  m_buffer_pointers.push_back(0);
  m_buffer_sizes.push_back(0);
  m_skip_counts.push_back(0u);
  //Set row idx for live cell markers
  for (auto i=0ull; i<query_config.get_num_rows_to_query(); ++i)
    m_live_cell_markers.set_row_idx(i, query_config.get_array_row_idx_for_query_row_idx(i));
  //END query idx
  assert(m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  m_END_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_END_IDX);
  assert(m_END_query_idx < m_fields.size() && m_END_query_idx == 0u); //END must always be the first field
  //first read
  begin_new_query_column_interval(tiledb_ctx, array_path.c_str(), &attribute_names);
}

SingleCellTileDBIterator::~SingleCellTileDBIterator() {
  if (m_cell)
    delete m_cell;
  m_cell = 0;
  if (m_owned_tiledb_array)
    tiledb_array_finalize(m_owned_tiledb_array);
  m_owned_tiledb_array = 0;
  m_tiledb_array = 0;
#ifdef DO_PROFILING
  std::cerr << "Iterator stats";
  for (auto i=0u; i<GenomicsDBIteratorStatsEnum::NUM_STATS; ++i)
    std::cerr << "," << m_num_cells_traversed_stats[i];
  std::cerr << "\n";
  std::cerr << "Useless cell segment lengths histogram\n";
  for (auto i=0ull; i<m_useless_cell_interval_lengths_histogram.size(); ++i)
    std::cerr << i << "  " <<m_useless_cell_interval_lengths_histogram[i]<<"\n";
  auto num_live_list_entries = 0ull;
  auto num_free_list_entries = 0ull;
  for (const auto& field : m_fields) {
    num_free_list_entries += field.get_free_buffer_list_length();
    num_live_list_entries += field.get_live_buffer_list_length();
  }
  std::cerr << "Buffer_lists_lengths final "<<num_free_list_entries
            <<" "<<num_live_list_entries <<"\n";
  std::cerr << "Histogram:\n";
  for (auto val : m_num_cells_traversed_in_find_intersecting_intervals_mode_histogram)
    std::cerr << val << "\n";
#ifdef COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW
  std::cerr << "COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW histogram\n";
  for (auto i=0ull; i<m_histogram_cell_counts_since_last_cell_from_same_row.size(); ++i)
    std::cerr << COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW_HISTOGRAM_BIN_SIZE*(i+1ull)
              << " " << m_histogram_cell_counts_since_last_cell_from_same_row[i] << "\n";
  m_histogram_cell_counts_since_last_cell_from_same_row.clear();
  m_cell_counts_since_last_cell_from_same_row.clear();
#endif //COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW_HISTOGRAM_BIN_SIZE
#ifdef PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL
  //Final window
  if (m_num_observed_cells_in_curr_window > 0u)
    std::cerr << m_observed_cells_in_curr_window[0].second << '\t'
              << INT64_MAX-1 << '\t' << m_num_observed_cells_in_curr_window << '\n';
#endif //PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL
#endif //DO_PROFILING
}

inline bool SingleCellTileDBIterator::keep_advancing_in_find_intersecting_intervals_mode() const {
  return (!m_done_reading_from_TileDB
          && (m_num_queried_rows_initialized < m_query_config->get_num_rows_to_query()));
}

void SingleCellTileDBIterator::begin_new_query_column_interval(TileDB_CTX* tiledb_ctx, const char* array_path,
    std::vector<const char*>* attribute_names) {
  do {
    //This can happen iff in the previous iteration intersecting intervals were searched OR
    //full array scan is requested and this is the first time a read is being performed OR
    //if only a single pass simple traversal is desired
    if (m_in_find_intersecting_intervals_mode
        || (m_first_read_from_TileDB && (m_query_config->get_num_column_intervals() == 0u
                                         || m_query_config->bypass_intersecting_intervals_phase())))
      move_into_simple_traversal_mode();
    else {
      reset_for_next_query_interval();
      //If first read, do the read without worrying about m_query_column_interval_idx
      if (!m_first_read_from_TileDB) {
        //Move to next interval
        ++m_query_column_interval_idx;
        //If no more query intervals remaining - exit
        if (m_query_column_interval_idx >= m_query_config->get_num_column_intervals()) {
          m_done_reading_from_TileDB = true;
          return;
        }
      }
      m_in_simple_traversal_mode = m_query_config->bypass_intersecting_intervals_phase();
      m_in_find_intersecting_intervals_mode = !m_query_config->bypass_intersecting_intervals_phase();
    }
    auto smallest_row_idx_in_array = m_query_config->get_smallest_row_idx_in_array();
    //Range of the query
    int64_t query_range[4] = { smallest_row_idx_in_array,
                               static_cast<int64_t>(m_query_config->get_num_rows_in_array()+smallest_row_idx_in_array-1), //inclusive bounds
                               0, INT64_MAX-1
                             };
    if (m_query_config->get_num_column_intervals() > 0u) {
      assert(m_query_column_interval_idx < m_query_config->get_num_column_intervals());
      query_range[2] = m_query_config->get_column_begin(m_query_column_interval_idx);
      if (m_in_simple_traversal_mode) {
        query_range[3] = m_query_config->get_column_end(m_query_column_interval_idx);
      }
    }
    auto status = -1;
    if (m_tiledb_array == 0) {
      assert(array_path && attribute_names && tiledb_ctx);
      /* Initialize the array in READ mode. */
      status = tiledb_array_init(
                 tiledb_ctx,
                 &m_owned_tiledb_array,
                 array_path,
                 TILEDB_ARRAY_READ,
                 reinterpret_cast<const void*>(query_range),
                 &((*attribute_names)[0]),
                 attribute_names->size());
      if (status != TILEDB_OK)
        throw GenomicsDBIteratorException(std::string("Error while initializing TileDB array object")
                                          + "\nTileDB error message : "+tiledb_errmsg);
      m_tiledb_array = m_owned_tiledb_array;
    } else {
      //TileDB_Array object is provided by caller and this is the first read, reset attributes to query
      if (m_first_read_from_TileDB) {
        assert(attribute_names);
        status = tiledb_array_reset_attributes(m_tiledb_array,
                                               &((*attribute_names)[0]),
                                               attribute_names->size());
        if (status == TILEDB_OK && !m_query_config->get_query_filter().empty()) {
          status = tiledb_array_apply_filter(m_tiledb_array, m_query_config->get_query_filter().c_str());
        }
        if (status != TILEDB_OK)
          throw GenomicsDBIteratorException(std::string("Error while initializing attributes for the TileDB array object")
                                            + "\nTileDB error message : "+tiledb_errmsg);
      }
      status = tiledb_array_reset_subarray(m_tiledb_array, reinterpret_cast<const void*>(query_range));
      if (status != TILEDB_OK)
        throw GenomicsDBIteratorException(std::string("Error in tiledb_array_reset_subarray()")
                                          + "\nTileDB error message : "+tiledb_errmsg);
    }
    read_from_TileDB(false);
    m_first_read_from_TileDB = false;
    if (!m_done_reading_from_TileDB) {
      //This is useful when a subset of rows is queried - skip rows not part of query
      advance_to_next_useful_cell(0u);
      //Dealing with new query interval begin
      if (m_in_find_intersecting_intervals_mode) {
#ifdef DO_PROFILING
        auto num_useless_cells_traversed =
          m_num_cells_traversed_stats[GenomicsDBIteratorStatsEnum::NUM_USELESS_CELLS_IN_FIND_INTERSECTING_INTERVALS_MODE];
#endif
        //Array done or all initialized - exit loop
        while (keep_advancing_in_find_intersecting_intervals_mode()) {
#ifdef DO_PROFILING
          assert(m_num_queried_rows_initialized < m_query_config->get_num_rows_to_query());
          m_num_cells_traversed_in_find_intersecting_intervals_mode_histogram[m_num_queried_rows_initialized]
          += (m_num_cells_traversed_stats[GenomicsDBIteratorStatsEnum::NUM_USELESS_CELLS_IN_FIND_INTERSECTING_INTERVALS_MODE]
              - num_useless_cells_traversed);
          num_useless_cells_traversed =
            m_num_cells_traversed_stats[GenomicsDBIteratorStatsEnum::NUM_USELESS_CELLS_IN_FIND_INTERSECTING_INTERVALS_MODE];
#endif
          handle_current_cell_in_find_intersecting_intervals_mode();
          if (keep_advancing_in_find_intersecting_intervals_mode())
            advance_to_next_useful_cell(1u);
        }
#ifdef DO_PROFILING
        auto num_live_list_entries = 0ull;
        auto num_free_list_entries = 0ull;
        for (const auto& field : m_fields) {
          num_free_list_entries += field.get_free_buffer_list_length();
          num_live_list_entries += field.get_live_buffer_list_length();
        }
        logger.info("Buffer_lists_lengths end_of_intersecting_intervals_mode {} {}", num_free_list_entries, num_live_list_entries);
#endif
      }
    }
    //If done reading current column interval, move to next query column interval, if any remaining OR
    //if in find intersecting intervals mode, switch to simple traversal in the next iteration
  } while ((m_done_reading_from_TileDB && m_PQ_live_cell_markers.empty()
            && (m_query_column_interval_idx+1u) < m_query_config->get_num_column_intervals())
           || m_in_find_intersecting_intervals_mode);
}

void SingleCellTileDBIterator::move_into_simple_traversal_mode() {
  //Move into simple traversal mode
  m_done_reading_from_TileDB = false;
  m_in_find_intersecting_intervals_mode = false;
  m_in_simple_traversal_mode = true;
  //Query all attributes
  m_query_attribute_idx_vec.resize(m_fields.size());
  for (auto i=0u; i<m_fields.size(); ++i)
    m_query_attribute_idx_vec[i] = i;
}

void SingleCellTileDBIterator::reset_for_next_query_interval() {
  m_done_reading_from_TileDB = false;
  m_in_find_intersecting_intervals_mode = false;
  m_in_simple_traversal_mode = false;
  m_at_new_query_column_interval = true;
  //Reset markers
  assert(m_PQ_live_cell_markers.empty());
  m_live_cell_markers.reset();
  m_num_queried_rows_initialized = 0u;
  //Move all buffers to free list - moving to the next query interval
  for (auto& field : m_fields)
    field.move_all_buffers_from_live_list_to_free_list();
  //Query all attributes
  m_query_attribute_idx_vec.resize(m_fields.size());
  for (auto i=0u; i<m_fields.size(); ++i)
    m_query_attribute_idx_vec[i] = i;
}

void SingleCellTileDBIterator::read_from_TileDB(const bool skip_cells) {
  //Zero out all buffer sizes
  memset(&(m_buffer_sizes[0]), 0, m_buffer_sizes.size()*sizeof(size_t));
  //Zero out skip counts
  m_skip_counts.assign(m_skip_counts.size(), 0u);
  //Only set non-0 buffer sizes for attributes for which we wish to get more data
  for (auto i=0ull; i<m_query_attribute_idx_vec.size(); ++i) {
    auto query_idx = m_query_attribute_idx_vec[i];
    assert(static_cast<size_t>(query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    //Get a free buffer and move it to the live list - buffer gets added to tail of live list
    genomicsdb_columnar_field.move_buffer_to_live_list(genomicsdb_columnar_field.get_or_allocate_free_buffer());
    //Always read new data into tail
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail();
    assert(genomicsdb_buffer_ptr);
    assert(static_cast<size_t>(query_idx) < m_query_attribute_idx_to_tiledb_buffer_idx.size());
    auto buffer_idx = m_query_attribute_idx_to_tiledb_buffer_idx[query_idx];
    assert(buffer_idx < m_buffer_pointers.size());
    auto skip_count = skip_cells ? m_query_attribute_idx_num_cells_to_increment_vec[i] : 0ull;
    //For variable length field, first the offsets buffer
    if (genomicsdb_columnar_field.is_variable_length_field()) {
      m_buffer_pointers[buffer_idx] = reinterpret_cast<void*>(genomicsdb_buffer_ptr->get_offsets_pointer());
      m_buffer_sizes[buffer_idx] = genomicsdb_buffer_ptr->get_offsets_size_in_bytes();
      m_skip_counts[buffer_idx] = skip_count;
      ++buffer_idx;
      assert(buffer_idx < m_buffer_pointers.size());
    }
    m_buffer_pointers[buffer_idx] = reinterpret_cast<void*>(genomicsdb_buffer_ptr->get_raw_pointer());
    m_buffer_sizes[buffer_idx] = genomicsdb_buffer_ptr->get_data_vector_size_in_bytes();
    m_skip_counts[buffer_idx] = skip_count;
  }
  auto status = tiledb_array_skip_and_read(m_tiledb_array, &(m_buffer_pointers[0]), &(m_buffer_sizes[0]), &(m_skip_counts[0u]));
  if (status != TILEDB_OK)
    throw GenomicsDBIteratorException(std::string("Error while reading from TileDB array ")
                                      + "\nTileDB error message : "+tiledb_errmsg);
#ifdef DEBUG
  auto num_done_fields = 0u;
  //tiledb_array_skip_and_read must have skipped all cells that needed to be skipped
  for (auto val : m_skip_counts)
    assert(val == 0u);
#endif
  //Set number of live entries in each buffer
  for (auto i=0ull; i<m_query_attribute_idx_vec.size(); ++i) {
    auto query_idx = m_query_attribute_idx_vec[i];
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail();
    auto buffer_idx = m_query_attribute_idx_to_tiledb_buffer_idx[query_idx];
    //tiledb_array_skip_and_read() must have taken care of this
    assert(m_skip_counts[buffer_idx] == 0u);
    m_query_attribute_idx_num_cells_to_increment_vec[i] = 0u; //index i, not query_idx
    //For variable length field, first the offsets buffer
    auto filled_buffer_size = m_buffer_sizes[buffer_idx];
    auto num_live_entries = 0u;
    if (genomicsdb_columnar_field.is_variable_length_field()) {
      num_live_entries = filled_buffer_size/sizeof(size_t);
      genomicsdb_buffer_ptr->set_num_live_entries(num_live_entries);
      //Append size of the variable length field
      genomicsdb_buffer_ptr->set_offset(num_live_entries, m_buffer_sizes[buffer_idx+1u]);
      //tiledb_array_skip_and_read() must have taken care of this
      assert(m_skip_counts[buffer_idx+1u] == 0u);
    } else {
      num_live_entries = filled_buffer_size/
                         (genomicsdb_columnar_field.get_fixed_length_field_size_in_bytes());
      genomicsdb_buffer_ptr->set_num_live_entries(num_live_entries);
    }
    //Marks data as valid/invalid
    genomicsdb_columnar_field.set_valid_vector_in_live_buffer_list_tail_ptr();
    if (num_live_entries == 0u) {
#ifdef DEBUG
      ++num_done_fields;
#endif
      m_done_reading_from_TileDB = true;
      genomicsdb_columnar_field.move_buffer_to_free_list(genomicsdb_buffer_ptr);
    }
  }
#ifdef DEBUG
  //The way the code is written, the final call to read_from_TileDB() will attempt to
  //fetch data for all queried fields. TileDB should return 0 for all fields. Hence, either all
  //queried fields return with 0 size or have some data to process
  assert(num_done_fields == 0u || (num_done_fields == m_query_attribute_idx_vec.size() && m_done_reading_from_TileDB));
#endif
}

void SingleCellTileDBIterator::initialize_live_cell_marker_from_tail(const size_t row_query_idx, const int64_t* coords,
    const int64_t END_field_value) {
  const auto row_idx = coords[0];
  const auto coords_column = coords[1];
  assert(m_query_config->is_queried_array_row_idx(row_idx));
  assert(m_live_cell_markers.get_row_idx(row_query_idx) == row_idx);
  //Set initialized
  m_live_cell_markers.set_initialized(row_query_idx, true);
  m_live_cell_markers.set_valid(row_query_idx, true);
  m_live_cell_markers.set_column_interval(row_query_idx, coords_column, END_field_value);
  //Track offsets for each field and keep the buffer alive
  for (auto i=0u; i<m_fields.size(); ++i) {
    auto& genomicsdb_columnar_field = m_fields[i];
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail();
    assert(genomicsdb_buffer_ptr);
    m_live_cell_markers.set_field_marker(i, row_query_idx, genomicsdb_buffer_ptr,
                                         genomicsdb_columnar_field.get_curr_index_in_live_list_tail());
    genomicsdb_buffer_ptr->increment_num_live_entries();
  }
}

void SingleCellTileDBIterator::decrement_num_live_entries_in_live_cell_marker(const size_t row_query_idx) {
  assert(m_live_cell_markers.is_valid(row_query_idx) && m_live_cell_markers.is_initialized(row_query_idx));
  for (auto i=0u; i<m_fields.size(); ++i) {
    auto* buffer_ptr = m_live_cell_markers.get_buffer(i, row_query_idx);
    buffer_ptr->decrement_num_live_entries(1u);
    if (buffer_ptr->get_num_live_entries() == 0u)
      m_fields[i].move_buffer_to_free_list(buffer_ptr);
  }
  m_live_cell_markers.set_valid(row_query_idx, false);
  m_live_cell_markers.set_initialized(row_query_idx, false);
}

void SingleCellTileDBIterator::handle_current_cell_in_find_intersecting_intervals_mode() {
  assert(m_in_find_intersecting_intervals_mode && !m_done_reading_from_TileDB);
  auto& coords_columnar_field = m_fields[m_fields.size()-1u];
  auto* coords = reinterpret_cast<int64_t*>(
                   coords_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
                     coords_columnar_field.get_live_buffer_list_tail(),
                     coords_columnar_field.get_curr_index_in_live_list_tail()
                   )
                 );
  auto& END_columnar_field = m_fields[m_END_query_idx];
  auto& END_field_value = *(reinterpret_cast<int64_t*>(
                              END_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
                                END_columnar_field.get_live_buffer_list_tail(),
                                END_columnar_field.get_curr_index_in_live_list_tail()
                              )
                            ));
  auto row_idx = coords[0];
  auto row_query_idx = m_query_config->get_query_row_idx_for_array_row_idx(row_idx);
  //should not be initialize already
  assert(!(m_live_cell_markers.is_initialized(row_query_idx)));
  //Reference - value gets modified
  auto& coords_column = coords[1];
  assert(m_query_column_interval_idx < m_query_config->get_num_column_intervals());
  auto query_interval_begin = m_query_config->get_column_begin(m_query_column_interval_idx);
  //Only deal with intervals that begin before query_interval_begin and intersect it
  assert(is_duplicate_cell_at_end_position_that_begins_before_query_interval(coords_column,
         END_field_value, query_interval_begin));
  std::swap(coords_column, END_field_value);
  initialize_live_cell_marker_from_tail(row_query_idx, coords, END_field_value);
  //Add to PQ
  m_PQ_live_cell_markers.push(ColumnRowQueryIdxPair(coords[1], row_query_idx));
  ++m_num_queried_rows_initialized;
}

const SingleCellTileDBIterator& SingleCellTileDBIterator::operator++() {
  //Intervals intersecting the column begin still exist in PQ
  if (!m_PQ_live_cell_markers.empty()) {
    auto row_query_idx = m_PQ_live_cell_markers.top().second;
    m_PQ_live_cell_markers.pop();
    //Decrease #live entries and move to free list if needed
    decrement_num_live_entries_in_live_cell_marker(row_query_idx);
    m_at_new_query_column_interval = false;
    if (!m_PQ_live_cell_markers.empty())
      return *this;
    //Done processing PQ cells
    //Note that begin_new_query_column_interval() moves into simple traversal mode
    //immediately after processing cells that intersect with query column interval begin
    //If m_done_reading_from_TileDB, this means that no cells in the current query interval
    //were available in the simple traversal mode. Move to next column interval
    if (m_done_reading_from_TileDB)
      begin_new_query_column_interval();
  } else {
    assert(m_in_simple_traversal_mode);
    m_at_new_query_column_interval = false;
    auto increment_done = advance_to_next_useful_cell(1u);
    //TileDB couldn't provide > 1 cell, must move to next column interval
    if (!increment_done)
      begin_new_query_column_interval();
  }
  return *this;
}

bool SingleCellTileDBIterator::advance_to_next_useful_cell(const uint64_t min_num_cells_to_increment) {
  uint64_t num_cells_incremented = 0ull;
  auto increment_done = advance_coords_and_END_till_useful_cell_found(min_num_cells_to_increment,
                        num_cells_incremented);
  //Don't bother with rest of the fields since TileDB couldn't provide >= min_num_cells_to_increment
  if (!increment_done)
    return false;
  if (num_cells_incremented > 0u) {
#ifdef DO_PROFILING
    increment_num_cells_traversed_stats(num_cells_incremented);
#endif
    advance_fields_other_than_coords_END(num_cells_incremented);
  }

  // Evaluate filter expression if necessary
  m_cell_evaluated_with_filter_expression = true;
  if (m_query_config->get_query_filter().size() > 0) {
    const auto& attribute_ids = m_query_config->get_query_attributes_schema_idxs();
    std::vector<void *> buffers;
    std::vector<size_t> buffer_sizes;
    std::vector<int64_t> positions;
    for (auto i=0; i<m_fields.size(); i++) {
      auto& genomicsdb_columnar_field = m_fields[i];
      auto buffer = get_buffer_and_index(i);
      if (genomicsdb_columnar_field.is_variable_length_field()) {
        buffers.push_back(const_cast<GenomicsDBBuffer *>(buffer.first)->get_offsets_pointer());
        buffer_sizes.push_back(buffer.first->get_offsets_size_in_bytes());
        positions.push_back((int64_t)buffer.second);
      }
      buffers.push_back(const_cast<uint8_t *>(buffer.first->get_raw_pointer()));
      buffer_sizes.push_back(buffer.first->get_data_vector_size_in_bytes());
      positions.push_back((int64_t)buffer.second);
    }
    if (tiledb_array_evaluate_cell(m_tiledb_array, buffers.data(), buffer_sizes.data(), positions.data()) != TILEDB_OK) {
      m_cell_evaluated_with_filter_expression = false;
    }
  }

  return true;
}

//Advance iterator till a useful cell is found
//A cell is "useful" iff
//its row is part of the queried rows && (
//(in simple traversal mode - END >= begin) ||
//(in find intersecting intervals mode - !initialized))
bool SingleCellTileDBIterator::advance_coords_and_END_till_useful_cell_found(
  const uint64_t min_num_cells_to_increment,
  uint64_t& num_cells_incremented) {
  num_cells_incremented = 0ull;
  auto hitting_useless_cells = true;
  const auto& coords_columnar_field = m_fields[m_fields.size()-1u];
  assert(m_END_query_idx < m_fields.size());
  const auto& END_columnar_field = m_fields[m_END_query_idx];
  auto num_cells_to_advance_in_next_iteration = min_num_cells_to_increment;
  const auto subset_of_rows_queried = m_query_config->get_num_rows_to_query()
    < m_query_config->get_num_rows_in_array();
  while (hitting_useless_cells && !m_done_reading_from_TileDB) {
    //TODO: opportunity for vectorization
    //in simple mode, create bitvector Columnar[COORDS[1]] <= Columnar[END]
    //find first bit != 0, num_cells_to_advance == index of first bit

#if defined(DO_PROFILING) && defined(PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL)
    update_sliding_window_to_profile_num_cells_to_traverse(coords_columnar_field);
#endif
    auto increment_done = advance_coords_and_END(num_cells_to_advance_in_next_iteration);
    //TileDB couldn't provide > num_cells_to_advance_in_next_iteration
    if (!increment_done)
      return false;
    num_cells_incremented += num_cells_to_advance_in_next_iteration;
    num_cells_to_advance_in_next_iteration = 1u; //single cell stepping
    const auto* coords = reinterpret_cast<const int64_t*>(
                           coords_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
                             coords_columnar_field.get_live_buffer_list_tail(),
                             coords_columnar_field.get_curr_index_in_live_list_tail()
                           )
                         );
    //std::cerr << "COORDS "<<coords[0] <<" "<<coords[1] << "\n";
    //TileDB doesn't have a good way of requesting a subset of rows
    //Skip over rows that are not part of the query if coords are part of the queried attributes
    //in this round
    if (subset_of_rows_queried && !(m_query_config->is_queried_array_row_idx(coords[0])))
      hitting_useless_cells = true;
    else {
      auto coords_column = coords[1];
      assert(END_columnar_field.get_live_buffer_list_tail()->get_num_unprocessed_entries() > 0u);
      auto END_field_value = *(reinterpret_cast<const int64_t*>(
                                 END_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
                                   END_columnar_field.get_live_buffer_list_tail(),
                                   END_columnar_field.get_curr_index_in_live_list_tail()
                                 )
                               ));
      if (m_in_simple_traversal_mode) {
        //keep incrementing iterator if hitting duplicates at end in simple traversal mode
        hitting_useless_cells = is_duplicate_cell_at_end_position(coords_column, END_field_value);
#if defined(DO_PROFILING) && defined(COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW)
        auto row_query_idx = m_query_config->get_query_row_idx_for_array_row_idx(coords[0]);
        auto histogram_bin_idx = (m_cell_counts_since_last_cell_from_same_row[row_query_idx]/
                                  COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW_HISTOGRAM_BIN_SIZE);
        if (histogram_bin_idx >= m_histogram_cell_counts_since_last_cell_from_same_row.size())
          m_histogram_cell_counts_since_last_cell_from_same_row.resize(histogram_bin_idx+1ull, 0ull);
        ++(m_histogram_cell_counts_since_last_cell_from_same_row[histogram_bin_idx]);
        for (auto i=0ull; i<m_cell_counts_since_last_cell_from_same_row.size(); ++i)
          ++(m_cell_counts_since_last_cell_from_same_row[i]);
        m_cell_counts_since_last_cell_from_same_row[row_query_idx] = 0ull;
#endif
      } else {
        assert(m_in_find_intersecting_intervals_mode);
        auto row_idx = coords[0];
        assert(row_idx >= m_query_config->get_smallest_row_idx_in_array());
        auto row_query_idx = m_query_config->get_query_row_idx_for_array_row_idx(row_idx);
        assert(m_live_cell_markers.get_row_idx(row_query_idx) == row_idx);
        //if the row is already initialized in the find intersecting intervals mode, can skip
        //this cell
        if (m_live_cell_markers.is_initialized(row_query_idx))
          hitting_useless_cells = true;
        else {
          assert(m_query_column_interval_idx < m_query_config->get_num_column_intervals());
          auto query_interval_begin = m_query_config->get_column_begin(m_query_column_interval_idx);
          //Interval that begins before the queried column interval begin and intersects it
          //This is a valid interval for the query and must be handled by
          //handle_current_cell_in_find_intersecting_intervals_mode
          if (is_duplicate_cell_at_end_position_that_begins_before_query_interval(coords_column, END_field_value,
              query_interval_begin)) {
            //do not set initialized to true here, handle_current_cell_in_find_intersecting_intervals_mode will
            //do it
            hitting_useless_cells = false;
          } else {
            m_live_cell_markers.set_initialized(row_query_idx, true); //Non-intersecting cell, mark initialized
            ++m_num_queried_rows_initialized;
            hitting_useless_cells = keep_advancing_in_find_intersecting_intervals_mode();
          }
        }
      }
    }
  }
  //Increment was successful iff didn't hit the end of the array
  return !m_done_reading_from_TileDB;
}

//Advance iterator for coords and END fields
bool SingleCellTileDBIterator::advance_coords_and_END(const uint64_t num_cells_to_advance) {
  auto curr_query_column_interval_idx = m_query_column_interval_idx;
  if (num_cells_to_advance == 0u)
    return true;
  //2 = coords and END
  m_query_attribute_idx_vec.resize(2u); //no heap operations here
  m_query_attribute_idx_num_cells_to_increment_vec.resize(2u);
  m_query_attribute_idx_vec[0u] = m_fields.size()-1u; //coords
  m_query_attribute_idx_vec[1u] = m_END_query_idx;
  m_query_attribute_idx_num_cells_to_increment_vec[0u] = num_cells_to_advance;
  m_query_attribute_idx_num_cells_to_increment_vec[1u] = num_cells_to_advance;
  increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
  //some fields have exhausted buffers, need to fetch from TileDB
  while (!m_query_attribute_idx_vec.empty()) {
    read_from_TileDB(false);
    assert(curr_query_column_interval_idx == m_query_column_interval_idx);
    //read_from_TileDB() might have determined that no more cells exist for the current query
    //interval
    if (m_done_reading_from_TileDB)
      return false;
    increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
  }
  return true;
}

void SingleCellTileDBIterator::advance_fields_other_than_coords_END(const uint64_t num_cells_to_increment) {
  auto curr_query_column_interval_idx = m_query_column_interval_idx;
  //Increment iterator for other fields
  //All fields other than coords and END
  m_query_attribute_idx_vec.resize(m_fields.size()-2u);
  //END must be the first field
  assert(m_END_query_idx == 0u);
  for (auto i=0u; i<m_query_attribute_idx_vec.size(); ++i)
    m_query_attribute_idx_vec[i] = i+1u;  //END is the first field - ignore
  //-2 - ignore coords as well
  m_query_attribute_idx_num_cells_to_increment_vec.resize(m_fields.size()-2u);
  //For all fields, #cells to skip == num_cells_to_increment initially
  m_query_attribute_idx_num_cells_to_increment_vec.assign(
    m_query_attribute_idx_num_cells_to_increment_vec.size(), num_cells_to_increment);
  increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
  //some fields have exhausted buffers, need to fetch from TileDB
  while (!m_query_attribute_idx_vec.empty()) {
    read_from_TileDB(true);
    //TODO: either all fields are done reading or none are - is that right?
    //also, same column interval being queried
    assert(!m_done_reading_from_TileDB
           && m_query_column_interval_idx == curr_query_column_interval_idx);
#ifdef DEBUG
    //After the TileDB read with skip cells, there shouldn't be any attributes
    //whose cells need to be skipped after a call to read_from_TileDB()
    for (auto i=0u; i<m_query_attribute_idx_vec.size(); ++i)
      assert(m_query_attribute_idx_num_cells_to_increment_vec[i] == 0u);
#endif
    m_query_attribute_idx_vec.resize(0u); //no more attributes to query after skip
    //increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
  }
#ifdef DEBUG
  for (auto i=0u; i<m_fields.size(); ++i)
    assert(m_fields[i].get_live_buffer_list_tail()
           && m_fields[i].get_live_buffer_list_tail()->get_num_unprocessed_entries());
#endif
}

void SingleCellTileDBIterator::print(const int field_query_idx, std::ostream& fptr) const {
  assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[field_query_idx];
  auto x = get_buffer_and_index(field_query_idx);
  genomicsdb_columnar_field.print_data_in_buffer_at_index(fptr, x.first, x.second);
}

void SingleCellTileDBIterator::increment_iterator_within_live_buffer_list_tail_ptr_for_fields() {
  auto next_iteration_num_query_attributes = 0u;
  for (auto i=0u; i<m_query_attribute_idx_vec.size(); ++i) {
    assert(next_iteration_num_query_attributes <= i);
    assert(static_cast<size_t>(m_query_attribute_idx_vec[i]) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[m_query_attribute_idx_vec[i]];
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail();
    //TODO: either all fields are done reading or none are - is that right?
    if (genomicsdb_buffer_ptr) {
      auto unprocessed_entries_delta = std::min(m_query_attribute_idx_num_cells_to_increment_vec[i],
                                       genomicsdb_buffer_ptr->get_num_unprocessed_entries());
      genomicsdb_columnar_field.advance_curr_index_in_live_list_tail(unprocessed_entries_delta);
      genomicsdb_buffer_ptr->decrement_num_live_entries(unprocessed_entries_delta);
      genomicsdb_buffer_ptr->decrement_num_unprocessed_entries(unprocessed_entries_delta);
      if (genomicsdb_buffer_ptr->get_num_live_entries() == 0ull) //no more live entries, move buffer to free list
        genomicsdb_columnar_field.move_buffer_to_free_list(genomicsdb_buffer_ptr);
      //buffer completely processed, add attribute to next round of fetch from TileDB
      if (genomicsdb_buffer_ptr->get_num_unprocessed_entries() == 0ull) {
        m_query_attribute_idx_vec[next_iteration_num_query_attributes]
          = m_query_attribute_idx_vec[i];
        m_query_attribute_idx_num_cells_to_increment_vec[next_iteration_num_query_attributes]
          = m_query_attribute_idx_num_cells_to_increment_vec[i] - unprocessed_entries_delta;
        ++next_iteration_num_query_attributes;
      }
    }
  }
  m_query_attribute_idx_vec.resize(next_iteration_num_query_attributes); //no heap operations occur here
}

void SingleCellTileDBIterator::print_ALT(const int field_query_idx, std::ostream& fptr) const {
  assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[field_query_idx];
  auto x = get_buffer_and_index(field_query_idx);
  genomicsdb_columnar_field.print_ALT_data_in_buffer_at_index(fptr, x.first, x.second);
}

void SingleCellTileDBIterator::print_csv(const int field_query_idx, std::ostream& fptr) const {
  assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[field_query_idx];
  auto x = get_buffer_and_index(field_query_idx);
  genomicsdb_columnar_field.print_data_in_buffer_at_index_as_csv(fptr, x.first, x.second);
}

#ifdef DO_PROFILING
void SingleCellTileDBIterator::increment_num_cells_traversed_stats(const uint64_t num_cells_incremented) {
  m_num_cells_traversed_stats[GenomicsDBIteratorStatsEnum::TOTAL_CELLS_TRAVERSED] += num_cells_incremented;
  if (m_in_find_intersecting_intervals_mode)
    m_num_cells_traversed_stats[GenomicsDBIteratorStatsEnum::NUM_CELLS_TRAVERSED_IN_FIND_INTERSECTING_INTERVALS_MODE]
    += num_cells_incremented;
  auto num_useless_cells = (num_cells_incremented - 1ull);
  m_num_cells_traversed_stats[m_in_find_intersecting_intervals_mode
                              ? GenomicsDBIteratorStatsEnum::NUM_USELESS_CELLS_IN_FIND_INTERSECTING_INTERVALS_MODE
                              : GenomicsDBIteratorStatsEnum::NUM_USELESS_CELLS_IN_SIMPLE_TRAVERSAL_MODE] += num_useless_cells;
  auto idx = std::min<uint64_t>(num_useless_cells, m_useless_cell_interval_lengths_histogram.size()-1u);
  ++(m_useless_cell_interval_lengths_histogram[idx]);
}
#endif

#if defined(DO_PROFILING) && defined(PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL)
void SingleCellTileDBIterator::update_sliding_window_to_profile_num_cells_to_traverse(
  const GenomicsDBColumnarField& coords_columnar_field) {
  const auto* curr_coords = reinterpret_cast<const int64_t*>(
                              coords_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
                                coords_columnar_field.get_live_buffer_list_tail(),
                                coords_columnar_field.get_curr_index_in_live_list_tail()
                              )
                            );
  auto curr_row_idx = curr_coords[0];
  auto num_observed_cells_in_curr_window_for_row_idx_before_increment
    = m_row_idx_to_num_observed_cells_in_curr_window[curr_row_idx];
  ++(m_row_idx_to_num_observed_cells_in_curr_window[curr_row_idx]);
  ++m_num_observed_cells_in_curr_window;
  m_observed_cells_in_curr_window.emplace_back(curr_row_idx, curr_coords[1]);
  //no cells for this row were seen in the current window before this cell
  //Check if window is "complete"
  if (num_observed_cells_in_curr_window_for_row_idx_before_increment == 0u) {
    ++m_num_observed_row_idxs_in_curr_window;
    if (m_num_observed_row_idxs_in_curr_window == m_query_config->get_num_rows_to_query())
      std::cerr << m_observed_cells_in_curr_window[0].second << '\t'
                << curr_coords[1] << '\t' << m_num_observed_cells_in_curr_window << '\n';
    //Current sliding window is "complete" - move left edge of window
    //Complete == contains at least 1 cell from every row
    auto window_left_idx = 0ull;
    auto last_column = -1ll;
    //One row has no cells in window OR
    //Cells in the same column
    while (m_num_observed_row_idxs_in_curr_window == m_query_config->get_num_rows_to_query()
           || (window_left_idx < m_observed_cells_in_curr_window.size()
               && last_column == m_observed_cells_in_curr_window[window_left_idx].second)) {
      const auto& left_cell_pair = m_observed_cells_in_curr_window[window_left_idx];
      auto left_row_idx = left_cell_pair.first;
      --(m_row_idx_to_num_observed_cells_in_curr_window[left_row_idx]);
      //No cells for this row left in the window
      if (m_row_idx_to_num_observed_cells_in_curr_window[left_row_idx] == 0u)
        --m_num_observed_row_idxs_in_curr_window;
      ++window_left_idx;
      last_column = left_cell_pair.second;
    }
    //Delete cells from the left of the window
    if (window_left_idx > 0u) {
      m_observed_cells_in_curr_window.erase(m_observed_cells_in_curr_window.begin(),
                                            m_observed_cells_in_curr_window.begin()+window_left_idx);
      assert(window_left_idx <= m_num_observed_cells_in_curr_window);
      m_num_observed_cells_in_curr_window -= window_left_idx;
    }
  }
}
#endif

//GVCF iterator
GenomicsDBGVCFIterator::GenomicsDBGVCFIterator(TileDB_CTX* tiledb_ctx,
    const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
    const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
  : GenomicsDBGVCFIterator(tiledb_ctx,
                             0,
                             vid_mapper, variant_array_schema,
                             array_path, query_config, buffer_size)
{}

#define INIT_INVALID_HISTOGRAM_NUM_BINS 20u

GenomicsDBGVCFIterator::GenomicsDBGVCFIterator(TileDB_CTX* tiledb_ctx,
    const TileDB_Array* tiledb_array,
    const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
    const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
  : SingleCellTileDBIterator(tiledb_ctx, tiledb_array,
      vid_mapper, variant_array_schema, array_path, query_config, buffer_size),
  m_current_start_position(-1ll),
  m_current_end_position(-1ll),
  m_next_start_position(-1ll),
  m_cell(new GenomicsDBGVCFCell(this)),
  m_alleles_combiner(*this, query_config.get_num_rows_to_query()),
  m_gt_remapper(query_config.get_query_idx_for_known_field_enum(GVCF_GT_IDX), *this, m_alleles_combiner)
#ifdef DO_PROFILING
  , m_num_times_initialized(INIT_INVALID_HISTOGRAM_NUM_BINS, 0u),
  m_num_times_invalidated(INIT_INVALID_HISTOGRAM_NUM_BINS, 0u),
  m_bin_size((m_query_config->get_num_rows_to_query()+INIT_INVALID_HISTOGRAM_NUM_BINS-1)/INIT_INVALID_HISTOGRAM_NUM_BINS)
#endif
{
  m_REF_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_REF_IDX);
  m_ALT_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_ALT_IDX);
  //Constructor for SingleCellTileDBIterator invokes SingleCellTileDBIterator::begin_new_query_column_interval()
  //So, it's safe to call GenomicsDBGVCFIterator::begin_new_query_column_interval()
  begin_new_query_column_interval();
}

GenomicsDBGVCFIterator::~GenomicsDBGVCFIterator() {
  if(m_cell)
    delete m_cell;
  m_cell = 0;
#ifdef DO_PROFILING
  std::cerr << "INIT INVALID HISTOGRAM\n";
  for(auto i=0u;i<INIT_INVALID_HISTOGRAM_NUM_BINS;++i)
    std::cerr << i*m_bin_size << " "<<m_num_times_initialized[i] << " " <<m_num_times_invalidated[i] <<"\n";
#endif
}

template<>
void GenomicsDBGVCFIterator::update_num_deletions_and_MNVs<true>(const size_t row_query_idx) {
  const auto REF_ptr_length_pair = get_raw_pointer_and_length_for_query_idx(row_query_idx, m_REF_query_idx);
  const auto ALT_ptr_length_pair = get_raw_pointer_and_length_for_query_idx(row_query_idx, m_ALT_query_idx);
  STRING_VIEW REF(reinterpret_cast<const char*>(REF_ptr_length_pair.first), REF_ptr_length_pair.second);
  STRING_VIEW delimited_ALT(reinterpret_cast<const char*>(ALT_ptr_length_pair.first), ALT_ptr_length_pair.second);
  m_alleles_combiner.insert_allele_info(row_query_idx, REF, delimited_ALT,
      (m_live_cell_markers.get_begin(row_query_idx) < m_current_start_position));
}

template<>
void GenomicsDBGVCFIterator::update_num_deletions_and_MNVs<false>(const size_t row_query_idx) {
  m_alleles_combiner.remove_allele_info(row_query_idx);
}

void GenomicsDBGVCFIterator::begin_new_query_column_interval() {
  //Assumes SingleCellTileDBIterator::begin_new_query_column_interval is already invoked
  //Hence, the two loops execute only when SingleCellTileDBIterator::begin_new_query_column_interval
  //has determined that there is valid data to work on
  if(SingleCellTileDBIterator::end())
    return;
  m_alleles_combiner.reset_for_next_query_interval();
  const auto& coords_columnar_field = m_fields.back();
  const auto& END_columnar_field = m_fields[m_END_query_idx];
  m_current_start_position = (m_query_config->get_num_column_intervals() == 0u)
    ? -1ll : m_query_config->get_column_begin(m_query_column_interval_idx);
  m_query_interval_limit = (m_query_config->get_num_column_intervals() == 0u)
    ? m_variant_array_schema->dim_domains()[1].second
    : m_query_config->get_column_end(m_query_column_interval_idx);
  //m_current_start_position might have to be updated if
  //(a) there are no intersecting intervals for query interval begin AND
  //(b) there's valid data in simple traversal mode
  if(m_PQ_live_cell_markers.empty() && !m_done_reading_from_TileDB)
    m_current_start_position = (reinterpret_cast<const int64_t*>(
	coords_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
	  coords_columnar_field.get_live_buffer_list_tail(),
	  coords_columnar_field.get_curr_index_in_live_list_tail()
	  )
	))[1];
  //All elements of m_PQ_live_cell_markers begin before the query interval
  //and intersect the query interval - hence, the m_current_start_position
  //of all elements of PQ is simply the start of the query interval
  while(!m_PQ_live_cell_markers.empty()) {
    auto row_query_idx = m_PQ_live_cell_markers.top().second;
    m_PQ_live_cell_markers.pop();
    update_num_deletions_and_MNVs<true>(row_query_idx);
    m_end_set.insert(GVCFEndSetElementTy(m_live_cell_markers.get_end(row_query_idx), row_query_idx));
  }
  //There might be more elements which start at query interval begin
  //which are accessed in simple traversal mode
  m_alleles_combiner.reset_before_adding_new_sample_info_at_current_position();
  fill_end_set_in_simple_traversal_mode();
  m_alleles_combiner.finished_updating_allele_info_for_current_position();
}

const GenomicsDBGVCFIterator& GenomicsDBGVCFIterator::operator++() {
  assert(!m_end_set.empty());
  m_current_start_position = m_current_end_position+1;
  //Past the limit, clean-up end set and move to next interval
  if(m_current_start_position > m_query_interval_limit)
    m_current_start_position = INT64_MAX;
  auto iter = m_end_set.begin();
#ifdef DO_PROFILING
  auto num_times_invalidated = 0ull;
#endif
  //Remove all entries that finish at or before m_current_end_position
  for(;iter!=m_end_set.end();++iter) {
    auto min_end_element = *iter;
    const auto END_field_value = min_end_element.first;
    const auto row_query_idx = min_end_element.second;
    if(END_field_value < m_current_start_position) {
      update_num_deletions_and_MNVs<false>(row_query_idx);
      decrement_num_live_entries_in_live_cell_marker(row_query_idx);
#ifdef DO_PROFILING
      ++num_times_invalidated;
#endif
    }
    else
      break;
  }
#ifdef DO_PROFILING
  ++(m_num_times_invalidated[num_times_invalidated/m_bin_size]);
#endif
  if(iter == m_end_set.end()) { //all entries have END <= m_current_end_position
    m_end_set.clear();
    m_current_start_position = m_next_start_position;
  }
  else
    m_end_set.erase(m_end_set.begin(), iter);
  m_alleles_combiner.reset_before_adding_new_sample_info_at_current_position();
  if(m_current_start_position == m_next_start_position)
    fill_end_set_in_simple_traversal_mode();
  else {
    assert(!m_end_set.empty() || m_current_start_position == INT64_MAX);
    update_current_end_position();
  }
  //No more valid data - move to next query interval if needed
  if(m_end_set.empty() && m_done_reading_from_TileDB) {
    SingleCellTileDBIterator::begin_new_query_column_interval();
    GenomicsDBGVCFIterator::begin_new_query_column_interval();
  }
  m_alleles_combiner.finished_updating_allele_info_for_current_position();
  m_at_new_query_column_interval = false;
  return *this;
}

void GenomicsDBGVCFIterator::fill_end_set_in_simple_traversal_mode() {
  //We keep searching till we find a cell whose start position is > m_current_start_position
  //or no more cells exist for current query interval
  m_next_start_position = -1;
  const auto& coords_columnar_field = m_fields.back();
  const auto& END_columnar_field = m_fields[m_END_query_idx];
#ifdef DO_PROFILING
  auto num_times_initialized = 0ull;
#endif
  while(!m_done_reading_from_TileDB && m_next_start_position == -1) {
    const auto* coords = reinterpret_cast<const int64_t*>(
	coords_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
	  coords_columnar_field.get_live_buffer_list_tail(),
	  coords_columnar_field.get_curr_index_in_live_list_tail()
	  )
	);
    const auto END_field_value = *(reinterpret_cast<const int64_t*>(
	  END_columnar_field.get_raw_pointer_to_data_in_buffer_at_index(
	    END_columnar_field.get_live_buffer_list_tail(),
	    END_columnar_field.get_curr_index_in_live_list_tail()
	    )
	  ));
    const auto row_query_idx = m_query_config->get_query_row_idx_for_array_row_idx(coords[0]);
    if(coords[1] == m_current_start_position) {
      assert(m_live_cell_markers.get_row_idx(row_query_idx) == coords[0]);
      //Handle overlapping intervals - this row already has an entry in the set
      if(m_live_cell_markers.is_valid(row_query_idx)) {
	//Remove from m_end_set
	auto iter = m_end_set.find(GVCFEndSetElementTy(m_live_cell_markers.get_end(row_query_idx), row_query_idx));
	assert(iter != m_end_set.end());
	m_end_set.erase(iter);
        update_num_deletions_and_MNVs<false>(row_query_idx);
	//Decrement num live entries
	decrement_num_live_entries_in_live_cell_marker(row_query_idx);
      }
      //Add new cell, initialize live marker from tail of buffer and increment live entries
      initialize_live_cell_marker_from_tail(row_query_idx, coords, END_field_value);
      update_num_deletions_and_MNVs<true>(row_query_idx);
      m_end_set.insert(GVCFEndSetElementTy(END_field_value, row_query_idx));
      advance_to_next_useful_cell(1u);
#ifdef DO_PROFILING
      ++num_times_initialized;
#endif
    }
    else
      m_next_start_position = coords[1];
  }
#ifdef DO_PROFILING
  ++(m_num_times_initialized[num_times_initialized/m_bin_size]);
#endif
  update_current_end_position();
}

void GenomicsDBGVCFIterator::update_current_end_position() {
  if(m_next_start_position < 0) {
    if(m_query_config->get_num_column_intervals() == 0u) {
      auto& dim_domains = m_variant_array_schema->dim_domains();
      m_next_start_position = dim_domains[1].second;
    }
    else
      m_next_start_position = m_query_config->get_column_end(m_query_column_interval_idx);
    //to avoid wrap-around
    m_next_start_position = std::min(m_next_start_position, INT64_MAX-1)+1;
  }
  if(m_alleles_combiner.get_num_calls_with_deletions_or_MNVs() > 0u)
    m_current_end_position = m_current_start_position; //to put spanning deletions
  else {
    const auto min_END_field_value = m_end_set.empty() ? INT64_MAX
      : (*(m_end_set.begin())).first;
    m_current_end_position = std::min(min_END_field_value, m_next_start_position-1);
  }
}
