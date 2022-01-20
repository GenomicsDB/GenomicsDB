#include "tiledb_loader_rc.h"
#include "tiledb_utils.h"
#include <queue>
#include <functional>
#include <memory>
#include <regex>

SinglePosition2TileDBLoader::SinglePosition2TileDBLoader(const std::string& config_filename, const int idx, std::string gtf_name, std::string gi_name)
  : GenomicsDBImportConfig(), m_idx(idx) {

  GenomicsDBConfigBase::read_from_JSON(GenomicsDBConfigBase::read_from_file(config_filename, idx), idx);
  //m_vid_mapper = GenomicsDBConfigBase::get_vid_mapper();

  auto array_name = get_array_name(m_idx);
  auto workspace = get_workspace(m_idx);

  auto dim_names = std::vector<std::string>({"samples", "position"});
  auto dim_domains = std::vector<std::pair<int64_t, int64_t>>({ {0, INT64_MAX-1}, {0, INT64_MAX-1}});
  std::vector<std::string> attribute_names;
  std::vector<std::type_index> types;
  std::vector<int> num_vals;
  
  /*for(int i = 0; i < m_vid_mapper.get_num_fields(); i++) {
    auto& inf = m_vid_mapper.get_field_info(i);
    std::cout << inf.m_name << std::endl;

    attribute_names.push_back(inf.m_name);
    types.push_back(inf.get_tiledb_type().get_tuple_element_type_index(0u));
    num_vals.push_back(inf.m_length_descriptor.is_fixed_length_field()
                       ? inf.m_length_descriptor.get_num_elements()
                       : TILEDB_VAR_NUM);
  }*/

  attribute_names.push_back("START");
  types.push_back(typeid(int64_t));
  num_vals.push_back(1);

  attribute_names.push_back("END");
  types.push_back(typeid(int64_t));
  num_vals.push_back(1);

  attribute_names.push_back("SCORE");
  types.push_back(typeid(float));
  num_vals.push_back(1);

  attribute_names.push_back("NAME");
  types.push_back(typeid(char));
  num_vals.push_back(TILEDB_VAR_NUM);

  attribute_names.push_back("GENE");
  types.push_back(typeid(char));
  num_vals.push_back(TILEDB_VAR_NUM);

  // Add type for coords
  types.push_back(std::type_index(typeid(int64_t)));

  //std::vector<int> compression_types(types.size(), TILEDB_NO_COMPRESSION);
  //std::vector<int> compression_levels(types.size(), 0);
  
  std::vector<int> compression_types(types.size(), TILEDB_GZIP);
  std::vector<int> compression_levels(types.size(), 1);

  // m_import_config_ptr->get_tiledb_compression_type()

  //m_array_schema = VariantArraySchema(array_name, attribute_names, dim_names, dim_domains, types, num_vals, compression_types, compression_levels);
  //auto m_array_schema_ptr = &m_array_schema;
  //m_vid_mapper.build_tiledb_array_schema(m_array_schema_ptr, array_name, false, TILEDB_NO_COMPRESSION, 0, true);
  m_array_schema = VariantArraySchema(array_name, attribute_names, dim_names, dim_domains, types, num_vals, compression_types, compression_levels);

//=================================================================

  m_storage_manager = std::make_shared<VariantStorageManager>(workspace, m_segment_size, enable_shared_posixfs_optimizations());
  if (delete_and_create_tiledb_array()) {
    m_storage_manager->delete_array(array_name);
  }

  //Open array in write mode
  m_array_descriptor = m_storage_manager->open_array(array_name, &m_vid_mapper, "w");

  //Check if array already exists
  //Array does not exist - define it first
  if (m_array_descriptor < 0) {
    if (m_storage_manager->define_array(&m_array_schema,
                                        m_num_cells_per_tile,
                                        disable_delta_encode_offsets(),
                                        disable_delta_encode_coords(),
                                        enable_bit_shuffle_gt(),
                                        enable_lz4_compression_gt()) != TILEDB_OK) {
      throw LoadOperatorException(std::string("Could not define TileDB array") + "\nTileDB error message : " + tiledb_errmsg);
    }
    //Open array in write mode
    m_array_descriptor = m_storage_manager->open_array(array_name, &m_vid_mapper, "w");
  } else if (fail_if_updating()) {
    throw LoadOperatorException(std::string("Array ") + workspace + "/" + array_name + " exists and flag \"fail_if_updating\" is set to true in the loader JSON configuration");
  }
  if (m_array_descriptor < 0) {
    throw LoadOperatorException(std::string("Could not open TileDB array for loading") + "\nTileDB error message : " + tiledb_errmsg);
  }
  //m_storage_manager->update_row_bounds_in_array(m_array_descriptor, get_row_bounds().first, std::min(get_row_bounds().second, m_vid_mapper.get_max_callset_row_idx()));
  m_storage_manager->update_row_bounds_in_array(m_array_descriptor, get_row_bounds().first, std::min(get_row_bounds().second, m_vid_mapper.get_max_callset_row_idx()));
  m_storage_manager->write_column_bounds_to_array(m_array_descriptor, get_column_partition(m_idx).first, get_column_partition(m_idx).second);


//=================================================================

  // use gi file if exists, otherwise gtf/gff
  std::ifstream gi_file(gi_name);
  std::ifstream gtf_file(gtf_name);
  std::string type;

  if(gi_file.good()) {
    deserialize_transcript_map(gi_file);
  }
  else if(gtf_file.good()){
    std::regex reg_gtf("(.*)(gtf)($)");
    if(std::regex_match(gtf_name, reg_gtf)) {
      type = "gtf";

      read_uncompressed_gtf(gtf_file, type);
      gi_file.close();
      std::ofstream gi_out(gi_name);
      serialize_transcript_map(gi_out);
    }
    else {
      std::cout << "\t\t\t\tREGEX DOES NOT MATCH" << std::endl;
      logger.error("No valid GTF/GFF or GI file specified, cannot import transcriptomics style data");
    }
  }
  else {
    logger.error("No valid GTF/GFF or GI file specified, cannot import transcriptomics style data");
  }

  /*size_t then = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  read_uncompressed_gtf(gtf_file, type);
  size_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  double time = double(now - then) / 1000;
  std::cout << "Read uncompressed gtf took " << time << " s" << std::endl;

  then = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  std::ofstream serial_file;
  serial_file.open("/nfs/home/andrei/benchmarking_requirements/serialized", std::ios::out | std::ios::binary);
  serialize_transcript_map(serial_file);
  now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  time = double(now - then) / 1000;
  std::cout << "serialize took " << time << " s" << std::endl;

  serial_file.close();
  transcript_map.clear();

  then = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  std::ifstream deserial_file;
  deserial_file.open("/nfs/home/andrei/benchmarking_requirements/serialized", std::ios::in | std::ios::binary);
  deserialize_transcript_map(deserial_file);
  now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  time = double(now - then) / 1000;
  std::cout << "deserialize took " << time << " s" << std::endl;*/
}

// TODO: maybe also create end cell, add ability to specify here
// assume valid input
bool SinglePosition2TileDBLoader::info_to_cell(int64_t start, int64_t end, std::string name, std::string gene, float score, uint64_t row) {
  //construct cell
  std::vector<uint8_t> cell(16);
  *(reinterpret_cast<uint64_t*>(cell.data())) = row; // write row in cell

  int64_t pos = start;

  if(start > end) { // indicates end cell
    start = end;
    end = pos;
  }

  *(reinterpret_cast<uint64_t*>(cell.data()) + 1) = pos; // write position in cell

  // reserve space for cell size
  for(int j = 0; j < sizeof(size_t); j++) {
    cell.push_back(0);
  }

  // attributes
  for(int j = 0; j < m_array_schema.attribute_num(); j++) {
    std::string attribute_name = m_array_schema.attribute_name(j);
    if(attribute_name == "START") {
      cell.insert(cell.end(), {0, 0, 0, 0, 0, 0, 0, 0});
      *(reinterpret_cast<int64_t*>(cell.data() + cell.size() - 8)) = start;
    }
    if(attribute_name == "END") {
      cell.insert(cell.end(), {0, 0, 0, 0, 0, 0, 0, 0});
      *(reinterpret_cast<int64_t*>(cell.data() + cell.size() - 8)) = end;
    }
    if(attribute_name == "SCORE") {
      int s;
      memcpy(&s, &score, 4);

      cell.insert(cell.end(), {0, 0, 0, 0});
      *(reinterpret_cast<float*>(cell.data() + cell.size() - 4)) = score;
    }
    if(attribute_name == "NAME") {
      cell.insert(cell.end(), {0, 0, 0, 0});
      *(reinterpret_cast<uint32_t*>(cell.data() + cell.size() - 4)) = name.length();
      for(auto c : name) {
        cell.push_back(c);
      }
    }
    if(attribute_name == "GENE") {
      cell.insert(cell.end(), {0, 0, 0, 0});
      *(reinterpret_cast<uint32_t*>(cell.data() + cell.size() - 4)) = gene.length();
      for(auto c : gene) {
        cell.push_back(c);
      }
    }
  }
  // fill in cell size
  *(reinterpret_cast<size_t*>(cell.data() + 2*sizeof(int64_t))) = cell.size();
  //write
  m_storage_manager->write_cell_sorted(m_array_descriptor, reinterpret_cast<const void*>(cell.data()));
  return true;
}

// chrom is empty string if invalid
std::tuple<std::string, int64_t, int64_t, std::string, float> SinglePosition2TileDBLoader::parse_and_check(std::string str, VidMapper& vid_mapper) {
  std::string chrom, name;
  int64_t start, end;
  float score;

  std::vector<std::string> vec;
  size_t index;
  while((index = str.find("\t")) != std::string::npos) {
    vec.push_back(str.substr(0, index));
    str.erase(0, index + 1);
  }
  vec.push_back(str);

  if(vec.size() >= 5) {
    chrom = vec[0];
    name = vec[3];
    try {
      bool valid = vid_mapper.get_tiledb_position(start, chrom, std::stol(vec[1]));
      valid = valid && vid_mapper.get_tiledb_position(end, chrom, std::stol(vec[2]));
      score = std::stof(vec[4]);
      if(!valid) {
        return {"", 0, 0, "", 0.};
      }
    }
    catch(...) {
      return {"", 0, 0, "", 0.};
    }
  }
  else {
    return {"", 0, 0, "", 0.};
  }

  return {chrom, start, end, name, score};
}

/*cell_info SinglePosition2TileDBLoader::next_cell_info(std::ifstream& file, int type, int ind) {
  if(!type) {
    std::string str;
    while(std::getline(file, str)) {
      std::cout << "============================================= " << str << std::endl;

      auto[chrom, start, end, name, score] = parse_and_check(str, m_vid_mapper);

      if(chrom == "") {
        continue;
      }
      else {
        return {start, end, name, "", score, ind};
      }
    }
    return {-1, -1, "", "", 0, -1};
  }
  else {
    int mat_fields = 6;
    int start_ind = 0, end_ind = 1, name_ind = 2, sample_ind = 4, score_ind = 6;

    std::vector<std::string> fields;
    std::string temp;

    for(int i = 0; i < mat_fields; i++) {
      if(file >> temp) {
        fields.push_back(temp);
      }
      else {
        return {-1, -1, "", "", 0, -1};
      }
    }
    try {
      cell_info inf = {std::stol(fields[start_ind]), std::stol(fields[end_ind]), fields[name_ind], "placeholder", std::stof(fields[score_ind]), std::stoi(fields[sample_ind])};
      return inf;
    }
    catch (...) {
      return {-1, -1, "", "", 0, -1};
    }
  }
}*/

// get line using TileDBUtils api to work with cloud storage as well
// return value indicates if line was read
bool TranscriptomicsFileReader::generalized_getline(std::string& retval) {
  retval = "";

  while(m_chars_read < m_file_size || m_str_buffer.size()) {
    int idx = m_str_buffer.find('\n');
    if(idx != std::string::npos) {
      retval = retval + m_str_buffer.substr(0, idx); // exclude newline
      m_str_buffer.erase(0, idx + 1); // erase newline
      return true;
    }

    retval = retval + m_str_buffer;
    m_str_buffer.clear();

    int chars_to_read = std::min<ssize_t>(m_buffer_size, m_file_size - m_chars_read);

    if(chars_to_read) {
      TileDBUtils::read_file(m_fname, m_chars_read, m_buffer, chars_to_read);
       m_chars_read += chars_to_read;
    }

    m_str_buffer.insert(m_str_buffer.end(), m_buffer, m_buffer + chars_to_read);
  }

  return false;
}

cell_info BedReader::next_cell_info() {
  std::string str;
  while(generalized_getline(str)) {

    auto[chrom, start, end, name, score] = SinglePosition2TileDBLoader::parse_and_check(str, m_vid_mapper);

    if(chrom == "") {
      continue;
    }
    else {
      return {start, end, name, "NA", score, m_sample_idx, m_file_idx};
    }
  }
  return {-1, -1, "", "", 0, -1, m_file_idx}; 
}

// assuming format is gene, name, scores
cell_info MatrixReader::next_cell_info() {
  int mat_fields = 6;
  int start_ind = 0, end_ind = 1, name_ind = 2, sample_ind = 4, score_ind = 6;

  while(m_current_sample >= m_current_line.size()) { // read until next valid line or eof
    std::string line;

    // eof
    if(!generalized_getline(line)) {
      return {-1, -1, "", "", 0, -1, m_file_idx};
    }

    std::stringstream ss(line);

    ss >> m_current_gene >> m_current_name;

    if(m_transcript_map.count(m_current_gene)) {
      auto[s, e] = m_transcript_map[m_current_gene];
      m_current_start = s;
      m_current_end = e;
    }
    else {
      continue; // gene not in map, skip row
    }

    std::string str;
    while(ss >> str) {
      m_current_line.push_back(str);
    }
    m_current_sample = 0;
    idx_to_row_pos = 0;

    // done setting up current row in matrix
    
    // check if done with row in matrix
    if(idx_to_row_pos >= idx_to_row.size()) {
      continue;
    }

    // otherwise go to next relevent column
    try {
      auto p = idx_to_row[idx_to_row_pos++];
      if(p.first >= samples.size()) {
        logger.error("Error index in file {} is outside of matrix file's {} columns", p.first, samples.size());
        exit(1);
      }

      std::string str = m_current_line[p.first];

      float score = std::stof(str);
      cell_info inf = {m_current_start, m_current_end, m_current_name, m_current_gene, score, p.second, m_file_idx};
      return inf;
    }
    catch (...) {
      return {-1, -1, "", "", 0, -1, m_file_idx};
    }
  }
  return {-1, -1, "", "", 0, -1, m_file_idx};
}

/*void SinglePosition2TileDBLoader::read_compressed_gtf(std::string fname) {
  size_t filesize, buffer_size, allocated_buffer_size, chunk_size;
  void* buffer;
  TileDBUtils::gzip_read_buffer(fname, filesize, buffer, buffer_size, allocated_buffer_size, chunk_size);
}*/

void SinglePosition2TileDBLoader::read_uncompressed_gtf(std::istream& input, std::string format) {
  if(format != "gff" and format != "gtf") {
    logger.error("unrecognized file type {}", format);
    return;
  }
  
  std::string str;

  int ind = -1;
  while(std::getline(input, str)) {
    if(str[0] == '#') {
      continue;
    }

    ++ind;

    std::stringstream ss(str);
    std::string field;
    std::vector<std::string> fields;

    // read the 8 fields before attribute
    for(int i = 0; i < 8; i++) {
      ss >> field;
      fields.push_back(field);
    }

    long start, end;

    if(fields[2] != "transcript") {
      continue;
    }

    try {
      start = std::stol(fields[3]);
      end = std::stol(fields[4]);
    }
    catch (...) {
      continue;
    }

    std::string attributes;
    std::getline(ss, attributes);

    int lo = attributes.find("transcript_id ") + 14; // character after end of pattern
    int hi = attributes.find(";", lo);

    std::string tid = attributes.substr(lo, hi - lo);
    tid.erase(std::remove(tid.begin(), tid.end(), '\"'), tid.end()); // remove quotes

    auto sz = transcript_map.size();
    transcript_map[tid] = {start, end};
 
    if(sz == transcript_map.size()) {
      std::cout << tid << " is duplicate, start/end: " << start << ", " << end << std::endl;
    }
  }
}

void SinglePosition2TileDBLoader::serialize_transcript_map(std::ostream& output) {
  char version = 0;
  output.write(&version, 1);

  // write 5B number of entries
  int64_t size = transcript_map.size();
  output.write((char*)&size, 5);

  // write 2B string length, string, 5B start, 5B end
  for(auto& a : transcript_map) {
    std::string name = a.first;
    int16_t len = name.length();
    int64_t start = a.second.first;
    int64_t end = a.second.second;

    output.write((char*)&len, 2);
    output.write(name.c_str(), len);
    output.write((char*)&start, 5);
    output.write((char*)&end, 5);
  }
}

void SinglePosition2TileDBLoader::deserialize_transcript_map(std::istream& input) {
  char version;
  input.read(&version, 1);
  if(version) {
    std::cout << "version " << (int)version << " not supported" << std::endl;
  }

  int64_t size = 0;
  input.read(((char*)&size), 5);

  for(int i = 0; i < size; i++) {
    int16_t len;
    input.read((char*)&len, 2);
    char name[len];
    input.read(name, len);
    int64_t start = 0;
    int64_t end = 0;
    input.read((char*)&start, 5);
    input.read((char*)&end, 5);
    transcript_map[std::string(name, len)] = {start, end};
  }
}

// TODO sample name
void SinglePosition2TileDBLoader::read_all() {
  // construct transcript map

  auto remove_file = [](cell_info c) { std::get<6>(c) = -1; return c; };
  auto reverse_info = [](cell_info c) {
    int64_t start = std::get<0>(c);
    std::get<0>(c) = std::get<1>(c);
    std::get<1>(c) = start;
    return c;
  };

  auto workspace = get_workspace(m_idx);
  auto array_name = get_array_name(m_idx);

  auto comp = [](cell_info l, cell_info r) { return (std::get<0>(l) > std::get<0>(r)) || (std::get<0>(l) == std::get<0>(r) && std::get<5>(l) > std::get<5>(r)); };
  std::priority_queue<cell_info, std::vector<cell_info>, decltype(comp)> pq(comp);

  //std::vector<std::unique_ptr<std::ifstream>> files;
  std::vector<std::shared_ptr<TranscriptomicsFileReader>> files;
  // 0 for bed, 1 for matrix
  std::vector<int> file_types;

  for(int i = 0; i < m_vid_mapper.get_num_callsets(); i++) {
    auto ci = m_vid_mapper.get_callset_info(i);
    int f_ind = ci.m_file_idx;
    std::string fname = m_vid_mapper.get_file_info(f_ind).m_name;
    int64_t idx_in_file = ci.m_idx_in_file;
    int64_t row_idx = ci.m_row_idx;

    // keep track of which files have already appeared
    std::map<std::string, std::shared_ptr<TranscriptomicsFileReader>> previous_files;
    
    //auto fi = m_vid_mapper.get_file_info(i);

    int type;

    if(std::regex_match(fname, std::regex("(.*)(bed)($)"))) { // bed file
      if(previous_files.count(fname)) {
        logger.error("File {} appears twice in callset mapping file", fname);
        continue;
      }

      type = 0;
      file_types.push_back(0);
      files.push_back(std::make_shared<BedReader>(fname, files.size(), m_vid_mapper, row_idx));
      
      previous_files.insert({fname, files.back()});
    }
    else if(std::regex_match(fname, std::regex("(.*)(resort)($)"))) { // matrix
      if(!previous_files.count(fname)) { // create reader object for file
        type = 1;
        file_types.push_back(1);
        files.push_back(std::make_shared<MatrixReader>(fname, files.size(), m_vid_mapper, transcript_map));
        previous_files.insert({fname, files.back()});
      }

      MatrixReader& file = dynamic_cast<MatrixReader&>(*(previous_files[fname]));
      std::pair<int, int64_t> p = {idx_in_file, row_idx};

      if(file.idx_to_row.size()) { // check that callset is sorted
        auto[bi, br] = file.idx_to_row.back();

        if(bi >= idx_in_file) {
          logger.error("callsets for {} are not sorted by index in file or have duplicates", fname);
          exit(1);
        }

        if(br >= row_idx) {
          logger.error("callsets for {} are not sorted by row index or have duplicates", fname);
          exit(1);
        }
      }

      file.idx_to_row.push_back(p); // insert last read callset into index to row map in applicable matrix file reader
    }
    else {
      file_types.push_back(2);
      logger.error("Transcriptomics: Unknown file type for file {}", fname);
      exit(1);
    }


    // put first cell from each file in pq
    for(auto& file_ptr : files) {
      auto inf = file_ptr->next_cell_info();
      if(std::get<0>(inf) >= 0) {
        pq.push(remove_file(inf)); // set file index to -1 to indicate this is a start (and not to read the file for another cell yet)
        pq.push(reverse_info(inf)); // reverse start and end to indicate this is an end
      }
    }
  }

  while(pq.size()) {
    auto[start, end, name, gene, score, ind, file_ind] = pq.top();
    pq.pop();

    info_to_cell(start, end, name, gene, score, ind);

    if(file_ind >= 0) { // if cell is an end
      auto tup = files[ind]->next_cell_info(); // read next cell from originiating file
      if(std::get<0>(tup) >= 0 && std::get<6>(tup) != -1) { // check cell info is valid
        pq.push(tup);
      }
    }
  }

  if (m_storage_manager && m_array_descriptor >= 0) {
    m_storage_manager->close_array(m_array_descriptor, consolidate_tiledb_array_after_load());
  }

  read_array();
}

void check_rc(int rc) {
  if (rc) {
    printf("%s", &tiledb_errmsg[0]);
    printf("[Examples::%s] Runtime Error.\n", __FILE__);
    return;
  }
}

void SinglePosition2TileDBLoader::read_array() {
  auto array_name = get_array_name(m_idx);
  auto workspace = get_workspace(m_idx);

  // Initialize context with the default configuration parameters
  TileDB_CTX* tiledb_ctx;
  tiledb_ctx_init(&tiledb_ctx, NULL);

  char a1[200], a2[200], a3[200], a4[200], a5[200], a6[200], a7[200], a8[200], a9[200];
  void* buffers[] = {a1, a2, a3, a4, a5, a6, a7, a8, a9};
  size_t sizes[] = {sizeof(a1), sizeof(a2), sizeof(a3), sizeof(a4), sizeof(a5), sizeof(a6), sizeof(a7), sizeof(a8), sizeof(a9)};

  std::cout << "before initialize itor " << std::endl; 

  int64_t subarray[] = { 0, 0, 100, 349251121 };

  // Initialize array
  TileDB_ArrayIterator* tiledb_it;
  check_rc(tiledb_array_iterator_init(
      tiledb_ctx,                                       // Context
      &tiledb_it,                                       // Array object
      (workspace + "/" + array_name).c_str(),           // Array name
      TILEDB_ARRAY_READ,                                // Mode
      NULL,                                             // Whole domain
      NULL,                                             // All attributes
      0,                                                // Number of attributes
      buffers,                                          // buffers used internally
      sizes));                                          // buffer sizes

  std::cout << "after initialize itor" << std::endl;

  int* val;
  size_t size;
  while(!tiledb_array_iterator_end(tiledb_it)) {
    // START
    int64_t* start = 0;
    size_t start_size;
    tiledb_array_iterator_get_value(
          tiledb_it,             // Array iterator
          0,                     // Attribute id
          (const void**) &start, // Value
          &start_size);          // Value size (useful in variable-sized attributes)
    std::cout << "read start as " << *start << std::endl;
    std::cout << "start size is " << start_size << std::endl;
    // END
    int64_t* end = 0;
    size_t end_size;
    tiledb_array_iterator_get_value(
          tiledb_it,             // Array iterator
          1,                     // Attribute id
          (const void**) &end,   // Value
          &end_size);            // Value size (useful in variable-sized attributes)
    std::cout << "read end as " << *end << std::endl;
    std::cout << "end size is " << end_size << std::endl;
    // SCORE
    float* score = 0;
    size_t score_size;
    tiledb_array_iterator_get_value(
          tiledb_it,             // Array iterator
          2,                     // Attribute id
          (const void**) &score, // Value
          &score_size);          // Value size (useful in variable-sized attributes)
    std::cout << "read score as " << *score << std::endl;
    std::cout << "score size is " << score_size << std::endl;
    // NAME
    char* name = 0;
    size_t name_size;
    tiledb_array_iterator_get_value(
          tiledb_it,             // Array iterator
          3,                     // Attribute id
          (const void**) &name,  // Value
          &name_size);           // Value size (useful in variable-sized attributes)
    std::cout << "read name as " << std::endl;
    for(int j = 0; j < name_size; j++) {
      std::cout << name[j];
    }
    std::cout << std::endl << "name size is " << name_size << std::endl;
    // GENE
    char* gene = 0;
    size_t gene_size;
    tiledb_array_iterator_get_value(
          tiledb_it,             // Array iterator
          4,                     // Attribute id
          (const void**) &gene,  // Value
          &gene_size);           // Value size (useful in variable-sized attributes)
    std::cout << "read gene as " << std::endl;
    for(int j = 0; j < gene_size; j++) {
      std::cout << gene[j];
    }
    std::cout << std::endl << "gene size is " << gene_size << std::endl;

    // COORDS
    int64_t* coords = 0;
    size_t coords_size;
    tiledb_array_iterator_get_value(
          tiledb_it,              // Array iterator
          5,                      // Attribute id
          (const void**) &coords, // Value
          &coords_size);          // Value size (useful in variable-sized attributes)
    std::cout << "read coords as " << coords[0] << ", " << coords[1] << std::endl;
    std::cout << "coords size is " << coords_size << std::endl;

    std::cout << std::endl << std::endl;

    // Advance iterator
    tiledb_array_iterator_next(tiledb_it);
  }
 
  // Finalize array 
  tiledb_array_iterator_finalize(tiledb_it);

  // Finalize context
  tiledb_ctx_finalize(tiledb_ctx);

  /*
  // Prepare cell buffers
  int64_t buffer_start[4];
  int64_t buffer_end[4];
  float buffer_score[4];
  size_t buffer_name[4];
  char buffer_var_name[20];
  int64_t buffer_coords[8];
  void* buffers[] =
      { buffer_start, buffer_end, buffer_score, buffer_name, buffer_var_name, buffer_coords };
  size_t buffer_sizes[] =
  {
      sizeof(buffer_start),
      sizeof(buffer_end),
      sizeof(buffer_score),
      sizeof(buffer_name),
      sizeof(buffer_var_name),
      sizeof(buffer_coords)
  };

  std::cout << "read array buffer sizes" << std::endl;
  for(auto o : buffer_sizes) {
    std::cout << o << std::endl;
  }

  // Read from array
  check_rc(tiledb_array_read(tiledb_array, buffers, buffer_sizes));

  std::cout << "after read " << std::endl; */

  // Print cell values
  /*int64_t result_num = buffer_sizes[0] / sizeof(int);
  printf("coords\t a1\t   a2\t     (a3.first, a3.second)\n");
  printf("--------------------------------------------------\n");
  for(int i=0; i<result_num; ++i) {
    printf("(%" PRId64 ", %" PRId64 ")", buffer_coords[2*i], buffer_coords[2*i+1]);
    printf("\t %3d", buffer_a1[i]);
    size_t var_size = (i != result_num-1) ? buffer_a2[i+1] - buffer_a2[i]
                                          : buffer_sizes[2] - buffer_a2[i];
    printf("\t %4.*s", int(var_size), &buffer_var_a2[buffer_a2[i]]);
    printf("\t\t (%5.1f, %5.1f)\n", buffer_a3[2*i], buffer_a3[2*i+1]);
  }*/

  //int results = buffer_sizes[0] / sizeof(size_t);

  //for(int i = 0; i < results; i++) {
  //  std::cout << buffer_coords[2*i] << ", " << buffer_coords[2*i + 1] << std::endl;
  //}


  /*for(int i = 0; i < 4; i++) {
    std::cout << i << std::endl;
    std::cout << "start " << buffer_start[i] << std::endl;
    std::cout << "end " << buffer_end[i] << std::endl;
    int s;
    memcpy(&s, buffer_score + i, 4);
    std::cout << "score " << buffer_score[i] << " as int again: " << s << std::endl;
    std::cout << "name " << buffer_name << std::endl;
    std::cout << "coords " << buffer_coords[2*i] << ", " << buffer_coords[2*i + 1] << std::endl;
  }

  for(int i = 0; i < sizeof(buffer_var_name); i++) {
    std::cout << buffer_var_name[i] << std::endl;
  }

  // Finalize the array
  check_rc(tiledb_array_finalize(tiledb_array));

  // Finalize context
  check_rc(tiledb_ctx_finalize(tiledb_ctx));

  return;*/
}
