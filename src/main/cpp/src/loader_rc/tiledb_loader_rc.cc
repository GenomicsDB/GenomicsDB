#include "tiledb_loader_rc.h"

SinglePosition2TileDBLoader::SinglePosition2TileDBLoader(const std::string& config_filename, const int idx)
  : GenomicsDBImportConfig(), m_idx(idx) {

  std::cout << "before constructing GenomicsDBConfigBase" << std::endl;

  GenomicsDBConfigBase::read_from_JSON(GenomicsDBConfigBase::read_from_file(config_filename, idx), idx);
  //m_vid_mapper = GenomicsDBConfigBase::get_vid_mapper();

  std::cout << "Constructor" << std::endl;

  auto array_name = get_array_name(m_idx);
  auto workspace = get_workspace(m_idx);

  auto dim_names = std::vector<std::string>({"samples", "position"});
  auto dim_domains = std::vector<std::pair<int64_t, int64_t>>({ {0, INT64_MAX-1}, {0, INT64_MAX-1}});
  std::vector<std::string> attribute_names;
  std::vector<std::type_index> types;
  std::vector<int> num_vals;
  
  for(int i = 0; i < m_vid_mapper.get_num_fields(); i++) {
    auto& inf = m_vid_mapper.get_field_info(i);
    std::cout << inf.m_name << std::endl;

    attribute_names.push_back(inf.m_name);
    types.push_back(inf.get_tiledb_type().get_tuple_element_type_index(0u));
    num_vals.push_back(inf.m_length_descriptor.is_fixed_length_field()
                       ? inf.m_length_descriptor.get_num_elements()
                       : TILEDB_VAR_NUM);
  }

  std::vector<int> compression_types(types.size(), TILEDB_NO_COMPRESSION);
  std::vector<int> compression_levels(types.size(), 0);

  std::cout << "attributes and types " << std::endl;
  std::cout << attribute_names.size() << std::endl;
  std::cout << types.size() << std::endl;

  //m_array_schema = VariantArraySchema(array_name, attribute_names, dim_names, dim_domains, types, num_vals, compression_types, compression_levels);
  auto m_array_schema_ptr = &m_array_schema;
  m_vid_mapper.build_tiledb_array_schema(m_array_schema_ptr, array_name, false, TILEDB_NO_COMPRESSION, 0, true);

//=================================================================
  std::cout << "Line " << __LINE__ << std::endl;

  m_storage_manager = std::make_shared<VariantStorageManager>(workspace, m_segment_size, enable_shared_posixfs_optimizations());
  if (delete_and_create_tiledb_array()) {
    m_storage_manager->delete_array(array_name);
  }

  std::cout << "Line " << __LINE__ << std::endl;

  //Open array in write mode
  m_array_descriptor = m_storage_manager->open_array(array_name, &m_vid_mapper, "w");

  std::cout << "Line " << __LINE__ << std::endl;

  std::cout << "array descriptor " << m_array_descriptor << std::endl;

  std::cout << "attributes " << m_array_schema.attribute_num() << std::endl;

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
    std::cout << "Line " << __LINE__ << std::endl;
    m_array_descriptor = m_storage_manager->open_array(array_name, &m_vid_mapper, "w");
    std::cout << "Line " << __LINE__ << std::endl;
  } else if (fail_if_updating()) {
    throw LoadOperatorException(std::string("Array ") + workspace + "/" + array_name + " exists and flag \"fail_if_updating\" is set to true in the loader JSON configuration");
  }
  if (m_array_descriptor < 0) {
    throw LoadOperatorException(std::string("Could not open TileDB array for loading") + "\nTileDB error message : " + tiledb_errmsg);
  }
  std::cout << "Line " << __LINE__ << std::endl;
  m_storage_manager->update_row_bounds_in_array(m_array_descriptor, get_row_bounds().first, std::min(get_row_bounds().second, m_vid_mapper.get_max_callset_row_idx()));
  std::cout << "Line " << __LINE__ << std::endl;
  m_storage_manager->write_column_bounds_to_array(m_array_descriptor, get_column_partition(m_idx).first, get_column_partition(m_idx).second);

//=================================================================
}

void SinglePosition2TileDBLoader::read_all() {
  std::cout << "Hello" << std::endl;

  auto workspace = get_workspace(m_idx);
  auto array_name = get_array_name(m_idx);

  std::cout << "workspace " << workspace << std::endl;
  std::cout << "array " << array_name << std::endl;

  for(int i = 0; i < m_vid_mapper.get_num_callsets(); i++) {
    auto fi = m_vid_mapper.get_file_info(i);
    std::cout << fi.m_name << std::endl;
    std::ifstream file(fi.m_name);
    std::string str;

    while(std::getline(file, str)) {
      std::cout << str << std::endl;

      std::vector<std::string> vec;
      size_t index;
      while((index = str.find("\t")) != std::string::npos) {
        vec.push_back(str.substr(0, index));
        str.erase(0, index + 1);
      }
      vec.push_back(str);

      if(vec.size() >= 5) {
        std::cout << "chrom " << vec[0] << std::endl;
        std::cout << "start " << vec[1] << std::endl;
        std::cout << "end " << vec[2] << std::endl;
        std::cout << "name " << vec[3] << std::endl;
        std::cout << "score " << vec[4] << std::endl << std::endl;
      } 
    }
  }
}

