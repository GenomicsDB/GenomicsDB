#include "tiledb_loader_file_base.h"
#include "genomicsdb_config_base.h"
#include "variant_storage_manager.h"
#include "genomicsdb.h"
#include "tiledb_utils.h"

#include <htslib/sam.h>
#include <iostream>
#include <map>
#include <vector>

void read_sam_file(std::string filename);

/*struct OmicsFieldBase() {
  template <class U>
  virtual std::vector<U>& get_data() = 0;
  template <class U>
  virtual U& get_element_at(size_t idx) = 0;

  std::string name;
  enum OmicsType { omics_int64_t, omics_char, omics_string , omics_unknown_type };
}

template <class T>
struct OmicsField : OmicsFieldBase {
  OmicsField(std::string name) : name(name) {}
  OmicsField(std::string name, const std::vector<T>& data) : name(name), data(data) {}
  OmicsField(std::string name, std::vector<T>&& data) : name(name), data(data) {}
  template <class U>
  virtual std::vector<U>& get_data() override {
    assert(std::is_same(T, U));
    return data;
  }
  template <class U>
  virtual U& get_element_at(size_t idx) override {
    assert(std::is_same(T, U));
    return data[idx];
  }

  std::vector<T> data;
}*/

struct OmicsAttribute {
  OmicsAttribute(std::string name, OmicsCellAttributeType type) : name(name), type(type) {}

  std::string name;
  enum OmicsCellAttributeType { omics_char, omics_uint8_t, omics_int8_t,
                                omics_uint16_t, omics_int16_t, omics_uint32_t,
                                omics_int32_t, omics_uint64_t, omics_int64_t };
  OmicsCellAttributeType type;
}

struct OmicsCell {
  OmicsCell(std::map<std::string, std::vector<uint8_t>> attributes = {}) : attributes(attributes) {}

  std::map<std::string, std::vector<uint8_t>> attributes;
};

struct OmicsMultiCell {
  std::vector<uint8_t> as_cell(const VariantArraySchema& schema);
  
  int64_t coords[2];
  std::vector<OmicsCell> subcells;
};

class OmicsFileReader {
  public:
    OmicsFileReader(std::string filename, int file_idx, VidMapper& vid_mapper) : /*m_file(filename),*/ m_filename(filename), m_file_idx(file_idx), m_vid_mapper(vid_mapper) {
      m_buffer = new char[m_buffer_size];
      m_file_size = TileDBUtils::file_size(filename);
    }
    ~OmicsFileReader() {
      delete[] m_buffer;
    }
    // virtual std::pair<transcriptomics_cell, transcriptomics_cell> next_cell_info() = 0;
    /* static transcriptomics_cell create_end_cell(transcriptomics_cell c) {
      c.file_idx = -1;
      c.position = (c.position == c.start) ? c.end : c.start;
      return c;
    }*/
    OmicsMultiCell next_cell_info();
  protected:
    std::string m_filename;
    std::ifstream m_file;
    int m_file_idx; // index in vector of OmicsFileReaders, used to read again from same file after processing input
    ssize_t m_file_size = 0;
    ssize_t m_chars_read = 0;
    VidMapper& m_vid_mapper;
    bool generalized_getline(std::string& retval);

    const int m_buffer_size = 512;
    char* m_buffer;
    std::string m_str_buffer;
};

// intervals are represented as start and end cells, but have no special query
// support as intervals cannot always be prevented from overlapping
class GenericTileDBLoader : public GenomicsDBImportConfig {
  public:
    GenericTileDBLoader(
      const std::string& config_filename,
      const int idx,
      const bool superimpose = false
    );
    virtual void read_all() = 0;// import data from callsets
    virtual void read_array() = 0; // query
  protected:
    bool m_superimpose; // whether to contain data for multiple logical cells within one cell
    std::string m_config_filename;
    int m_idx;
    VariantArraySchema m_array_schema;
    int m_array_descriptor;
    std::shared_ptr<VariantStorageManager> m_storage_manager;
    // bool info_to_cell(const transcriptomics_cell& tc);
    //cell_info next_cell_info(std::ifstream& file, int type, int ind);
    // std::map<std::string, std::pair<long, long>> transcript_map;
    // void read_uncompressed_gtf(std::istream& input, std::string format);
    // void read_compressed_gtf(std::string filename);
    // void serialize_transcript_map(std::ostream& output);
    // void deserialize_transcript_map(std::istream&);
};
