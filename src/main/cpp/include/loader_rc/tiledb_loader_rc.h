#include "tiledb_loader_file_base.h"
#include "genomicsdb_config_base.h"
#include "variant_storage_manager.h"
#include "load_operators.h"
#include "tiledb_utils.h"

// flattened start, flattened end, name, gene, score, sample index, file index
typedef std::tuple<int64_t, int64_t, std::string, std::string, float, int, int> cell_info;

/*cell_info reverse_cell_info(cell_info inf) {
  int64_t start = std::get<0>(inf);
  std::get<0>(inf) = std::get<1>(inf);
  std::get<1>(inf) = start;
  return inf;
}*/

class TranscriptomicsFileReader {
  public:
    TranscriptomicsFileReader(std::string fname, int file_idx, VidMapper& vid_mapper) : /*m_file(fname),*/ m_fname(fname), m_file_idx(file_idx), m_vid_mapper(vid_mapper) {
      m_buffer = new char[m_buffer_size];
      m_file_size = TileDBUtils::file_size(fname);
    }
    ~TranscriptomicsFileReader() {
      delete[] m_buffer;
    }
    virtual cell_info next_cell_info() = 0;
  protected:
    std::string m_fname;
    std::ifstream m_file;
    int m_file_idx;
    ssize_t m_file_size = 0;
    ssize_t m_chars_read = 0;
    VidMapper& m_vid_mapper;
    bool generalized_getline(std::string& retval);

    const int m_buffer_size = 512;
    char* m_buffer;
    std::string m_str_buffer;
};

class BedReader : public TranscriptomicsFileReader {
  public:
    BedReader(std::string fname, int ind, VidMapper& vid_mapper, int sample_idx) : TranscriptomicsFileReader(fname, ind, vid_mapper), m_sample_idx(sample_idx) { }
    cell_info next_cell_info() override;
  protected:
    int64_t m_sample_idx;
};

class MatrixReader : public TranscriptomicsFileReader {
  public:
    MatrixReader(std::string fname, int ind, VidMapper& vid_mapper, std::map<std::string, std::pair<long, long>>& transcript_map) : TranscriptomicsFileReader(fname, ind, vid_mapper), m_transcript_map(transcript_map) {
      std::string str;
      generalized_getline(str);
      std::string tok;
      std::stringstream ss(str);

      ss >> tok >> tok; // consume gene and name headers

      // collect samples
      while(ss >> tok) {
        samples.push_back(tok);
      }
    }
    cell_info next_cell_info() override;
    std::vector<std::pair<int, int64_t>> idx_to_row; // maps from column in matrix (sample) to row in GenomicsDB
    int idx_to_row_pos = 0; // keeps position in above vector
  protected:
    std::vector<int> m_cols;
    std::map<std::string, std::pair<long, long>>& m_transcript_map;
    std::vector<std::string> m_current_line;
    std::vector<std::string> samples;
    long m_current_start;
    long m_current_end;
    std::string m_current_name;
    std::string m_current_gene;
    int m_current_sample;
};


class SinglePosition2TileDBLoader : public GenomicsDBImportConfig {
  public:
    SinglePosition2TileDBLoader(
      const std::string& config_filename,
      const int idx,
      std::string gtf_name = "",
      std::string gi_name = ""
    );
    void read_all();
    void read_array();
    static std::tuple<std::string, int64_t, int64_t, std::string, float> parse_and_check(std::string, VidMapper& vid_mapper);
  private:
    int m_idx;
    VariantArraySchema m_array_schema;
    int m_array_descriptor;
    std::shared_ptr<VariantStorageManager> m_storage_manager;
    bool info_to_cell(int64_t start, int64_t end, std::string name, std::string gene, float score, uint64_t row);
    cell_info next_cell_info(std::ifstream& file, int type, int ind);
    std::map<std::string, std::pair<long, long>> transcript_map;
    void read_uncompressed_gtf(std::istream& input, std::string format);
    void read_compressed_gtf(std::string fname);
    void serialize_transcript_map(std::ostream& output);
    void deserialize_transcript_map(std::istream&);
};
