#include "tiledb_loader_file_base.h"
#include "genomicsdb_config_base.h"
#include "variant_storage_manager.h"
#include "load_operators.h"

class SinglePosition2TileDBLoader : public GenomicsDBImportConfig {
  public:
    SinglePosition2TileDBLoader(
      const std::string& config_filename,
      const int idx
    );
    void read_all();
  private:
    int m_idx;
    VariantArraySchema m_array_schema;
    int m_array_descriptor;
    std::shared_ptr<VariantStorageManager> m_storage_manager;
};
