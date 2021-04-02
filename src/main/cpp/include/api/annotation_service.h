#include "htslib/regidx.h"
#include "htslib/tbx.h"

#include "genomicsdb.h"

#include <map>
#include <string>
#include <vector>

#include "genomicsdb_export_config.pb.h"

#define VERIFY(X) if(!(X)) throw GenomicsDBException(#X);

/**
  Use this service to add annotations from VCF datasources to variant's genomic ields
*/
class AnnotationService {
 public:
  // Default constructor
  AnnotationService();

  // Buffer which associates genomic location with a list of annotation values.
  // The dataSource and info field which point to this value are stored elsewhere
  // Example: value_buffer["chr1-123456-A-G"] = { "val1", "val2", val3 };
  std::map<std::string, std::vector<std::string>> value_buffer;

  // Value which separates dataSource and info field
  const std::string DATA_SOURCE_FIELD_SEPARATOR = "_AA_";

  // Read configuration from an ExportConfiguration
  void read_configuration(const std::string& str);

  //
  void annotate(genomic_interval_t &genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields) const;

  genomic_field_t get_genomic_field(const std::string &data_source, const std::string &info_attribute, const char *value, const int32_t value_length) const;

  // List of configured annotation data sources
  std::vector<genomicsdb_pb::AnnotationSource> m_annotate_sources;
};
