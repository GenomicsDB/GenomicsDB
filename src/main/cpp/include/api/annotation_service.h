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

  ~AnnotationService();

  // Add a string to the list of annotation values
  const char* copy_to_buffer(const char* value);

  // Value which separates dataSource and info field
  const std::string DATA_SOURCE_FIELD_SEPARATOR = "_AA_";

  // Read configuration from an ExportConfiguration
  void read_configuration(const std::string& str);

  //
  void annotate(genomic_interval_t &genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields);

  genomic_field_t get_genomic_field(const std::string &data_source, const std::string &info_attribute, const char *value, const int32_t value_length);

  // List of configured annotation data sources
  std::vector<genomicsdb_pb::AnnotationSource> m_annotate_sources;

  // Buffer for annotation values
  std::vector<std::string> buffer;
};
