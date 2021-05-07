#include "annotation_service.h"

#include "genomicsdb_export_config.pb.h"
#include "htslib/vcf.h"
#include <fstream>

/**
Default constructor for AnnotationService
*/
AnnotationService::AnnotationService() {
	// jDebug: an initial size is necessary, otherwise the first value gets corrupted for some reason.
	buffer = std::vector<std::string>(4);
}

AnnotationService::~AnnotationService() {
  // jDebug: Is there any destruction required for a std::vector?
  buffer.clear();
}

/**
  Add an annotation value to the buffer so it won't go out of scope
*/
void AnnotationService::copy_to_buffer(const char* value, int32_t length) {
	std::string value_string(value, length);
  buffer.push_back(std::move(value_string));
}

/**
  Read annotation configuration from a pointer to the ExportConfiguration. Create copies of the AnnotationSources.
*/
void AnnotationService::read_configuration(const std::string& str) {
  // Read the configuration
  genomicsdb_pb::ExportConfiguration export_config;
  export_config.ParseFromString(str);

  // Create a copy of each AnnotationSource
  for(auto i=0; i<export_config.annotation_source_size(); ++i) {
    genomicsdb_pb::AnnotationSource annotation_source = export_config.annotation_source(i);

    // std::filesystem::path vcf_file = annotation_source.filename();
    // if(!std::filesystem::exists(vcf_file)) {
    // std::ifstream vcf_file(annotation_source.filename());
    // FILE *file = fopen(annotation_source.filename().c_str());
    // if(!vcf_file.good()) {
    std::fstream fs;
    fs.open (annotation_source.filename().c_str(), std::fstream::in);
    if(!fs.is_open()) {
      std::__fs::filesystem::path cwd = std::__fs::filesystem::current_path();

      std::string message("VCF file does not exist: ");
      message.append(annotation_source.filename());
      message.append(". jDebug I am in ");
      message.append(cwd.string());
      throw GenomicsDBException(message);
    }
    fs.close();

    m_annotate_sources.push_back(annotation_source);
  }
}

/**
  Create a new genomic_field_t
 */
genomic_field_t AnnotationService::get_genomic_field(const std::string &data_source,
                                                     const std::string &info_attribute,
                                                     const char *value,
                                                     const int32_t value_length) {
  std::string genomic_field_name = data_source + DATA_SOURCE_FIELD_SEPARATOR + info_attribute;

  copy_to_buffer(value, value_length);
  genomic_field_t genomic_annotation(genomic_field_name, buffer.back().c_str(), buffer.back().length());
  return genomic_annotation;
}

/**
  Annotate a genomic location using all of the configured data sources. Annotations
  are stored in the genomic_fields vector with the name of the field set to the result of
  concatonating the dataSource (ie ClinVar) with a separator (see AnnotationService.DATA_SOURCE_FIELD_SEPARATOR),
  and the INFO field label.
 */
void AnnotationService::annotate(genomic_interval_t& genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields) {
  for(genomicsdb_pb::AnnotationSource annotation_source: m_annotate_sources) {
    // Open the VCF file
    htsFile * vcfInFile = hts_open(annotation_source.filename().c_str(), "r");
    VERIFY(vcfInFile != NULL && "Unable to open VCF file");

    // Check file type
    enum htsExactFormat format = hts_get_format(vcfInFile)->format;
    VERIFY(format == vcf && "File is not VCF");

    // Read the header
    bcf_hdr_t *hdr = bcf_hdr_read(vcfInFile);
    VERIFY(hdr && "Could not read VCF header");

    // I'm not sure what this does. Whatever it is, it takes a really long time.
    // Need to look in to how to remove this variable.
    regidx_t *reg_idx = NULL;
    // reg_idx = regidx_init(annotation_source.filename().c_str(), NULL, NULL, 0, NULL);
    // VERIFY(reg_idx && "Unable to read file");

    // Read the tbi index file
    tbx_t *tbx = tbx_index_load3(annotation_source.filename().c_str(), NULL, 0);
    VERIFY(tbx && "Could not load .tbi index");

    // Query using chromosome and position range
    std::string variantQuery = genomic_interval.contig_name + ":" + std::to_string(genomic_interval.interval.first) + "-" + std::to_string(genomic_interval.interval.second);

    hts_itr_t *itr = tbx_itr_querys(tbx, variantQuery.c_str());
    VERIFY(itr && "Problem opening iterator");

    const char **seq = NULL;
    int nseq;
    if (reg_idx) {
      seq = tbx_seqnames(tbx, &nseq);
    }

    kstring_t str = {0,0,0};
    bcf1_t *rec    = bcf_init1();
    int idx = 0;

    // jDebug: it isn't necessary to assign this to NULL is it?
    char* infoValue = NULL;

    int32_t infoValueLength = 0;

    // Iterate over each matching position in the VCF
    while (tbx_itr_next(vcfInFile, tbx, itr, &str) >= 0)
    {
      if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) {
        continue;
      }

      int readResult = vcf_parse1(&str, hdr, rec);
      VERIFY(readResult==0 && "Problem parsing current line of VCF");

      // jDebug: update this comment. It says "need to add code here", and i think you already have.
      // The query is constrained by chromosome and position but not by allele.  Need to add code here which
      // compares the matches alt and ref alleles to see if they match what we are looking for.
      bcf_unpack((bcf1_t*)rec, BCF_UN_ALL); // Using BCF_UN_INFO is probably a little faster

      if(ref.compare(rec->d.allele[0]) != 0) {
        // REF doesn't match
        continue;
      } else if(alt.compare(rec->d.allele[1])) {
        // ALT doesn't match
        // NOTE: This only looks at one allele. If there are multiple alternate alleles you need to update the code.
        continue;
      }

      // iteration over the list of INFO fields we are interested in
      for(auto info_attribute: annotation_source.attributes()) {
        // bcf_info_t *info = rec->d.info;

          // If the request is for "ID" then give the VCF row ID
          if(info_attribute == "ID") {
            genomic_fields.push_back(get_genomic_field(annotation_source.data_source(), info_attribute, rec->d.id, strlen(rec->d.id)));
            // jDebug: there is a problem with ID value.
            printf("jDebug.cc.ID: %s %s=%s\n", annotation_source.data_source().c_str(), info_attribute.c_str(), rec->d.id);
          } else {
            // Look for the info field
            int res = bcf_get_info_string(hdr, rec, info_attribute.c_str(), &infoValue, &infoValueLength);
            if(res > 0) {
              genomic_field_t jDebugGenomicField = get_genomic_field(annotation_source.data_source(), info_attribute, infoValue, infoValueLength);
              // genomic_fields.push_back(get_genomic_field(annotation_source.data_source(), info_attribute, infoValue, infoValueLength));
              genomic_fields.push_back(jDebugGenomicField); // jDebug: delete this line and restore the one above
              printf("jDebug.cc.a: %s %s=%s\n", annotation_source.data_source().c_str(), info_attribute.c_str(), infoValue);
              printf("jDebuc.cc.b:    %s=%s\n", jDebugGenomicField.name.c_str(), jDebugGenomicField.str_value().c_str());
            } else if (res == -2) {
              // jDebug: figure out what to do about this.
              // "clash between types defined in the header and encountered in the VCF record"
              // We need to modify this section so it determines the info field type and then uses the correct
              // bcf_get_info_* method.
              printf("JDEBUG: Need to deal with error -2, info_attribute=%s\n", info_attribute.c_str());
            } else {
              // jDebug: throw an exception? No, i think it is ok for field not to be found. 
              // info field not found
              printf("JDEBUG: Need to deal with ELSE error, info_attribute=%s\n", info_attribute.c_str());
            }
          }
      }
   }

   regidx_destroy(reg_idx);

   // jDebug: this function wasn't found so there's probably a leak: regitr_destroy(itr);
   // jDebug: need to find out if the iterator needs to be destroyed.
   bcf_itr_destroy(itr);

    int ret = hts_close(vcfInFile);
    VERIFY(ret == 0 && "Non-zero status when trying to close VCF file");
  }
}
