#include "annotation_service.h"

#include "genomicsdb_export_config.pb.h"
#include "htslib/vcf.h"

/**
Default constructor for AnnotationService
*/
AnnotationService::AnnotationService() {
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
    m_annotate_sources.push_back(annotation_source);
  }
}

/**
  Create a new genomic_field_t
 */
genomic_field_t AnnotationService::get_genomic_field(const std::string &data_source,
                                                     const std::string &info_attribute,
                                                     const char *value,
                                                     const int32_t value_length) const {
  std::string genomic_field_name = data_source + DATA_SOURCE_FIELD_SEPARATOR + info_attribute;
  genomic_field_t genomic_annotation(genomic_field_name, value, value_length);
  return genomic_annotation;
}

/**
  Annotate a genomic location using all of the configured data sources. Annotations
  are stored in the genomic_fields vector with the name of the field set to the result of
  concatonating the dataSource (ie ClinVar) with a separator (see AnnotationService.DATA_SOURCE_FIELD_SEPARATOR),
  and the INFO field label.
 */
void AnnotationService::annotate(genomic_interval_t& genomic_interval, std::string& ref, const std::string& alt, std::vector<genomic_field_t>& genomic_fields) const {
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
            // printf("jDebug.cc: %s ID=%s\n", annotation_source.data_source().c_str(), rec->d.id);
          } else {
            // Look for the info field
            int res = bcf_get_info_string(hdr, rec, info_attribute.c_str(), &infoValue, &infoValueLength);
            if(res > 0) {
              genomic_fields.push_back(get_genomic_field(annotation_source.data_source(), info_attribute, infoValue, infoValueLength));
              // printf("jDebug.cc: %s %s=%s\n", annotation_source.data_source().c_str(), info_attribute.c_str(), infoValue);
            } else if (res == -2) {
              // "clash between types defined in the header and encountered in the VCF record"
              // We need to modify this section so it determines the info field type and then uses the correct
              // bcf_get_info_* method.
            } else {
              // info field not found
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
