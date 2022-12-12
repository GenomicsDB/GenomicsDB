/**
 * @file genomicsdb_plink.cc
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2022 Omics Data Automation, Inc.
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
 *
 * @section DESCRIPTION
 *
 * Implementation of the query interface to GenomicsDB and plink
 *
 **/

#include "genomicsdb_plink.h"

#include <iostream>
#include <string>

#include "annotation_service.h"
#include "query_variants.h"
#include "tiledb_utils.h"
#include "variant_field_data.h"
#include "variant_operations.h"
#include "variant_query_config.h"
#include "vcf_adapter.h"
#include "vid_mapper_pb.h"

#define TO_VARIANT_QUERY_CONFIG(X) (reinterpret_cast<VariantQueryConfig *>(static_cast<void *>(X)))

// Prototypes to internal methods in genomicsdb.cc declared here instead of header to keep the api opaque
std::map<std::string, genomic_field_type_t> create_genomic_field_types(const VariantQueryConfig &query_config,
                                                   void *annotation_service, bool change_alt_to_string=false);
void GenomicsDB::generate_plink(const std::string& array,
                                genomicsdb_ranges_t column_ranges,
                                genomicsdb_ranges_t row_ranges,
                                unsigned char format,
                                int compression,
                                bool one_pass,
                                bool verbose,
                                double progress_interval,
                                const std::string& output_prefix,
                                const std::string& fam_list) {
   // Create Variant Config for given concurrency_rank
  VariantQueryConfig *base_query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);
  VariantQueryConfig query_config(*base_query_config);
  query_config.set_array_name(array);
  if (column_ranges.size() > 0) {
    query_config.set_query_column_ranges(column_ranges);
  }
  if (row_ranges.size() > 0) {
    query_config.set_query_row_ranges(row_ranges);
  }

  query_config.validate();

  GenomicsDBPlinkProcessor proc(&query_config, array, format, compression, verbose, progress_interval,
                                output_prefix, fam_list, m_concurrency_rank);

  proc.initialize(create_genomic_field_types(query_config, m_annotation_service, true));

  if(!one_pass) { // if one_pass is true, skip first pass and get participating samples from callset (will include samples without data)
    query_variants(array, &query_config, proc);
  }
  proc.advance_state();
  query_variants(array, &query_config, proc);
  proc.advance_state();
}

void GenomicsDB::generate_plink(unsigned char format,
                                int compression,
                                bool one_pass,
                                bool verbose,
                                double progress_interval,
                                const std::string& output_prefix,
                                const std::string& fam_list) {
  VariantQueryConfig* query_config = TO_VARIANT_QUERY_CONFIG(m_query_config);

  auto array = query_config->get_array_name(m_concurrency_rank);
  GenomicsDBPlinkProcessor proc(query_config, array, format, compression, verbose, progress_interval,
                                output_prefix, fam_list, m_concurrency_rank);

  proc.initialize(create_genomic_field_types(*query_config, m_annotation_service, true));

  if(!one_pass) { // if one_pass is true, skip first pass and get participating samples from callset (will include samples without data)
    query_variants(array, query_config, proc);
  }
  proc.advance_state();
  query_variants(array, query_config, proc);
  proc.advance_state();
}

GenomicsDBPlinkProcessor::GenomicsDBPlinkProcessor(VariantQueryConfig* qc,
                                                   const std::string& array,
                                                   unsigned char formats,
                                                   int compression,
                                                   bool verbose,
                                                   double progress_interval,
                                                   std::string prefix,
                                                   std::string fam_list,
                                                   int rank
                                                   ) : query_config(qc), vid_mapper(qc->get_vid_mapper()), array(array), compression(compression), verbose(verbose), progress_interval(progress_interval), prefix(prefix), fam_list(fam_list), rank(rank) {
  make_bgen = (bool)(formats & BGEN_MASK);
  make_bed = (bool)(formats & BED_MASK);
  make_tped = (bool)(formats & TPED_MASK);

  // For use with BGEN compression
  if(compression == 1) {
    TileDBUtils::create_codec(&codec, TILEDB_GZIP, Z_DEFAULT_COMPRESSION);
  }
  else {
    TileDBUtils::create_codec(&codec, TILEDB_ZSTD, 9);
  }

  // open various files
  if(make_tped) {
    tped_file.open(prefix + ".tped", std::ios::out);
  }
  if(make_bed) {
    bed_file.open(prefix + ".bed", std::ios::out | std::ios::binary);
    bim_file.open(prefix + ".bim", std::ios::out);
  }
  if(make_tped || make_bed) {
    fam_file.open(prefix + ".fam", std::ios::out);
  }
  if(make_bgen) {
    bgen_file.open(prefix + ".bgen", std::ios::out | std::ios::binary);
  }

  if(make_bed) {
    char magic_numbers[] = {0x6c, 0x1b, 0x01};
    bed_file.write(magic_numbers, 3); // BED: magic numbers
  }

  if(make_bgen) {
    int32_t zero = 0;
    int32_t offset = 20; // BGEN: offset of first variant data block relative to 5th byte. Always 20 here because free data area left empty
    bgen_file.write((char*)&offset, 4); // BGEN: first 4 bytes has offset
    bgen_file.write((char*)&offset, 4); // BGEN: beginning of header, next 4 bytes same in this case
    bgen_file.write((char*)&zero, 4); // BGEN: 4 bytes number of variants (M), filled in later
    bgen_file.write((char*)&zero, 4); // BGEN: 4 bytes number of samples (N), filled in later

    char bgen_magic_numbers[] = {'b', 'g', 'e', 'n'};
    bgen_file.write(bgen_magic_numbers, 4); // BGEN: 4 bytes bgen magic number

    int32_t flags = 0b10000000000000000000000000001000; // BGEN: layout 2, sample identifier block present
    flags = flags | compression;
    bgen_file.write((char*)&flags, 4); // BGEN: 4 bytes flags, end of header
  }

  // for use with progress bar
  for(auto& a : query_config->get_query_row_ranges(rank)) {
    total_rows += a.second - a.first + 1;
  }
  for(auto& b : query_config->get_query_column_ranges(rank)) {
    total_cols += b.second - b.first + 1;
  }
}

void GenomicsDBPlinkProcessor::bgen_variant_data_block(const std::string& rsid, const genomic_interval_t& genomic_interval, const std::vector<std::string>& vec, bool phased) {
  // BGEN: fill in genotype block size and min/max ploidy from previous iteration
  if(last_sample != -1) { // no need on first column
    if(samples_in_column < sample_map.size()) {
      min_ploidy = (min_ploidy > 2) ? 2 : min_ploidy;
      max_ploidy = (max_ploidy < 2) ? 2 : max_ploidy;
    }
    samples_in_column = 0;
    bgen_finish_gt();
  }
  min_ploidy = 64;
  max_ploidy = -1;

  // BGEN: variant data blocks
  int16_t zero = 0;
  int16_t one = 1;
  bgen_file.write((char*)&one, 2); // BGEN: 2 byte length of variant identifier, not stored in GenomicsDB so using dummy
  bgen_file.write((char*)&one, 1); // BGEN: dummy variant id
  int16_t rsid_len = rsid.length();
  bgen_file.write((char*)&rsid_len, 2); // BGEN: 2 byte length of rsid
  bgen_file.write(rsid.c_str(), rsid_len); // BGEN: rsid
  std::string chrom = genomic_interval.contig_name;
  int16_t chrom_len = chrom.length();
  bgen_file.write((char*)&chrom_len, 2); // BGEN: 2 byte chrom length
  bgen_file.write(chrom.c_str(), chrom_len); // BGEN: chrom
  uint32_t variant_pos = genomic_interval.interval.first;
  bgen_file.write((char*)&variant_pos, 4); // BGEN: 4 byte variant position
  int16_t K = vec.size();
  bgen_file.write((char*)&K, 2);// BGEN: 2 byte K for number of alleles
  // write alleles and lengths
  int32_t len;
  for(auto& a : vec) { // iterate over alleles
    len = a.length();
    bgen_file.write((char*)&len, 4); // BGEN: 4 byte length of allele
    bgen_file.write(a.c_str(), len); // BGEN: allele
  }

  // BGEN: genotype data block (layout 2, uncompressed)
  bgen_gt_size_offset = bgen_file.tellp();
  bgen_gt_size = 0;
  int32_t fourB_zero = 0;

  // BGEN: preallocate probability data storage
  int32_t N = sample_map.size();
  codec_buf.append((char*)&N, 4); // BGEN: 4 byte N
  codec_buf.append((char*)&K, 2); // BGEN: 2 byte K
  bgen_min_ploidy_offset = bgen_file.tellp();
  codec_buf.append((char*)&fourB_zero, 1); // BGEN: 1 byte min ploidy (placeholder)
  bgen_max_ploidy_offset = bgen_file.tellp();
  codec_buf.append((char*)&fourB_zero, 1); // BGEN: 1 byte max ploidy (placeholder)
  bgen_ploidy_info_offset = bgen_file.tellp();
  char default_sample_info = 0b10000010; // default missingness and ploidy information: set to missing/diploid unspecified
  for(int j = 0; j < N; j++) {
    codec_buf.append(&default_sample_info, 1); // BGEN: default sample information within this variant: because missing is set to 1, no need to backfill skipped cells
  }
  codec_buf.append((char*)&phased, 1); // BGEN: 1 byte phased
  char B = 8; // precision at one byte to avoid alignment difficulties
  codec_buf.append(&B, 1); // BGEN: 1 byte unsigned bits of precision
  bgen_probability_offset = bgen_file.tellp();
}

void GenomicsDBPlinkProcessor::process(const Variant& variant) {
  auto& calls = variant.get_calls();

  bool phased = true;
  std::string gt_string;

  for(auto& c : calls) {
    auto fields = get_genomic_fields_for(array, &c, query_config);
    
    for(auto& f : fields) {
      if(f.name == "GT") {
        gt_string = f.to_string(get_genomic_field_type(f.name));

        if(gt_string.find('/') != std::string::npos && gt_string.find('.') == std::string::npos) { // unphased if contains a slash or is ./. (might be able to use GL/PL, which are unphased probabilities)
          phased = false;
        }
        break;
      }
    }
    if(!phased) { break; }
  }

  for(auto& c : calls) {
    auto fields = get_genomic_fields_for(array, &c, query_config);

    int64_t coords[] = {(int64_t)c.get_row_idx(), (int64_t)c.get_column_begin()};
    int64_t end_position = c.get_column_end();

    std::string contig_name;
    int64_t contig_position;
    if (!vid_mapper.get_contig_location(coords[1], contig_name, contig_position)) {
      std::cerr << "Could not find genomic interval associated with Variant(Call) at "
        << coords[1] << std::endl;
      continue;
    }

    contig_position++;
    genomic_interval_t genomic_interval(std::move(contig_name),
                                        std::make_pair(contig_position, contig_position+end_position-coords[1]));

    std::string sample_name;
    if (!vid_mapper.get_callset_name(coords[0], sample_name)) {
      sample_name = "NONE";
    }

    process(sample_name, coords, genomic_interval, fields, phased);
  }
}

void GenomicsDBPlinkProcessor::process(const std::string& sample_name,
                                       const int64_t* coords,
                                       const genomic_interval_t& genomic_interval,
                                       const std::vector<genomic_field_t>& genomic_fields,
                                       const bool phased) {
  last_phased = phased;

  if(progress_interval > 0) {
    progress_bar(coords);
  }

  std::string ref_string, alt_string, gt_string, id_string, pl_string, gl_string, pq_string;
  for(auto& f : genomic_fields) {
    if(f.name == "ALT") {
      std::string combined_alt = f.recombine_ALT_value(get_genomic_field_type(f.name));
      if(combined_alt.size()) {
        alt_string = combined_alt.substr(1, combined_alt.length() - 2);
      }
    }
    else if(f.name == "REF") {
      ref_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "GT") {
      gt_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "ID") {
      id_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "PL") {
      pl_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "GL") {
      gl_string = f.to_string(get_genomic_field_type(f.name));
    }
    else if(f.name == "PQ") {
      pq_string = f.to_string(get_genomic_field_type(f.name));
    }
  }

  if(!alt_string.size()) {
    logger.error("No ALT field for sample: {} row: {} column: {}", sample_name, coords[0], coords[1]);
    exit(1);
  }
  if(!ref_string.size()) {
    logger.error("No REF field for sample: {} row: {} column: {}", sample_name, coords[0], coords[1]);
    exit(1);
  }
  if(!gt_string.size()) {
    logger.error("No GT field for sample: {} row: {} column: {}", sample_name, coords[0], coords[1]);
    exit(1);
  }

  // collect alleles
  std::vector<std::string> vec = {ref_string};
  size_t index;
  while((index = alt_string.find(", ")) != std::string::npos) {
    vec.push_back(alt_string.substr(0, index));
    alt_string.erase(0, index + 2);
  }
  vec.push_back(alt_string);

  std::vector<int> gt_vec;
  auto iter = gt_string.begin();

  // parse GT
  try {
    while((iter = find_if(gt_string.begin(), gt_string.end(), [] (char c) { return c == '|' || c == '/'; })) != gt_string.end()) {
      index = iter - gt_string.begin();
      gt_vec.push_back(std::stoi(gt_string.substr(0, index)));
      gt_string.erase(0, index + 1);
    }
    gt_vec.push_back(std::stoi(gt_string));
  }
  catch (...) { // probably ./., treat as missing cell
    return;
  }

  for(auto g : gt_vec) {
    if(g < 0 || g >= vec.size()) {
      if(verbose) {
        logger.error("GT field for sample: {} row: {} column: {},  contains index {}, which is out of bounds (1 ref allele, {} alt allele(s))", sample_name, coords[0], coords[1], g, vec.size() - 1);
       }
      return;
    }
  }

  int16_t ploidy = gt_vec.size();

  if(!state) {
    sample_map_initialized = true; // possible to skip first pass, sample map will be populated from callset

    sample_map.insert(std::make_pair(coords[0], std::make_pair(-1, sample_name)));
    last_coord = coords[1];
    return;
  }

  if(ploidy != 2 && (make_bed || make_tped)) { 
    logger.error("The tped/bed formats do not support ploidies other than 2.");
    make_bed = 0;
    make_tped = 0;
    if(make_bgen) {
      logger.info("Continuing bgen generation");
    } 
  }

  if(state == 1) {
    std::string rsid;
    std::string rsid_row;
    if(id_string.size()) {
      rsid = id_string;
    }
    else {
      rsid = genomic_interval.contig_name + ":" + std::to_string(genomic_interval.interval.first);
    }
    rsid_row = genomic_interval.contig_name + " " + rsid + " 0 " + std::to_string(genomic_interval.interval.first);

    int sind = sample_map[coords[0]].first;

    // backfill if needed
    bool backfill = coords[1] - last_coord && last_coord != -1;
    int add_to_prev = backfill ? sample_map.size() - (last_sample + 1) : 0;
    int add_to_current = backfill ? sind : sind - (last_sample + 1);

    if(last_coord != coords[1]) { // first in variant
      num_variants++;

      for(int i = 0; i < add_to_prev; i++) { // backfill samples missing from previous variant
        if(make_tped) {
          tped_file << " 0 0";
        }
        if(make_bed) {
          write_to_bed(1);
        }
        if(make_bgen) {
          // BGEN: backfill probability data
          bgen_empty_cell(2, last_alleles, phased);
        }
      }
      if(make_bed) {
        flush_to_bed();
      }

      if(make_tped) {
        if(last_sample != -1) { // first line should not have newline
          tped_file << std::endl;
        }
        tped_file << rsid_row;
      }

      if(make_bgen) {
        bgen_variant_data_block(rsid, genomic_interval, vec, phased);
      }
    }
    else {
      samples_in_column++;
    }

    for(int i = 0; i < add_to_current; i++) { // backfill samples missing from current variant
      if(make_bed) {
        write_to_bed(1);
      }

      if(make_tped) {
        tped_file << " 0 0";
      }

      if(make_bgen) {
        // BGEN: backfill probability data
        bgen_empty_cell(2, vec.size(), phased);
      }
    }

    // safe to update now that backfilling is over
    last_alleles = vec.size();

    if(make_bed) {
      char gt_code;
      if(!(gt_vec[0] || gt_vec[1])) { // homozygous for first allele
        gt_code = 0;
      }
      else if(gt_vec[0] + gt_vec[1] == 1) { // heterozygous
        gt_code = 2;
      }
      else if(gt_vec[0] == 1 && gt_vec[1] == 1) { // homozygous for second allele
        gt_code = 3;
      }
      else { // missing genotype
        gt_code = 1;
      }

      if(last_coord != coords[1]) {
        bim_file << rsid_row;
        bim_file << " " << vec[0] << " " << vec[1] << std::endl;
      }
      write_to_bed(gt_code);
    }

    if(make_tped) {
      tped_file << " " << vec[gt_vec[0]] << " " << vec[gt_vec[1]];
    }

    last_sample = sample_map[coords[0]].first;
    //last_variant = vind;
    last_coord = coords[1];

    // convert PL to BGEN format
    std::vector<double> pl_vec;
    bool pl_dot = false;

    try {
      for(auto& tok : split(pl_string)) {
        if(tok == ".") pl_dot = true;
        int val = std::stoi(tok);
        if(val >= 0) {
          pl_vec.push_back(val);
        }
        else {
          throw std::runtime_error("PL value is negative");
        }
      }
    }
    catch(...) {
      pl_vec.clear();
      if(!pl_dot && verbose) {
        logger.error("PL field for sample: {} row: {} column: {} contains a negative value or otherwise did not parse, full string: {}", sample_name, coords[0], coords[1], pl_string);
      }
    }

    std::vector<char> pl_probs;

    double pl_total = 0;
    double epsilon = .05;
    for(auto& a : pl_vec) {
      double prob = std::pow(10, double(a)/-10);
      pl_total += prob;
      pl_probs.push_back(char(std::numeric_limits<unsigned char>::max() * prob));
    }

    // parse GL
    std::vector<double> gl_vec;
    bool gl_dot = false;

    try {
      for(auto& tok : split(gl_string)) {
        if(tok == ".") gl_dot = true;
        double val = std::stod(tok);
        if(val <= 0) {
          gl_vec.push_back(val);
        }
        else {
          throw std::runtime_error("GL value is positive");
        }
      }
    }
    catch(...) {
      gl_vec.clear();
      if(!gl_dot && verbose) {
        logger.error("GL field for sample: {} row: {} column: {} contains a strictly positive value or otherwise did not parse, full string {}", sample_name, coords[0], coords[1], gl_string);
      }
    }

    std::vector<char> gl_probs;

    double gl_total = 0;
    for(auto& a : gl_vec) {
      double prob = std::pow(10, a);
      gl_total += prob;
      gl_probs.push_back(char(std::numeric_limits<unsigned char>::max() * prob));
    }

    std::vector<char> probs;

    // subtract 1 representing reference
    auto num_genotypes = KnownFieldInfo::get_number_of_genotypes(vec.size() - 1, ploidy);

    if(gl_probs.size()) { // prefer gl as it is more precise (pl is rounded)
      if(gl_probs.size() == num_genotypes) {
        if(std::abs(1 - gl_total) > epsilon && verbose) {
          logger.warn("GL probabilities at sample: {} row: {} column: {} sum to {}, expected near 1. Generated BGEN may be invalid", sample_name, coords[0], coords[1], gl_total);
        }
        probs = gl_probs;
      }
      else {
        if(verbose) {
          logger.error("GL length at sample: {} row: {} column: {} is {}, expected {}. Defaulting to using GT to construct probabilities", sample_name, coords[0], coords[1], gl_vec.size(), num_genotypes);
        }
      }
    }
    else if(pl_probs.size()) {
      if(pl_probs.size() == num_genotypes) {
        if(std::abs(1 - pl_total) > epsilon && verbose) {
          logger.warn("PL probabilities at sample: {} row: {} column: {} sum to {}, expected near 1. Generated BGEN may be invalid", sample_name, coords[0], coords[1], pl_total);
        }
        probs = pl_probs;
      }
      else {
        if(verbose) {
          logger.error("PL length at sample: {} row: {} column: {} is {}, expected {}. Defaulting to using GT to construct probabilities", sample_name, coords[0], coords[1], pl_vec.size(), num_genotypes);
        }
      }
    }

    double pq;
    if(pq_string.length()) {
      try {
        pq = std::pow(10, (double)std::stoi(pq_string)/-10);
      }
      catch(...) {
        pq_string.clear();
      }
    }

    auto write_phased_probability = [&] (const std::vector<int>& v, size_t ind) {
      char p = gt_vec[v[0]] == v[1] ? -1 : 0;
      codec_buf.push_back(p);
    };

    auto write_unphased_probability = [&] (const std::vector<int>& v, size_t ind) {
      char p;

      if(!probs.size()) {
        std::vector<int> counts(vec.size(), 0);
        for(auto& g : gt_vec) {
          counts[g]++;
        }
        p = counts == v ? -1 : 0;
      }
      else {
        if(ind < probs.size()) {
          p = probs[ind];
        }
        else {
          // NOTE: this should never be triggered, GL/PL length is checked above
          logger.error("BGEN generation error: GL/PL probabilies only have {} term(s), halting BGEN generation", probs.size());
          make_bgen = 0;
        }
      }
      codec_buf.push_back(p);
    };

    if(make_bgen) {
      // store sample as not missing/ploidy info
      if(ploidy < 64) {
         min_ploidy = ploidy < min_ploidy ? ploidy : min_ploidy;
        max_ploidy = ploidy > max_ploidy ? ploidy : max_ploidy;
        char p = ploidy;
        codec_buf[8 + sample_map[coords[0]].first] = p;

        if(phased) { // phased
          bgen_enumerate_phased(ploidy, vec.size(), write_phased_probability);
        }
         else { // unphased
          bgen_enumerate_unphased(ploidy, vec.size(), write_unphased_probability);
        }
      }
    }
  }
}

void GenomicsDBPlinkProcessor::advance_state() {
  if(!state) {
    if(!sample_map_initialized) { // populate from callset
      std::string str;
      auto num_rows = query_config->get_num_rows_to_query();
      int row;
      for(int i = 0; i < num_rows; i++) {
        row = query_config->get_array_row_idx_for_query_row_idx(i);
        vid_mapper.get_callset_name(row, str);
        sample_map.insert(std::make_pair(row, std::make_pair(-1, str)));
      }
    }

    // associate samples with sorted position
    int i = -1;
    for(auto& a : sample_map) {
      ++i;
      a.second.first = (uint64_t)i;
    }
    // associate variants with sorted position
    i = -1;
    last_sample = -1;
    last_coord = -1;
  
    if(make_bed || make_tped) {
      // Find samples coincident with entries in fam files (if specified) and associate information with sample name
      std::map<std::string, std::string> fam_entries;
    
      if(fam_list.length()) {
        std::ifstream fam_list_file(fam_list);
        std::string fname;
        while(std::getline(fam_list_file, fname)) {
          std::ifstream file(fname);
          std::string entry;
          while(std::getline(file, entry)) {
            std::string fid, wfid, fthid, mthid;
            char sex, pt;
            std::stringstream(entry) >> fid >> wfid >> fthid >> mthid >> sex >> pt;
            fam_entries.insert(std::make_pair(wfid, entry));
          }
        }
      }
      for(auto& s : sample_map) {
        if(fam_entries.count(s.second.second)) {
          fam_file << fam_entries[s.second.second] << std::endl;
        }
        else {
          fam_file << s.second.second << " " << s.second.second << " 0 0 0 0" << std::endl;
        }
      }
    }

    if(make_bgen) {
      // BGEN: fill in N in header
      bgen_file.seekp(12);
      int32_t N = sample_map.size();
      bgen_file.write((char*)&N, 4); // BGEN: 4 byte N
      // Fill in M at end, as first pass, which counts variants, may be skipped

      // BGEN: write sample identifier block
      bgen_file.seekp(24);
      int32_t lsi = 0;
      bgen_file.write((char*)&lsi, 4); // BGEN: 4 byte total length of sample identifier block, filled in last
      bgen_file.write((char*)&N, 4); // BGEN: 4 byte N, deliberately duplicated here as per spec
      lsi = 8; // total length includes 8 bytes metainfo

      int16_t len;
      for(auto& s : sample_map) { // BGEN: write each sample id, can potentially be combined with above foreach loop
        len = s.second.second.length();
        lsi += len + 2; // total length increased by length of identifier and 2 byte length field
       
        bgen_file.write((char*)&len, 2);
        bgen_file.write(s.second.second.c_str(), len);
      }

      bgen_file.seekp(24);
      bgen_file.write((char*)&lsi, 4); // BGEN: 4 byte total length of sample identifier block, now with correct value
       bgen_file.seekp(0);
      int32_t offset = lsi + 20;
      bgen_file.write((char*)&offset, 4);
      // BGEN: update initial offset to include size of sample block
      bgen_file.seekp(24 + lsi); // seek to end of sample identifier block
    }
  }

  if(state == 1) {
    // if skipped some samples at end, fill with reample_map.insert(std::make_pair(sample_name, -1));
    int add_to_current = sample_map.size() - (last_sample + 1);

    for(int i = 0; i < add_to_current; i++) {
      if(make_tped) {
        tped_file << " 0 0";
      }
      if(make_bed) {
        write_to_bed(1);
      }
      if(make_bgen) {
        // BGEN: backfill last variant
        bgen_empty_cell(2, last_alleles, last_phased);
      }
    }
    if(make_bed) {
      flush_to_bed();
    }

    if(make_tped) {
      tped_file << std::endl;
    }

    if(sample_map.size() && make_bgen) {
      bgen_finish_gt();
    }

    if(make_bgen) {
      // BGEN: fill in M in header
      bgen_file.seekp(8);
      int32_t M = num_variants;
      bgen_file.write((char*)&M, 4); // BGEN: 4 byte M
    }
  }
  state++;
}
