/**
 * @file genomicsdb_plink.h
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
 * Interface to GenomicsDB plink support
 *
 **/

#pragma once

#include "genomicsdb.h"

// TODO: FIXME: Remove the following includes
#include "genomicsdb_logger.h"
#include "tiledb.h"
#include "tiledb_utils.h"
#include "variant_query_config.h"
#include "vid_mapper.h"

// TODO: FIXME: Remove references to internal classes, so the api remains opaque

class GENOMICSDB_EXPORT GenomicsDBPlinkProcessor : public GenomicsDBVariantProcessor {
  public:
    // the formats variable encodes which formats to produce by the values of specific bits: 0 - bgen, 1 - bed, 2 - tped, e.g. formats == 0b110 encodes bed and tped
    GenomicsDBPlinkProcessor(VariantQueryConfig* qc,
                             const std::string& array,
                             unsigned char formats = 7,
                             int compression = 1,
                             bool verbose = false,
                             double progress_interval = -1,
                             std::string prefix = "output",
                             std::string fam_list = "",
                             int rank = 0);

    ~GenomicsDBPlinkProcessor() {
      TileDBUtils::finalize_codec(codec);
    }

    virtual void process(const Variant& variant);
    virtual void process(const std::string& sample_name,
                         const int64_t* coordinates,
                         const genomic_interval_t& genomic_interval,
                         const std::vector<genomic_field_t>& genomic_fields,
                         const bool phased);

    void advance_state();
    const char BGEN_MASK = 1;
    const char  BED_MASK = 2;
    const char TPED_MASK = 4;
  private:
    const std::string& array;
    const VidMapper& vid_mapper;
    bool make_bgen, make_tped, make_bed;
    int compression = 0; // 0 for none, 1 for zlib, 2 for zstd
    bool verbose = false;
    // flattened coordinate to place in sorted map, phased status of column for bgen purposes (entire column considered unphased if any are unphased)
    std::map<uint64_t, std::pair<int, bool>> variant_map;
    double progress_interval;
    std::string fam_list;
    std::string prefix;
    VariantQueryConfig* query_config;
    // row to place in sorted map and sample name
    std::map<uint64_t, std::pair<int, std::string>> sample_map;
    bool sample_map_initialized = false;
    // fam is identical to tfam, used with bed, tped respectively
    std::fstream tped_file, fam_file, bim_file, bed_file, bgen_file;
    int state = 0;
    int last_sample = -1;
    int num_variants = 0;
    int last_coord = -1;
    int last_alleles = -1;
    bool last_phased;
    int rank;
    int total_rows = 0;
    int total_cols = 0;

    size_t tm = 0;
    void progress_bar(const int64_t* coords) {
      size_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      if (now - progress_interval > tm) {
        // flatten column and row ranges
        int row = 0, col = 0;
        for(auto& a : query_config->get_query_row_ranges(rank)) {
          if(coords[0] <= a.second) {
            row += coords[0] - a.first + 1;
            break;
          }
          else {
            row += a.second - a.first + 1;
          }
        }
        for(auto& b : query_config->get_query_column_ranges(rank)) {
          if(coords[1] <= b.second) {
            col += coords[1] - b.first + 1;
            break;
          }
          else {
            col += b.second - b.first + 1;
          }
        }

        long num = (long)col * total_rows + row + ((bool)state)*((long)total_rows * total_cols);
        long den = (long)total_rows * total_cols * 2;

        logger.info("Plink progress {} / {} = {:.2f}%", num, den, 100*(double)num/den);

        tm = now;
      }
    }

    // populate variant data block in bgen file for variant with given rsid spanning genomic_interval
    // vec is REF and ALT combined (REF is first entry)
    void bgen_variant_data_block(const std::string& rsid, const genomic_interval_t& genomic_interval, const std::vector<std::string>& vec, bool phased);

    // return vector of tokens that were separated by sep in str
    std::vector<std::string> split(std::string str, std::string sep = ",") {
      std::vector<std::string> retval;
      size_t index;
      if(str.length() >= 2) {
        if(str[0] == '[') {
           str = str.substr(1, str.length() - 2);
        }
      }
      while((index = str.find(sep)) != std::string::npos) {
        retval.push_back(str.substr(0, index));
        str.erase(0, index + 1);
      }
      retval.push_back(str);
      return retval;
    };

    // BED variables/functions
    char bed_buf = 0;
    char bed_buf_state = 0;
    void flush_to_bed() {
      if(bed_buf_state) {
        bed_file.write(&bed_buf, 1);
        bed_buf = 0;
        bed_buf_state = 0;
      }
    }

    void write_to_bed(char x) {
      bed_buf += x << bed_buf_state * 2;
      ++bed_buf_state %= 4;    
      if(!bed_buf_state){
        bed_file.write(&bed_buf, 1);
        bed_buf = 0;
      }
    }
    // BGEN variables
    char min_ploidy, max_ploidy;
    int32_t bgen_gt_size;
    uint32_t samples_in_column = 0;
    void *codec;
    std::string codec_buf;

    // FIXME: hard coded for B = 8
    // callback expects GT vector
    void bgen_enumerate_phased(int ploidy, int alleles, std::function<void(const std::vector<int>&, int)> callback, bool drop_last = true) {
      int size = 0;
      int ind = -1;
      int limit = alleles - drop_last;
      for(int i = 0; i < ploidy; i++) {
        for(int j = 0; j < limit; j++) {
          callback({i, j}, ++ind);
        }
      }
    }

    // get_probabilitiles expects allele counts
    void bgen_enumerate_unphased(int ploidy, int alleles, std::function<void(const std::vector<int>&, size_t)> callback, bool drop_last = true) {
      int size = 0;
      int ind = -1;
      std::vector<int> allele_counts(alleles);

      std::function<void(int, int)> enumerate_unphased;
      enumerate_unphased = [&] (int used, int depth) {
        int limit = ploidy - used - ((depth == alleles - 1) && drop_last); // if the highest depth (rightmost) do not iterate to highest possible count in order to drop last ((0, 0, ..., X) is the last)
        if(depth) {
          for(int i = 0; i <= limit; i++) { 
            allele_counts[depth] = i;
            enumerate_unphased(used + i, depth - 1);
          }
        }
        else {
          allele_counts[depth] = ploidy - used;
          ++ind;
          callback(allele_counts, (size_t)ind);
        }
      };

      enumerate_unphased(0, alleles - 1);    
    }

    void bgen_empty_cell(int ploidy, int alleles, bool phased) {
      auto write_zero = [&] (const std::vector<int>& v, int) {
        char z = 0;
        codec_buf.push_back(0);
      };
      if(phased) {
        bgen_enumerate_phased(ploidy, alleles, write_zero);
      }
      else {
        bgen_enumerate_unphased(ploidy, alleles, write_zero);
      }
    }

    // fill in size, min/max ploidy of last column
    void bgen_finish_gt() {
      // write min ploidy
      codec_buf[7] = min_ploidy;      

      // write max ploidy
      codec_buf[8] = max_ploidy;

      size_t uncompressed_size = codec_buf.size(), data_size = codec_buf.size();

      if(compression) {
        char* data;
        TileDBUtils::compress(codec, (unsigned char*)codec_buf.c_str(), codec_buf.length(), (void**)&data, data_size);
        size_t total_size = data_size + 4;
        bgen_file.write((char*)&total_size, 4); // BGEN: size of previous gt probability data plus D field
        bgen_file.write((char*)&uncompressed_size, 4);
        bgen_file.write(data, data_size);
      }
      else {
        bgen_file.write((char*)&uncompressed_size, 4); // BGEN: size of previous gt probability data plus D field
        bgen_file.write(codec_buf.c_str(), codec_buf.length());
      }

      codec_buf.clear();
      bgen_gt_size = 0;
      min_ploidy = 64;
      max_ploidy = -1;
    }

    // locations in file
    int bgen_gt_size_offset;
    int bgen_min_ploidy_offset;
    int bgen_max_ploidy_offset;
    int bgen_ploidy_info_offset;
    int bgen_probability_offset;
};
