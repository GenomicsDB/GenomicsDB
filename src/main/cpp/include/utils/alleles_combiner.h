/**
 * The MIT License (MIT)
 * Copyright (c) 2020-21 Omics Data Automation, Inc.
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
*/

#ifndef ALLELES_COMBINER_H
#define ALLELES_COMBINER_H

#include "headers.h"
#include "lut.h"
#include "gt_common.h"

/*
 * This module deals with everything related to REF and ALT alleles
 * a. Merging REF and ALT alleles across samples
 * b. Keep track of deletions/MNVs so that the iterator can step correctly (spanning deletions)
 * I've tried to optimize this module as much as possible since it was one of the few modules
 * that showed up individually as a big chunk of the total time (~12%). The following points are considered
 * during optimization
 * a. Single pass over REF and ALT fields - the idea is to extract all the information once
 * and not have to parse strings repeatedly. Keeping track of deletion/MNV & NON_REF indexes is one example
 * where we avoid parsing strings in downstream processing (reoredering fields, spanning deletions etc).
 * b. While merging REF-ALT alleles, the simplest method is the following:
 *   i. Set the merged REF to the longest REF across samples
 *  ii. Append suffix to each ALT allele - e.g merged REF allele is AGCT and sample contains A->AT. The ALT
 *      allele gets a suffix GCT and becomes ATGCT (the variant is represented as AGCT->ATGCT)
 * iii. Put the modified alleles into a map to track unique alleles
 *   This method can be improved in a couple of ways:
 *   i. Since each ALT allele gets a suffix, a copy from GenomicsDBBuffer is required (into a string for example).
 *      Hence, I store a "normalized" version of each allele which doesn't require any suffix. Here are the rules
 *      for normalized alleles
 *      SNV: single REF base, single ALT base. So, if a sample has SNV ATGC->TTGC, the trailing suffix gets dropped
 *      Symbolic alleles and insertions: single REF base. AGCT->ATGCT gets converted to A->AT
 *      Deletions: ALT has a single base. ATGC->AC gets converted to ATG->A
 *      MNVs: Last element of the REF and ALT must be different. So ATGC->TTAC gets converted to ATG->TTA
 *      By using this convention, no copies need to be made. We can use STRING_VIEW objects to refer to
 *      locations inside GenomicsDBBuffer.
 *  ii. We can exploit the properties of the normalized forms of each type of allele (SNV, deletion, MNV etc) 
 *      to make searching faster than a single map.
 *      SNVs: ALT can be A,T,G,C - so a simple lookup table of size 255 (ASCII) is enough
 *      Symbolic alleles and insertions: map<STRING_VIEW> - ALT allele as the key
 *      Deletions: map<unsigned> where the key is the number of bases from REF deleted
 *      MNVS: map<STRING_VIEW> with ALT allele as key.
 * Apart from this, the module keeps track of deletion and indexes for each sample. This will be used for
 * spanning deletions. Again, this information is kept in a way that minimizes dynamic allocation overhead
 * (see below)
 */

//This class stores the normalized form of each allele. A vector of such objects stores the merged alleles
//for the current position
class MergedAllelesVecEntry {
  public:
    MergedAllelesVecEntry(const size_t ref_length, const STRING_VIEW& allele_arg, const bool a_is_symbolic)
    : is_symbolic_allele(a_is_symbolic), REF_length(ref_length), allele(allele_arg){
    }
    bool is_symbolic_allele;
    size_t REF_length;
    STRING_VIEW allele;
};

class AlleleConfig;
class AllelesCombiner {
  public:
    AllelesCombiner(const size_t num_queried_rows);
    /*
     * Called by GenomicsDBGVCFIterator.
     * This is called right before data for samples with variants beginning at the current position
     * (current as defined by the iterator) are added. This function resets some data structures
     * used for tracking deletion indexes.
     */
    void reset_before_adding_new_sample_info_at_current_position();
    void reset_for_next_query_interval();
    /*
     * Update state based on REF and delimited ALT string for sample corresponding to row_query_idx
     * Last arg true if the interval begins before the current position (current as defined by the
     * GenomicsDBGVCFIterator)
     */
    void insert_allele_info(const size_t row_query_idx, const STRING_VIEW& REF,
        const STRING_VIEW& delimited_ALT_str, const bool begins_before_curr_start_position);
    /*
     * Should be called after the allele info for all samples for the current position is updated
     * Adds spanning deletions and NON_REF allele if necessary to the merged alleles vector
     */
    void finished_updating_allele_info_for_current_position();
    /*
     * Sample/callset corresponding to row_query_idx is being invalidated. Make necessary state
     * updates
     */
    void remove_allele_info(const size_t row_query_idx);
    /*
     * Control flags
     */
    inline void set_contains_deletion(const size_t row_query_idx, const bool v) {
      assert(row_query_idx < m_contains_deletion.size());
      m_contains_deletion[row_query_idx] = v;
    }
    inline void set_contains_MNV(const size_t row_query_idx, const bool v) {
      assert(row_query_idx < m_contains_MNV.size());
      m_contains_MNV[row_query_idx] = v;
    }
    inline bool contains_deletion_or_MNV(const size_t row_query_idx) const {
      assert(row_query_idx < m_contains_deletion.size());
      assert(row_query_idx < m_contains_MNV.size());
      return m_contains_deletion[row_query_idx] || m_contains_MNV[row_query_idx];
    }
    inline bool contains_NON_REF_allele(const size_t row_query_idx) const {
      assert(row_query_idx < m_NON_REF_idx_vec.size());
      return (m_NON_REF_idx_vec[row_query_idx] != UNDEFINED_ATTRIBUTE_IDX_VALUE);
    }
    inline uint64_t get_num_calls_with_deletions_or_MNVs() const {
      return m_num_calls_with_deletions_or_MNVs;
    }
    inline uint64_t get_num_calls_with_NON_REF_allele() const {
      return m_num_calls_with_NON_REF_allele;
    }
    /*
     * The following two functions deal with samples that begin at the current location (current defined by the
     * GenomicsDBGVCFIterator) and have deletions/MNVs.
     * Used for computing spanning deletions using min PL value.
     */
    inline const std::vector<int64_t>& get_row_query_idx_with_deletions_or_MNVs_at_current_location_vec() const {
      return m_row_query_idx_with_deletion_MNV_at_current_location_vec;
    }
    //Ideally this should return std::span, but that's not yet part of the standard
    //See the data structure below for understanding what this function does
    inline std::pair<const unsigned*, size_t> get_deletion_MNV_indexes_for_row_query_idx_at_index(const size_t idx) const {
      assert(idx < m_row_query_idx_with_deletion_MNV_at_current_location_vec.size());
      assert(idx+1u < m_index_in_deletions_MNVs_vec.size());
      const auto vec_index = m_index_in_deletions_MNVs_vec[idx];
      assert(m_index_in_deletions_MNVs_vec[idx+1u] > vec_index); //only rows with deletions/MNVs must be entered into this structure
      return std::pair<const unsigned*, size_t>(&(m_deletion_MNV_allele_idx_vec[vec_index]),
          m_index_in_deletions_MNVs_vec[idx+1u]-vec_index);
    }
    /*
     * Utility function that returns merged REF and ALT
     */
    inline const std::vector<MergedAllelesVecEntry>& get_merged_normalized_REF_ALT_vec() const {
      return m_merged_alleles_vec;
    }
    /*
     * Utility function that returns a vector of STRING_VIEW containing the merged alleles as per the  
     * VCF spec. Alleles have the appropriate suffixes appended. The buffer string is where the alleles
     * and suffixes are copied to from GenomicsDBBuffer. Element 0 is REF.
     * Keep in mind that the STRING_VIEW objects in alleles_vec point to locations in the buffer. So, any
     * modifications of buffer (resize, append etc) will leave dangling pointers in alleles_vec.
     */
    void get_merged_VCF_spec_alleles_vec(std::string& buffer, std::vector<STRING_VIEW>& alleles_vec) const;
    inline const CombineAllelesLUT& get_merged_alleles_lut() const {
      return m_alleles_LUT;
    }
  private:
    std::pair<unsigned, bool> handle_single_base_ALT(const char base);
    std::pair<unsigned, bool> handle_insertion_or_symbolic_allele(const STRING_VIEW& ALT);
    void handle_allele(const size_t row_query_idx, const STRING_VIEW& REF, AlleleConfig& allele_config);
  private:
    unsigned m_spanning_deletion_allele_idx; //'*' index
    uint64_t m_num_calls_with_deletions_or_MNVs;
    uint64_t m_num_calls_with_NON_REF_allele;
    //One entry per queried row in all these vectors
    std::vector<bool> m_contains_deletion;
    std::vector<bool> m_contains_MNV;
    std::vector<bool> m_is_REF_block; //contains single ALT allele <NON_REF>
    //Index of NON_REF allele in input ALT vec, UNDEFINED_ATTRIBUTE_IDX_VALUE if NON_REF doesn't exist
    std::vector<unsigned> m_NON_REF_idx_vec;
    //Size equal to number of merged alleles (including REF)
    std::vector<MergedAllelesVecEntry> m_merged_alleles_vec;
    //The following structures are used to track unique alleles during merging. The unsigned value in the maps/table
    //is the index of the ALT allele in the merged ALT allele list
    //The most common ALT alleles appear to be SNVs and NON_REF (which is stored as char '&' in TileDB/GenomicsDB).
    //So, we can optimize the lookup for single base ALT alleles by keeping a a table of size 256 (ASCII) and using
    //the ALT base to index into it
    unsigned m_single_base_ALT_allele_to_merged_idx[256];
    //This map is used for normalized alt alleles whose REF is a single base - this includes insertions and symbolic
    //alleles. The key is the ALT allele.
    std::unordered_map<STRING_VIEW, unsigned> m_single_base_REF_ALT_allele_to_index;
    //Normalized deletion alleles can be uniquely identified by their REF length alone
    //Example TAGC,->T TA->T, TAG->T
    //So m_deletion_REF_length_to_index is simply a map<REF_length, index>. The REF is tracked in
    //m_merged_alleles_vec
    std::unordered_map<unsigned, unsigned> m_deletion_REF_length_to_index;
    //Used only for MNVs since both REF and ALT are multi-base strings (map<ALT, unsigned>)
    //The property is that for a normalized MNV, the last base is different from the REF. So, just use 
    //ALT in the key of the map and track REF in m_merged_alleles_vec
    //You can make a case that for MNVs, REF and ALT should differ in every position (else split the variant
    //across positions) - https://genome.sph.umich.edu/wiki/Variant_classification
    //However, looks like GATK can produce MNVs that may not differ in all positions
    std::unordered_map<STRING_VIEW, unsigned> m_MNV_ALT_to_index;
    //LUT mapping of input alleles to merged alleles for all samples
    CombineAllelesLUT m_alleles_LUT;
    //Data structures for dealing with spanning deletions/MNVs
    //A sample may have multiple deletions/MNVs at a location.
    //For spanning deletions/MNVs, we need to track which allele indexes are deletions/MNVs
    //so that we can find the allele with minimum PL value. This allele will be used as the
    //spanning allele in the subsequent positions.
    //The easiest structure to achieve this is a vector<vector<unsigned>> where the outer
    //vector contains 1 entry per sample while the inner vector contains the allele indexes
    //of deletions/MNVs for that sample. However, this leads to non-contiguous accesses and more
    //importantly dynamic re-allocations when the outer vector is cleared/resized.
    //In the spirit of columnar structures and avoiding allocations/frees, I flatten the vector of
    //vectors structure into a single vector and store the number of alleles per sample in a separate
    //vector with one element per sample. This is similar to how TileDB stores variable length fields
    //Also, we only store the samples which actually have a deletion/MNV - not all samples.
    //Contains row indexes of samples with new data at current iteration and containing deletions/MNVs. By tracking
    //these indexes, we can determine which rows must insert spanning deletions in the subsequent locations
    std::vector<int64_t> m_row_query_idx_with_deletion_MNV_at_current_location_vec;
    //The size of the vector below is m_row_query_idx_with_deletion_MNV_at_current_location_vec.size()+1
    //Each entry contains the index in m_deletion_MNV_allele_idx_vec where the deletion/MNV indexes
    //for this particular row/sample begin. The number of deletions/MNVs for each row can be computed
    //by subtracting current index from the next index
    std::vector<size_t> m_index_in_deletions_MNVs_vec;
    //Contains index of alleles which are deletions/MNVs in the input allele list
    std::vector<unsigned> m_deletion_MNV_allele_idx_vec;
};

#endif
