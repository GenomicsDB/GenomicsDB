/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
 * Copyright (c) 2019 Omics Data Automation, Inc.
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

#include <catch2/catch.hpp>

#include "vid_mapper.h"
#include "genomicsdb_multid_vector_field.h"
#include "variant_operations.h"

using Catch::Equals;

void initialize_LUT_reordered(CombineAllelesLUT& lut);

TEST_CASE("multid_vector 2D_test", "[2D_test]")
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(FieldElementTypeDescriptor(std::type_index(typeid(int)), BCF_HT_INT));
  field_info.set_vcf_type(FieldElementTypeDescriptor(std::type_index(typeid(char)), BCF_HT_STR));
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(2u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_A);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "|");
  length_descriptor.set_vcf_delimiter(1u, ",");
  field_info.modify_field_type_if_multi_dim_field();
  GenomicsDBMultiDVectorField two_d_vector(field_info);
  std::string value = "1,2|3,4,5";
  auto total_size = two_d_vector.parse_and_store_numeric(value.c_str(), value.length());
  //Size of dim 0 in first 8 bytes
  auto dim_0_data_size = *(reinterpret_cast<const uint64_t*>(&(two_d_vector.get_rw_data(0u)[0u])));
  CHECK(total_size[0u] == dim_0_data_size
      +sizeof(uint64_t) //size - 8 bytes
      +sizeof(uint64_t) //#entries - 8 bytes
      +3*sizeof(uint64_t)); //3 offsets
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data(0u)[0]), &field_info);
  CHECK(dim_0_idx.get_current_dim_index() == -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_size_of_current_index() == 8u);
  auto data_ptr = dim_0_idx.get_ptr<int>();
  CHECK(data_ptr[0u] == 1);
  CHECK(data_ptr[1u] == 2);
  {
    //A[0][i]
    auto dim_1_idx = dim_0_idx;
    //A[0][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 0u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 2u);
    CHECK(dim_1_idx.get_size_of_current_index() == sizeof(int));
    CHECK(dim_1_idx.get_element<int>() == 1);
    //A[0][1]
    dim_1_idx.advance_index_in_current_dimension();
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 1u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 2u);
    CHECK(dim_1_idx.get_size_of_current_index() == sizeof(int));
    CHECK(dim_1_idx.get_element<int>() == 2);
  }
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_size_of_current_index() == 12u);
  data_ptr = dim_0_idx.get_ptr<int>();
  CHECK(data_ptr[0u] == 3);
  CHECK(data_ptr[1u] == 4);
  CHECK(data_ptr[2u] == 5);
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  two_d_vector.run_operation(print_op, std::vector<const uint8_t*>(1u, &(two_d_vector.get_rw_data(0u)[0])));
  CHECK(fptr.str() == value);
  {
    //Test re-ordering
    CombineAllelesLUT tmp_lut;
    ////A_inp mapped to A_merged
    ////B_inp mapped to C_merged
    ////C_inp [NON_REF] mapped to D_merged
    ////Implies that any allele corresponding to B_merged will get mapped to NON_REF allele
    ///[ 1, 2 ] corresponds to B_inp, [3, 4, 5] corresponds to C_inp (NON_REF)
    ///B_merged = C_inp = [3,4,5] (NON_REF)
    ///C_merged = B_inp = [1,2]
    ///D_merged = C_inp = [3,4,5] (NON_REF)
    ///A is ignored since length descriptor is BCF_VL_A (ALT alleles only)
    initialize_LUT_reordered(tmp_lut);
    std::vector<uint8_t> reordered_data;
    remap_allele_specific_annotations(two_d_vector.get_rw_data(0u), reordered_data,
        0u, tmp_lut, 4u, true, 2u,
        field_info);
    //Size of dim 0 in first 8 bytes
    auto dim_0_data_size = *(reinterpret_cast<const uint64_t*>(&(reordered_data[0u])));
    CHECK((total_size[0u]
	   +3*sizeof(int) //the extra 3 integers because B_merged gets mapped to C_inp (NON_REF) and the [3, 4, 5] is duplicated
           +sizeof(uint64_t)) == //one extra entry so one more offset value
          (dim_0_data_size
	   +sizeof(uint64_t) //size - 8 bytes
	   +sizeof(uint64_t) //#entries - 8 bytes
           +4*sizeof(uint64_t))); //4 offsets
    //Index
    GenomicsDBMultiDVectorIdx dim_0_idx(&(reordered_data[0]), &field_info);
    CHECK(dim_0_idx.get_current_dim_index() == -1);
    //A[0] - B_merged == C_inp (NON_REF)
    dim_0_idx.advance_to_index_in_next_dimension(0u);
    CHECK(dim_0_idx.get_current_dim_index() == 0);
    CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
    CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 3u);
    CHECK(dim_0_idx.get_size_of_current_index() == 3*sizeof(int)); //[3, 4, 5]
    auto data_ptr = dim_0_idx.get_ptr<int>();
    CHECK(data_ptr[0u] == 3);
    CHECK(data_ptr[1u] == 4);
    CHECK(data_ptr[2u] == 5);
    //A[1] - C_merged == B_inp
    dim_0_idx.advance_index_in_current_dimension();
    CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
    CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 3u);
    CHECK(dim_0_idx.get_size_of_current_index() == 2*sizeof(int)); //[1, 2]
    data_ptr = dim_0_idx.get_ptr<int>();
    CHECK(data_ptr[0u] == 1);
    CHECK(data_ptr[1u] == 2);
    //A[2] - D_merged == C_inp (NON_REF)
    dim_0_idx.advance_index_in_current_dimension();
    CHECK(dim_0_idx.get_current_index_in_current_dimension() == 2u);
    CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 3u);
    CHECK(dim_0_idx.get_size_of_current_index() == 3*sizeof(int)); //[3, 4, 5]
    data_ptr = dim_0_idx.get_ptr<int>();
    CHECK(data_ptr[0u] == 3);
    CHECK(data_ptr[1u] == 4);
    CHECK(data_ptr[2u] == 5);
  }
}

TEST_CASE("multid_vector 3D_test", "[3D_test]")
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(FieldElementTypeDescriptor(std::type_index(typeid(int)), BCF_HT_INT));
  field_info.set_vcf_type(FieldElementTypeDescriptor(std::type_index(typeid(char)), BCF_HT_STR));
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(3u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(2u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "$");
  length_descriptor.set_vcf_delimiter(1u, "|");
  length_descriptor.set_vcf_delimiter(2u, ",");
  field_info.modify_field_type_if_multi_dim_field();
  GenomicsDBMultiDVectorField three_d_vector(field_info);
  std::string value = "1,2|3,4,5$6,7,8,9|10,11|12";
  auto total_size = three_d_vector.parse_and_store_numeric(value.c_str(), value.length());
  //Size of dim 0 data in first 8 bytes
  auto dim_0_data_size = *(reinterpret_cast<const uint64_t*>(&(three_d_vector.get_rw_data(0u)[0u])));
  CHECK(total_size[0u] == dim_0_data_size
      +sizeof(uint64_t) //size - 8 bytes
      +sizeof(uint64_t) //#entries - 8 bytes
      +3*sizeof(uint64_t)); //3 offsets
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(three_d_vector.get_rw_data(0u)[0]), &field_info);
  CHECK(dim_0_idx.get_current_dim_index() == -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_size_of_current_index() == 60u);
  {
    //A[0][i]
    auto dim_1_idx = dim_0_idx;
    //A[0][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 0u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 2u);
    CHECK(dim_1_idx.get_size_of_current_index() == 8u);
    auto data_ptr = dim_1_idx.get_ptr<int>();
    CHECK(data_ptr[0] == 1);
    CHECK(data_ptr[1] == 2);
    //A[0][1]
    dim_1_idx.advance_index_in_current_dimension();
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 1u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 2u);
    CHECK(dim_1_idx.get_size_of_current_index() == 12u);
    data_ptr = dim_1_idx.get_ptr<int>();
    CHECK(data_ptr[0] == 3);
    CHECK(data_ptr[1] == 4);
    CHECK(data_ptr[2] == 5);
    {
      //A[0][1][i]
      auto dim_2_idx = dim_1_idx;
      dim_2_idx.advance_to_index_in_next_dimension(0u);
      //A[0][1][0]
      CHECK(dim_2_idx.get_current_dim_index() == 2);
      CHECK(dim_2_idx.get_current_index_in_current_dimension() == 0u);
      CHECK(dim_2_idx.get_num_entries_in_current_dimension() == 3u);
      CHECK(dim_2_idx.get_size_of_current_index() == sizeof(int));
      CHECK(dim_2_idx.get_element<int>() == 3);
      //A[0][1][1]
      dim_2_idx.advance_index_in_current_dimension();
      CHECK(dim_2_idx.get_current_dim_index() == 2);
      CHECK(dim_2_idx.get_current_index_in_current_dimension() == 1u);
      CHECK(dim_2_idx.get_num_entries_in_current_dimension() == 3u);
      CHECK(dim_2_idx.get_size_of_current_index() == sizeof(int));
      CHECK(dim_2_idx.get_element<int>() == 4);
      //A[0][1][2]
      dim_2_idx.advance_index_in_current_dimension();
      CHECK(dim_2_idx.get_current_dim_index() == 2);
      CHECK(dim_2_idx.get_current_index_in_current_dimension() == 2u);
      CHECK(dim_2_idx.get_num_entries_in_current_dimension() == 3u);
      CHECK(dim_2_idx.get_size_of_current_index() == sizeof(int));
      CHECK(dim_2_idx.get_element<int>() == 5);
    }
  }
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_size_of_current_index() == 76u);
  {
    //A[1][i]
    auto dim_1_idx = dim_0_idx;
    //A[1][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 0u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 3u);
    CHECK(dim_1_idx.get_size_of_current_index() == 16u);
    auto data_ptr = dim_1_idx.get_ptr<int>();
    CHECK(data_ptr[0] == 6);
    CHECK(data_ptr[1] == 7);
    CHECK(data_ptr[2]== 8);
    CHECK(data_ptr[3] == 9);
    //A[0][1]
    dim_1_idx.advance_index_in_current_dimension();
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 1u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 3u);
    CHECK(dim_1_idx.get_size_of_current_index() == 8u);
    data_ptr = dim_1_idx.get_ptr<int>();
    CHECK(data_ptr[0] == 10);
    CHECK(data_ptr[1] == 11);
    //A[0][2]
    dim_1_idx.advance_index_in_current_dimension();
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 2u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 3u);
    CHECK(dim_1_idx.get_size_of_current_index() == 4u);
    data_ptr = dim_1_idx.get_ptr<int>();
    CHECK(data_ptr[0] == 12);
  }
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  three_d_vector.run_operation(print_op, std::vector<const uint8_t*>(1u, &(three_d_vector.get_rw_data(0u)[0])));
  CHECK(fptr.str() == value);
}

TEST_CASE("multid_vector 2D_test_with_missing_values", "[2D_test_with_missing_values]")
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(FieldElementTypeDescriptor(std::type_index(typeid(int)), BCF_HT_INT));
  field_info.set_vcf_type(FieldElementTypeDescriptor(std::type_index(typeid(char)), BCF_HT_STR));
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(2u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "|");
  length_descriptor.set_vcf_delimiter(1u, ",");
  field_info.modify_field_type_if_multi_dim_field();
  GenomicsDBMultiDVectorField two_d_vector(field_info);
  std::string value = "|3,,NaN";
  auto total_size = two_d_vector.parse_and_store_numeric(value.c_str(), value.length());
  //Size of dim 0 data in first 8 bytes
  auto dim_0_data_size = *(reinterpret_cast<const uint64_t*>(&(two_d_vector.get_rw_data(0u)[0u])));
  CHECK(total_size[0u] == dim_0_data_size
      +sizeof(uint64_t) //size - 8 bytes
      +sizeof(uint64_t) //#entries - 8 bytes
      +3*sizeof(uint64_t)); //3 offsets
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data(0u)[0]), &field_info);
  CHECK(dim_0_idx.get_current_dim_index() == -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_size_of_current_index() == sizeof(float));
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_size_of_current_index() == 3*sizeof(int));
  auto data_ptr = dim_0_idx.get_ptr<int>();
  CHECK(data_ptr[0u] == 3);
  CHECK(data_ptr[1u] == get_bcf_missing_value<int>());
  CHECK(data_ptr[2u] == get_bcf_missing_value<int>());
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  two_d_vector.run_operation(print_op, std::vector<const uint8_t*>(1u, &(two_d_vector.get_rw_data(0u)[0])));
  std::string serialization_return_value = "|3,,"; //NaN gets replaced with missing
  CHECK_THAT(fptr.str(),  Equals(serialization_return_value));
}

TEST_CASE("multid_vector 2D_float_test", "[2D_float_test]")
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(FieldElementTypeDescriptor(std::type_index(typeid(float)), BCF_HT_REAL));
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(2u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "|");
  length_descriptor.set_vcf_delimiter(1u, ",");
  field_info.modify_field_type_if_multi_dim_field();
  GenomicsDBMultiDVectorField two_d_vector(field_info);
  std::string value = "29659.00|1.00|2.00|3.00";
  auto total_size = two_d_vector.parse_and_store_numeric(value.c_str(), value.length());
  //Size of dim 0 in first 8 bytes
  auto dim_0_data_size = *(reinterpret_cast<const uint64_t*>(&(two_d_vector.get_rw_data(0u)[0u])));
  CHECK(total_size[0u] == dim_0_data_size
      +sizeof(uint64_t) //size - 8 bytes
      +sizeof(uint64_t) //#entries - 8 bytes
      +5*sizeof(uint64_t)); //5 offsets
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data(0u)[0]), &field_info);
  CHECK(dim_0_idx.get_current_dim_index() == -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 4u);
  CHECK(dim_0_idx.get_size_of_current_index() == sizeof(float));
  auto data_ptr = dim_0_idx.get_ptr<float>();
  //  EXPECT_NEAR(data_ptr[0u] == 29659.0f == 1e-5f);
  CHECK(data_ptr[0u] == Approx(29659.0f).margin(1e-5f)); 
  {
    //A[0][i]
    auto dim_1_idx = dim_0_idx;
    //A[0][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    CHECK(dim_1_idx.get_current_dim_index() == 1);
    CHECK(dim_1_idx.get_current_index_in_current_dimension() == 0u);
    CHECK(dim_1_idx.get_num_entries_in_current_dimension() == 1u);
    CHECK(dim_1_idx.get_size_of_current_index() == sizeof(float));
    CHECK(dim_1_idx.get_element<float>() == 29659.0);
  }
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 4u);
  CHECK(dim_0_idx.get_size_of_current_index() == sizeof(float));
  data_ptr = dim_0_idx.get_ptr<float>();
  CHECK(data_ptr[0u] == 1.0);
  //A[2]
  dim_0_idx.advance_index_in_current_dimension();
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 2u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 4u);
  CHECK(dim_0_idx.get_size_of_current_index() == sizeof(float));
  data_ptr = dim_0_idx.get_ptr<float>();
  CHECK(data_ptr[0u] == 2.0);
  //A[3]
  dim_0_idx.advance_index_in_current_dimension();
  CHECK(dim_0_idx.get_current_dim_index() == 0);
  CHECK(dim_0_idx.get_current_index_in_current_dimension() == 3u);
  CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 4u);
  CHECK(dim_0_idx.get_size_of_current_index() == sizeof(float));
  data_ptr = dim_0_idx.get_ptr<float>();
  CHECK(data_ptr[0u] == 3.0);
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  two_d_vector.run_operation(print_op, std::vector<const uint8_t*>(1u, &(two_d_vector.get_rw_data(0u)[0])));
  //EXPECT_THAT(fptr.str(), MatchesRegex("29659\\.?.*|\\s*1\\.?.*|\\s*2\\.?.*|\\s*3\\.?.*"));
  //CHECK(fptr.str(), "29659|1|2|3");
}

TEST_CASE("multid_vector 2D_tuple_test", "[2D_tuple_test]")
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  FieldElementTypeDescriptor genomicsdb_type(2u);
  genomicsdb_type.set_tuple_element_type(0u, std::type_index(typeid(int)), BCF_HT_INT);
  genomicsdb_type.set_tuple_element_type(1u, std::type_index(typeid(float)), BCF_HT_REAL);
  field_info.set_type(genomicsdb_type);
  field_info.set_vcf_type(FieldElementTypeDescriptor(std::type_index(typeid(char)), BCF_HT_STR));
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(2u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "|");
  length_descriptor.set_vcf_delimiter(1u, ",");
  field_info.modify_field_type_if_multi_dim_field();
  GenomicsDBMultiDVectorField two_d_vector(field_info);
  std::string value = "1,1.5,2,2.5|3,3.5,4,4.5,5,5.5";
  for(auto k=0u;k<2u;++k) //once from value and once from print at the end of the first iteration
  {
    auto total_size = two_d_vector.parse_and_store_numeric(value.c_str(), value.length());
    //Size of dim 0 in first 8 bytes
    auto dim_0_data_size = *(reinterpret_cast<const uint64_t*>(&(two_d_vector.get_rw_data(0u)[0u])));
    CHECK(total_size[0u] == dim_0_data_size
        +sizeof(uint64_t) //size - 8 bytes
        +sizeof(uint64_t) //#entries - 8 bytes
        +3*sizeof(uint64_t)); //3 offsets
    //Integer element of the tuple
    {
      //Index
      GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data(0u)[0]), &field_info);
      CHECK(dim_0_idx.get_current_dim_index() == -1);
      //A[0]
      dim_0_idx.advance_to_index_in_next_dimension(0u);
      CHECK(dim_0_idx.get_current_dim_index() == 0);
      CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
      CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
      CHECK(dim_0_idx.get_size_of_current_index() == 8u);
      auto data_ptr = dim_0_idx.get_ptr<int>();
      CHECK(data_ptr[0u] == 1);
      CHECK(data_ptr[1u] == 2);
      //A[1]
      dim_0_idx.advance_index_in_current_dimension();
      CHECK(dim_0_idx.get_current_dim_index() == 0);
      CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
      CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
      CHECK(dim_0_idx.get_size_of_current_index() == 12u);
      data_ptr = dim_0_idx.get_ptr<int>();
      CHECK(data_ptr[0u] == 3);
      CHECK(data_ptr[1u] == 4);
      CHECK(data_ptr[2u] == 5);
      //VCF printer
      std::ostringstream fptr;
      FieldInfo int_tuple_element_field_info;
      int_tuple_element_field_info.set_type(FieldElementTypeDescriptor(std::type_index(typeid(int)), BCF_HT_INT));
      int_tuple_element_field_info.m_length_descriptor.resize(2u);
      int_tuple_element_field_info.m_length_descriptor.set_length_descriptor(0u, BCF_VL_R);
      int_tuple_element_field_info.m_length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
      //Delimiters
      int_tuple_element_field_info.m_length_descriptor.set_vcf_delimiter(0u, "|");
      int_tuple_element_field_info.m_length_descriptor.set_vcf_delimiter(1u, ",");
      int_tuple_element_field_info.modify_field_type_if_multi_dim_field();
      GenomicsDBMultiDVectorField int_tuple_element_vector(int_tuple_element_field_info);
      GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, int_tuple_element_field_info);
      int_tuple_element_vector.run_operation(print_op, std::vector<const uint8_t*>(1u, &(two_d_vector.get_rw_data(0u)[0])));
      CHECK(fptr.str() == "1,2|3,4,5");
    }
    //FP element of the tuple
    {
      //Index
      GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data(1u)[0]), &field_info);
      CHECK(dim_0_idx.get_current_dim_index() == -1);
      //A[0]
      dim_0_idx.advance_to_index_in_next_dimension(0u);
      CHECK(dim_0_idx.get_current_dim_index() == 0);
      CHECK(dim_0_idx.get_current_index_in_current_dimension() == 0u);
      CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
      CHECK(dim_0_idx.get_size_of_current_index() == 8u);
      auto data_ptr = dim_0_idx.get_ptr<float>();
      CHECK(data_ptr[0u] == 1.5);
      CHECK(data_ptr[1u] == 2.5);
      //A[1]
      dim_0_idx.advance_index_in_current_dimension();
      CHECK(dim_0_idx.get_current_dim_index() == 0);
      CHECK(dim_0_idx.get_current_index_in_current_dimension() == 1u);
      CHECK(dim_0_idx.get_num_entries_in_current_dimension() == 2u);
      CHECK(dim_0_idx.get_size_of_current_index() == 12u);
      data_ptr = dim_0_idx.get_ptr<float>();
      CHECK(data_ptr[0u] == 3.5);
      CHECK(data_ptr[1u] == 4.5);
      CHECK(data_ptr[2u] == 5.5);
    }
    {
      std::ostringstream fptr;
      GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
      std::vector<const uint8_t*> ptr_vec(2u);
      ptr_vec[0] = &(two_d_vector.get_rw_data(0u)[0]);
      ptr_vec[1] = &(two_d_vector.get_rw_data(1u)[0]);
      two_d_vector.run_operation(print_op, ptr_vec);
      value = fptr.str();
    }
  }
}
