#
# src/test/cpp/CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2019 Omics Data Automation, Inc.
# Copyright (c) 2023-2024 da̅tma, inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

include_directories(Catch INTERFACE Catch2/single_include)
include_directories(${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS})

set(CPP_TEST_SOURCES
  src/ctest_main.cc
  src/test_config.cc
  src/test_non_diploid_mapper.cc
  src/test_multid_vector.cc
  src/test_pb.cc
  src/test_logger.cc
  src/test_columnar_gvcf_iterator.cc
  src/test_htslib_plugin.cc
  src/test_alleles_combiner.cc
  src/test_gt_remapper.cc
  src/test_bgen.cc
  src/test_query_variants.cc
  src/test_mem_utils.cc
  )

if(LIBDBI_FOUND)
  set(CPP_TEST_SOURCES
    ${CPP_TEST_SOURCES}
    src/test_mapping_data_loader.cc)
endif()

function(target_link_libraries_for_GenomicsDB_tests target)
   build_GenomicsDB_links(${target})
   # Need to explicitly add in pthread/dl libs for static linkage
   if(BUILD_DISTRIBUTABLE_LIBRARY)
     target_link_libraries(${target} pthread ${CMAKE_DL_LIBS})
   endif()
endfunction()

add_executable(ctests ${CPP_TEST_SOURCES})
target_link_libraries_for_GenomicsDB_tests(ctests)

#test_bgen has to be fixed, until then...
if(BUILD_DISTRIBUTABLE_LIBRARY OR BUILD_FOR_PYTHON OR DISABLE_MPI)
  set(INCLUDE_TEST_BGEN 0)
else()
  set(INCLUDE_TEST_BGEN 1)
endif()
target_compile_definitions(ctests PRIVATE -DINCLUDE_TEST_BGEN=${INCLUDE_TEST_BGEN} -DGENOMICSDB_CTESTS_DIR="${CMAKE_CURRENT_BINARY_DIR}/inputs/" -DGENOMICSDB_TESTS_SRC_DIR="${CMAKE_SOURCE_DIR}/tests/")
add_test(ctests ctests -d yes)

add_executable(api_tests
  src/ctest_main.cc
  src/test_genomicsdb_api.cc)
target_link_libraries_for_GenomicsDB_tests(api_tests)
target_compile_definitions(api_tests PRIVATE -DGENOMICSDB_CTESTS_DIR="${CMAKE_CURRENT_BINARY_DIR}/inputs/")
add_test(api_tests api_tests -d yes)

add_custom_command(TARGET api_tests
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/../inputs/test.tgz
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../inputs/test.tgz)

add_custom_target(all_ctests COMMAND ${CMAKE_CTEST_COMMAND} -V DEPENDS ctests api_tests)

build_GenomicsDB_executable(test_genomicsdb_demo)
