#
# src/test/cpp/CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2019 Omics Data Automation, Inc.
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
  src/test_non_diploid_mapper.cc
  src/test_multid_vector.cc
  src/test_pb.cc
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
add_test(ctests ctests -d yes)

add_executable(api_tests
  src/ctest_main.cc
  src/test_genomicsdb_api.cc)
target_link_libraries_for_GenomicsDB_tests(api_tests)
target_compile_definitions(api_tests PRIVATE -DGENOMICSDB_CTESTS_DIR="${CMAKE_CURRENT_BINARY_DIR}/inputs/")
add_test(api_tests api_tests -d yes)

add_custom_command(TARGET api_tests
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/../inputs/test.tgz
  WORKING_DRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../inputs/test.tgz)

add_custom_target(all_ctests COMMAND ${CMAKE_CTEST_COMMAND} -V DEPENDS ctests api_tests)
