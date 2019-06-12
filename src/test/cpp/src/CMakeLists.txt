#
# src/test/cpp/src/CMakeLists.txt
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

include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/src/test/cpp/Catch2/single_include)

set(CPP_TEST_SOURCES
  ctest_main.cc
  test_non_diploid_mapper.cc
  test_multid_vector.cc)

if(LIBDBI_FOUND)
  set(CPP_TEST_SOURCES
    ${CPP_TEST_SOURCES}
    test_mapping_data_loader.cc)
endif()

add_executable(ctests ${CPP_TEST_SOURCES})
build_GENOMICSDB_executable_common(ctests)
add_test(ctests ctests)

add_custom_target(all_ctests COMMAND ${CMAKE_CTEST_COMMAND} -V DEPENDS ctests)