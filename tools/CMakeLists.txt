#
# tools/CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2019-2023 Omics Data Automation, Inc.
# Copyright (c) 2023 da̅tma, inc.
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

include_directories("include")

# CMake/GCC doesn't add "-pie" by default for executables (CMake issue #14983)
if(NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pie")
endif()

if(NOT BUILD_DISTRIBUTABLE_LIBRARY)
  build_GenomicsDB_executable(create_genomicsdb_workspace)
  build_GenomicsDB_executable(vcf2genomicsdb_init)
  build_GenomicsDB_executable(vcf2genomicsdb)
  build_GenomicsDB_executable(consolidate_genomicsdb_array)
  if(NOT DISABLE_MPI)
    build_GenomicsDB_executable(gt_mpi_gather)
    build_GenomicsDB_executable(vcfdiff)
    build_GenomicsDB_executable(vcf_histogram)
  endif()
  add_test(NAME tools_tests COMMAND ${CMAKE_SOURCE_DIR}/tests/test_tools.sh ${CMAKE_SOURCE_DIR}/tests/inputs/vcfs ${CMAKE_INSTALL_PREFIX} ${DISABLE_MPI})
endif()
