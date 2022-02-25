#
# FindGPerftools.cmake
#
# The MIT License
#
# Copyright (c) 2022 Omics Data Automation, Inc.
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
# Determine compiler flags for gperftools
# Once done this will define
# GPERFTOOLS_FOUND - gperftools found
#

if (USE_GPERFTOOLS)
  find_path(GPERFTOOLS_INCLUDE_DIR NAMES gperftools/profiler.h HINTS "${GPERFTOOLS_DIR}" "${GPERFTOOLS_DIR}/include")
  find_library(GPERFTOOLS_PROFILER_LIBRARY NAMES profiler HINTS "${GPERFTOOLS_DIR}" "${GPERFTOOLS_DIR}/lib64" "${GPERFTOOLS_DIR}/lib")
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GPerftools "Could not find gperftools headers and/or libraries ${DEFAULT_MSG}" GPERFTOOLS_PROFILER_LIBRARY GPERFTOOLS_INCLUDE_DIR)
endif()

if (USE_GPERFTOOLS_HEAP)
  find_path(GPERFTOOLS_INCLUDE_DIR NAMES gperftools/tcmalloc.h HINTS "${GPERFTOOLS_DIR}" "${GPERFTOOLS_DIR}/include")
  find_library(GPERFTOOLS_HEAP_PROFILER_LIBRARY NAMES tcmalloc HINTS "${GPERFTOOLS_DIR}" "${GPERFTOOLS_DIR}/lib64" "${GPERFTOOLS_DIR}/lib")
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GPerftools "Could not find gperftools heap headers and/or libraries ${DEFAULT_MSG}" GPERFTOOLS_HEAP_PROFILER_LIBRARY GPERFTOOLS_INCLUDE_DIR)
endif()
