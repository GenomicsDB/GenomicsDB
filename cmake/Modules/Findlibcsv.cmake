#
# Findlibcsv.cmake
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
# Determine compiler flags for libcsv
# Once done this will define
# LIBCSV_FOUND - libcsv found

find_path(LIBCSV_INCLUDE_DIR NAMES csv.h HINTS "${LIBCSV_DIR}/include" "${LIBCSV_DIR}")
find_library(LIBCSV_LIBRARY NAMES csv HINTS "${LIBCSV_DIR}/lib" "${LIBCSV_DIR}/.libs" "${LIBCSV_DIR}")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libcsv "Could not find libcsv headers and/or libraries ${DEFAULT_MSG}" LIBCSV_LIBRARY LIBCSV_INCLUDE_DIR)
