#
# cmake Module for libuuid
#
# The MIT License
#
# Copyright (c) 2023 dātma, inc™
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

# Determine compiler flags for libuuid
# Once done this will define
# LIBUUID_FOUND - libuuid found

find_path(LIBUUID_INCLUDE_DIR NAMES uuid/uuid.h HINTS "${LIBUUID_DIR}/include" "${LIBUUID_DIR}")

if(APPLE)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(libuuid "Could not find libuuid headers ${DEFAULT_MSG}" LIBUUID_INCLUDE_DIR)
else()
  if(BUILD_DISTRIBUTABLE_LIBRARY)
    find_library(LIBUUID_LIBRARY NAMES libuuid.a uuid HINTS "${LIBUUID_DIR}/lib64" "${LIBUUID_DIR}/lib" "${LIBUUID_DIR}")
  else()
    find_library(LIBUUID_LIBRARY NAMES uuid HINTS "${LIBUUID_DIR}/lib64" "${LIBUUID_DIR}/lib" "${LIBUUID_DIR}")
  endif()
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(libuuid "Could not find libuuid headers and/or libraries ${DEFAULT_MSG}" LIBUUID_INCLUDE_DIR LIBUUID_LIBRARY)
endif()
