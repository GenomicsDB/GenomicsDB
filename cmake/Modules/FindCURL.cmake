#
# cmake/Modules/FindCURL.cmake
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
# Find the native CURL headers and libraries.
#
#  CURL_INCLUDE_DIRS   - where to find curl/curl.h, etc.
#  CURL_LIBRARIES      - List of libraries when using curl.
#  CURL_FOUND          - True if curl found.
#  CURL_VERSION_STRING - the version of curl found (since CMake 2.8.8)

# Look for the header file.
find_path(CURL_INCLUDE_DIR NAMES curl/curl.h HINTS ${CURL_PREFIX_DIR} ${CURL_PREFIX_DIR}/include)
mark_as_advanced(CURL_INCLUDE_DIR)

if(BUILD_DISTRIBUTABLE_LIBRARY AND NOT APPLE)
  find_library(CURL_LIBRARY NAMES libcurl.a
    HINTS ${CURL_PREFIX_DIR} ${CURL_PREFIX_DIR}/lib64 ${CURL_PREFIX_DIR}/lib)
else()
  find_library(CURL_LIBRARY NAMES curl libcurl
    HINTS ${CURL_PREFIX_DIR} ${CURL_PREFIX_DIR}/lib64 ${CURL_PREFIX_DIR}/lib
    )
endif()
mark_as_advanced(CURL_LIBRARY)

if(CURL_INCLUDE_DIR)
  foreach(_curl_version_header curlver.h curl.h)
    if(EXISTS "${CURL_INCLUDE_DIR}/curl/${_curl_version_header}")
      file(STRINGS "${CURL_INCLUDE_DIR}/curl/${_curl_version_header}" curl_version_str REGEX "^#define[\t ]+LIBCURL_VERSION[\t ]+\".*\"")

      string(REGEX REPLACE "^#define[\t ]+LIBCURL_VERSION[\t ]+\"([^\"]*)\".*" "\\1" CURL_VERSION_STRING "${curl_version_str}")
      unset(curl_version_str)
      break()
    endif()
  endforeach()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CURL
                                  REQUIRED_VARS CURL_LIBRARY CURL_INCLUDE_DIR
                                  VERSION_VAR CURL_VERSION_STRING)

if(CURL_FOUND)
  set(CURL_LIBRARIES ${CURL_LIBRARY})
  set(CURL_INCLUDE_DIRS ${CURL_INCLUDE_DIR})
endif()
