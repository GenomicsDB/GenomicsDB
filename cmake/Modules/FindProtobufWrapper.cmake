#
# CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2019-2021 Omics Data Automation, Inc.
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

add_custom_target(protobuf_ep)

message(CHECK_START "Try finding Protobuf Library")

if(PROTOBUF_ROOT_DIR)
  set(PROTOBUF_PREFIX ${PROTOBUF_ROOT_DIR})
elseif(DEFINED ENV{PROTOBUF_ROOT_DIR})
  set(PROTOBUF_PREFIX $ENV{PROTOBUF_ROOT_DIR})
endif()
set(PROTOBUF_BIN_DIR ${PROTOBUF_PREFIX}/bin)
set(PROTOBUF_LIB_DIR ${PROTOBUF_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(PROTOBUF_INCLUDE_DIR ${PROTOBUF_PREFIX}/include)

if(EXISTS ${PROTOBUF_BIN_DIR}/protoc AND EXISTS ${PROTOBUF_LIB_DIR} AND EXISTS ${PROTOBUF_LIB_DIR} AND EXISTS ${PROTOBUF_INCLUDE_DIR})
  find_package(Protobuf PATHS ${PROTOBUF_PREFIX})
endif()

if(NOT Protobuf_FOUND AND NOT PROTOBUF_PREFIX)
  message(CHECK_FAIL "not found")
  message(FATAL_ERROR "Try invoking cmake with -DPROTOBUF_ROOT_DIR=/path/to/protobuf or -DCMAKE_INSTALL_PREFIX=/path/to/protobuf or set environment variable PROTOBUF_ROOT_DIR before invoking cmake")
elseif(Protobuf_FOUND)
  message(CHECK_PASS "found ${Protobuf_VERSION}")
else()
   # Try building from source
  message(CHECK_PASS "not found, building Protobuf ${GENOMICSDB_PROTOBUF_VERSION} as an external project")
  include(ExternalProject)
  ExternalProject_Add(protobuf_build
    PREFIX ${PROTOBUF_PREFIX}
    URL ${PROTOBUF_URL}
    SOURCE_SUBDIR cmake
    CMAKE_ARGS
    -Dprotobuf_BUILD_TESTS=OFF
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=${PROTOBUF_PREFIX}
    )
  add_dependencies(protobuf_ep protobuf_build)
endif()

set(PROTOBUF_PROTOC_EXECUTABLE ${PROTOBUF_BIN_DIR}/protoc)
set(PROTOBUF_INCLUDE_DIRS ${PROTOBUF_INCLUDE_DIR})
set(PROTOBUF_LIBRARIES ${PROTOBUF_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}protobuf${CMAKE_STATIC_LIBRARY_SUFFIX})
