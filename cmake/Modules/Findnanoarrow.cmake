#
# CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2024 dātma, inc™
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

add_custom_target(nanoarrow_ep)

message(CHECK_START "Finding nanarrow Build and Library")

include(GNUInstallDirs)

if(NANOARROW_ROOT_DIR)
  set(NANOARROW_PREFIX "${NANOARROW_ROOT_DIR}")
else()
  message(FATAL_ERROR "No NANOARROW_ROOT_DIR specified to cmake")
endif()
set(NANOARROW_INCLUDE_DIR "${NANOARROW_PREFIX}/include")
set(NANOARROW_LIB_DIR "${NANOARROW_PREFIX}/lib")
set(NANOARROW_LIBRARY "${NANOARROW_LIB_DIR}/libnanoarrow.a")

# Workaround for issues with semicolon separated substrings to ExternalProject_Add
# https://discourse.cmake.org/t/how-to-pass-cmake-osx-architectures-to-externalproject-add/2262
if(CMAKE_OSX_ARCHITECTURES)
  if(CMAKE_OSX_ARCHITECTURES MATCHES "x86_64" AND CMAKE_OSX_ARCHITECTURES MATCHES "arm64")
    set(CMAKE_OSX_ARGS -DCMAKE_OSX_ARCHITECTURES=arm64$<SEMICOLON>x86_64)
  else()
    set(CMAKE_OSX_ARGS -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES})
  endif()
endif()

if(EXISTS ${NANOARROW_INCLUDE_DIR} AND EXISTS ${NANOARROW_LIB_DIR} AND EXISTS ${NANOARROW_LIBRARY})
  message(STATUS "Found nanoarrow headers: ${NANOARROW_INCLUDE_DIR}")
  message(STATUS "Found nanoarrow library: ${NANOARROW_LIBRARY}")
  message(CHECK_PASS "found ${NANOARROW_ROOT_DIR}")
else()
  message(CHECK_PASS "not found, building nanoarrow ${NANOARROW_VERSION} as an external project")
  include(ExternalProject)
  ExternalProject_Add(nanoarrow_build
    PREFIX ${NANOARROW_PREFIX}
    URL ${NANOARROW_URL}
    CMAKE_ARGS ${CMAKE_OSX_ARGS}
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=${NANOARROW_PREFIX}
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DNANOARROW_ARROW_STATIC=ON
    -DCMAKE_C_VISIBILITY_PRESET=hidden
    -DCMAKE_CXX_VISIBILITY_PRESET=hidden
    -DCMAKE_VISIBILITY_INLINES_HIDDEN=TRUE
    )
  add_dependencies(nanoarrow_ep nanoarrow_build)
endif()

