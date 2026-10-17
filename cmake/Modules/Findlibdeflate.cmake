#
# Findlibdeflate.cmake
#
# The MIT License
#
# Copyright (c) 2026 Tim Fennell
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
# Finds the static libdeflate in LIBDEFLATE_ROOT_DIR, building it there first if it is missing. htslib
# is built with it so that BGZF blocks are compressed and decompressed by libdeflate rather than by
# the system's zlib, which on many Linux hosts is several times slower. libdeflate picks its CPU-specific
# code paths at run time, so it keeps the compiler's default instruction set.
# Once done this will define
# LIBDEFLATE_INCLUDE_DIR - directory with libdeflate.h
# LIBDEFLATE_LIB_DIR - directory with libdeflate.a
# LIBDEFLATE_LIBRARY - the static library
# libdeflate_ep - target that builds the library when it is missing

include(FindPackageHandleStandardArgs)

set(LIBDEFLATE_INCLUDE_DIR "${LIBDEFLATE_ROOT_DIR}/include")
set(LIBDEFLATE_LIB_DIR "${LIBDEFLATE_ROOT_DIR}/lib")
set(LIBDEFLATE_LIBRARY "${LIBDEFLATE_LIB_DIR}/libdeflate.a")

add_custom_target(libdeflate_ep)
if(NOT EXISTS "${LIBDEFLATE_LIBRARY}")
  message(STATUS "Building libdeflate ${GENOMICSDB_LIBDEFLATE_VERSION} in ${LIBDEFLATE_ROOT_DIR}")
  # Workaround for issues with semicolon separated substrings to ExternalProject_Add
  # https://discourse.cmake.org/t/how-to-pass-cmake-osx-architectures-to-externalproject-add/2262
  if(CMAKE_OSX_ARCHITECTURES)
    string(REPLACE ";" "$<SEMICOLON>" _libdeflate_osx_architectures "${CMAKE_OSX_ARCHITECTURES}")
    set(_libdeflate_osx_args -DCMAKE_OSX_ARCHITECTURES=${_libdeflate_osx_architectures})
  endif()
  include(ExternalProject)
  ExternalProject_Add(libdeflate_build
    PREFIX ${CMAKE_BINARY_DIR}/dependencies/libdeflate
    URL https://github.com/ebiggers/libdeflate/releases/download/v${GENOMICSDB_LIBDEFLATE_VERSION}/libdeflate-${GENOMICSDB_LIBDEFLATE_VERSION}.tar.gz
    URL_HASH ${GENOMICSDB_LIBDEFLATE_URL_HASH}
    CMAKE_ARGS ${_libdeflate_osx_args}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_INSTALL_PREFIX=${LIBDEFLATE_ROOT_DIR}
      -DCMAKE_INSTALL_LIBDIR=lib
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON
      -DLIBDEFLATE_BUILD_SHARED_LIB=OFF
      -DLIBDEFLATE_BUILD_GZIP=OFF
    BUILD_BYPRODUCTS ${LIBDEFLATE_LIBRARY}
    )
  add_dependencies(libdeflate_ep libdeflate_build)
endif()

find_package_handle_standard_args(libdeflate DEFAULT_MSG LIBDEFLATE_INCLUDE_DIR LIBDEFLATE_LIBRARY)
