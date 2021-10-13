#
# FindTileDB.cmake
#
# The MIT License
#
# Copyright (c) 2016 MIT and Intel Corporation
# Copyright (c) 2019-2020 Omics Data Automation, Inc.
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
# Determine compiler flags for TileDB
# Once done this will define
# TILEDB_FOUND - TileDB found

# Update git submodule if necessary
if((DEFINED TILEDB_SOURCE_DIR) AND (NOT "${TILEDB_SOURCE_DIR}" STREQUAL "") AND (NOT EXISTS "${TILEDB_SOURCE_DIR}/CMakeLists.txt"))
  MESSAGE(STATUS "Installing submodule TileDB at ${CMAKE_SOURCE_DIR}/dependencies")
  execute_process(
    COMMAND git submodule update --recursive --init ${CMAKE_SOURCE_DIR}/dependencies/TileDB
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE submodule_update_exit_code)
  if(NOT(submodule_update_exit_code EQUAL 0))
    message(FATAL_ERROR "Failure to update git submodule TileDB")
  endif()
  set(TILEDB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/TileDB" CACHE PATH "Path to TileDB source directory" FORCE)
endif()

#TileDB c-api returns a tiledb_errmsg which is what GenomicsDB should rely on
set(TILEDB_VERBOSE False)
set(TILEDB_DISABLE_TESTING True)

#Zlib
find_package(ZLIB REQUIRED)

#ZStd
set(ENABLE_ZSTD TRUE CACHE BOOL "Enable ZStd" FORCE)

if(BUILD_DISTRIBUTABLE_LIBRARY)
    set(OPENSSL_USE_STATIC_LIBS True)
endif()
find_package(OpenSSL REQUIRED) #now performed inside TileDB

#libuuid
find_package(libuuid REQUIRED)

include(FindPackageHandleStandardArgs)

set(TILEDB_DISABLE_TESTS True CACHE BOOL "Disable TileDB Testing")

#Build if TileDB source directory specified
if(TILEDB_SOURCE_DIR)
    #OpenMP
    if(DISABLE_OPENMP)
        set(USE_OPENMP False CACHE BOOL "Disable OpenMP" FORCE)
    else()
        set(USE_OPENMP True CACHE BOOL "Enable OpenMP" FORCE)
    endif()

    find_path(TILEDB_INCLUDE_DIR NAMES tiledb.h HINTS "${TILEDB_SOURCE_DIR}/core/include/c_api")
    find_package_handle_standard_args(TileDB "Could not find TileDB headers ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR)
    add_subdirectory(${TILEDB_SOURCE_DIR} dependencies/TileDB EXCLUDE_FROM_ALL)
else()
    find_path(TILEDB_INCLUDE_DIR NAMES "tiledb.h" HINTS "${TILEDB_INSTALL_DIR}/include")
    find_library(TILEDB_LIBRARY NAMES libtiledb.a tiledb HINTS "${TILEDB_INSTALL_DIR}" "${TILEDB_INSTALL_DIR}/lib64" "${TILEDB_INSTALL_DIR}/lib")
    find_package_handle_standard_args(TileDB "Could not find TileDB headers and/or libraries ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR TILEDB_LIBRARY)
endif()
