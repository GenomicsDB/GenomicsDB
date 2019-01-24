#
# FindTileDB.cmake
#
# The MIT License
#
# Copyright (c) 2016 MIT and Intel Corporation
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
# Determine compiler flags for TileDB
# Once done this will define
# TILEDB_FOUND - TileDB found

set(TILEDB_VERBOSE True)
set(TILEDB_DISABLE_TESTING True)

#Zlib
find_package(ZLIB REQUIRED)

#OpenSSL
if(OPENSSL_PREFIX_DIR AND NOT OPENSSL_ROOT_DIR)
    set(OPENSSL_ROOT_DIR "${OPENSSL_PREFIX_DIR}")
endif()
if(BUILD_DISTRIBUTABLE_LIBRARY)
    set(OPENSSL_USE_STATIC_LIBS True)
endif()
find_package(OpenSSL REQUIRED) #now performed inside TileDB

#libuuid
find_package(libuuid REQUIRED)

include(FindPackageHandleStandardArgs)


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
    add_subdirectory(${TILEDB_SOURCE_DIR} EXCLUDE_FROM_ALL)
else()
    find_path(TILEDB_INCLUDE_DIR NAMES "tiledb.h" HINTS "${TILEDB_INSTALL_DIR}/include")
    find_library(TILEDB_LIBRARY NAMES libtiledb.a tiledb HINTS "${TILEDB_INSTALL_DIR}" "${TILEDB_INSTALL_DIR}/lib64" "${TILEDB_INSTALL_DIR}/lib")
    find_package_handle_standard_args(TileDB "Could not find TileDB headers and/or libraries ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR TILEDB_LIBRARY)
endif()
