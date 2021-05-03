#
# FindTileDB.cmake
#
# The MIT License
#
# Copyright (c) 2020-2021 Omics Data Automation, Inc.
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
# Determine compiler flags for htslib
# Once done this will define
# HTSLIB_FOUND - htslib found

include(FindPackageHandleStandardArgs)

# Update git submodule if necessary
if((DEFINED HTSLIB_SOURCE_DIR) AND (NOT "${HTSLIB_SOURCE_DIR}" STREQUAL "") AND (NOT EXISTS "${HTSLIB_SOURCE_DIR}/htslib/vcf.h"))
  MESSAGE(STATUS "Installing submodule htslib at ${CMAKE_SOURCE_DIR}/dependencies")
  execute_process(
    COMMAND git submodule update --init ${CMAKE_SOURCE_DIR}/dependencies/htslib
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE submodule_update_exit_code)
  if(NOT(submodule_update_exit_code EQUAL 0))
    message(FATAL_ERROR "Failure to update git submodule htslib")
  endif()
  set(HTSLIB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/htslib" CACHE PATH "Path to htslib source directory" FORCE)
endif()

#Build if htslib source directory specified
if(HTSLIB_SOURCE_DIR)
    set(HTSLIB_Debug_CFLAGS " -Wall -fPIC -DDEBUG  -g3 -gdwarf-3 -DVCF_ALLOW_INT64=1")  #will be picked if compiling in debug mode
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(HTSLIB_Debug_CFLAGS "${HTSLIB_Debug_CFLAGS} -Wno-expansion-to-defined -Wno-nullability-completeness")
    endif()
    set(HTSLIB_Coverage_CFLAGS "${HTSLIB_Debug_CFLAGS}")
    set(HTSLIB_Release_CFLAGS " -Wall -fPIC -O3 -DVCF_ALLOW_INT64=1")
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      set(HTSLIB_Release_CFLAGS "${HTSLIB_Release_CFLAGS} -Wno-expansion-to-defined -Wno-nullability-completeness")
    endif()
    if(APPLE)
      set(HTSLIB_LDFLAGS " -L /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem")
    endif()
    set(HTSLIB_Debug_LDFLAGS "-g3 -gdwarf-3 ${HTSLIB_LDFLAGS}")
    set(HTSLIB_Coverage_LDFLAGS "${HTSLIB_Debug_LDFLAGS}")
    set(HTSLIB_Release_LDFLAGS "${HTSLIB_LDFLAGS}")
    include(ExternalProject)
    #Cross compiling for MacOSX
    if((NOT (CMAKE_SYSTEM_NAME STREQUAL CMAKE_HOST_SYSTEM_NAME)) AND APPLE)
        set(HTSLIB_OSXCROSS_COMPILE_FLAGS LIBS=${OSXCROSS_LIBS} CPPFLAGS=${OSXCROSS_CPPFLAGS} --host=${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM})
    endif()
    set(HTSLIB_${CMAKE_BUILD_TYPE}_CFLAGS "${HTSLIB_${CMAKE_BUILD_TYPE}_CFLAGS} -I${OPENSSL_INCLUDE_DIR} -I${CURL_INCLUDE_DIRS}")
    if(APPLE AND BUILD_DISTRIBUTABLE_LIBRARY)
      set(HTSLIB_ENV "MACOSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}")
    endif()
    ExternalProject_Add(
        htslib
        DOWNLOAD_COMMAND ""
        SOURCE_DIR "${HTSLIB_SOURCE_DIR}"
        UPDATE_COMMAND autoreconf -i ${HTSLIB_SOURCE_DIR}
        PATCH_COMMAND ""
        CONFIGURE_COMMAND ${HTSLIB_ENV} ${HTSLIB_SOURCE_DIR}/configure
            CFLAGS=${HTSLIB_${CMAKE_BUILD_TYPE}_CFLAGS}
            LDFLAGS=${HTSLIB_${CMAKE_BUILD_TYPE}_LDFLAGS}
            CC=${CMAKE_C_COMPILER} AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB}
            ${HTSLIB_OSXCROSS_COMPILE_FLAGS}
            --disable-lzma --disable-bz2 --disable-s3 --disable-gcs
        BUILD_COMMAND ${CMAKE_COMMAND} -E make_directory cram
            COMMAND ${CMAKE_COMMAND} -E make_directory test
            COMMAND ${CMAKE_COMMAND} -E make_directory test/fuzz
            COMMAND ${CMAKE_COMMAND} -E copy ${HTSLIB_SOURCE_DIR}/version.sh .
            COMMAND ${HTSLIB_ENV} $(MAKE) -f ${HTSLIB_SOURCE_DIR}/Makefile VPATH=${HTSLIB_SOURCE_DIR} SOURCE_DIR=${HTSLIB_SOURCE_DIR} AR=${CMAKE_AR}
        #BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
        )
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/vcf.h HINTS "${HTSLIB_SOURCE_DIR}" CMAKE_FIND_ROOT_PATH_BOTH NO_DEFAULT_PATH)
    ExternalProject_Get_Property(htslib BINARY_DIR)
    set(HTSLIB_DIR_IN_BUILD_DIR "${BINARY_DIR}")
    set(HTSLIB_LIBRARY "${BINARY_DIR}/libhts.a")
    find_package_handle_standard_args(htslib "Could not find htslib headers ${DEFAULT_MSG}" HTSLIB_INCLUDE_DIR)
else()
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/vcf.h HINTS "${HTSLIB_INSTALL_DIR}")
    find_library(HTSLIB_LIBRARY NAMES libhts.a hts HINTS "${HTSLIB_INSTALL_DIR}")
    find_package_handle_standard_args(htslib "Could not find htslib headers and/or libraries ${DEFAULT_MSG}" HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)
endif()
