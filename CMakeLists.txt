#
# CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2019-2022 Omics Data Automation, Inc.
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

project(GenomicsDB)
cmake_minimum_required(VERSION 3.17)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  message(FATAL_ERROR CMAKE_CXX_COMPILER_ID " Compiler not yet supported. Exiting.")
endif()

#External dependencies. Can be overridden to point to custom locations
set(HTSLIB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/htslib" CACHE PATH "Path to htslib source directory")
set(TILEDB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/TileDB" CACHE PATH "Path to TileDB source directory")
set(RAPIDJSON_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/dependencies/RapidJSON/include" CACHE PATH "Path to RapidJSON include directory")
set(SPDLOG_INCLUDE_DIR  "${CMAKE_SOURCE_DIR}/dependencies/spdlog/include" CACHE PATH "Path to SpdLog include directory")
#The fmt headers are in dependencies/spdlog/include/spdlog/fmt/bundled - source files will include "fmt/bundled/format.h"
set(FMT_ROOT_DIR "${SPDLOG_INCLUDE_DIR}" CACHE PATH "Path to fmt library")
set(USING_FMT_BUNDLED_WITH_SPDLOG True CACHE BOOL "Disable if using fmt library that's not bundled with spdlog")

#GIT_COMMIT_HASH derived from git tag to be used as part of GENOMICSDB_VERSION
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

#Build parameters/options
set(GENOMICSDB_RELEASE_VERSION "1.4.4-SNAPSHOT" CACHE STRING "GenomicsDB release version") #used in Maven builds
set(GENOMICSDB_VERSION "${GENOMICSDB_RELEASE_VERSION}-${GIT_COMMIT_HASH}" CACHE STRING "GenomicsDB full version string")

set(DISABLE_MPI False CACHE BOOL "Disable use of any MPI compiler/libraries")
set(DISABLE_OPENMP False CACHE BOOL "Disable OpenMP")

set(BUILD_DISTRIBUTABLE_LIBRARY False CACHE BOOL "Build the TileDB/GenomicsDB library with minimal runtime dependencies")
set(BUILD_JAVA False CACHE BOOL "Build Java/JNI interface for combined VCF records")

set(HTSLIB_INSTALL_DIR "" CACHE PATH "Path to htslib install directory")

#Sync the GENOMICSDB_PROTOBUF_VERSION with protobuf.version in pom.xml!
set(GENOMICSDB_PROTOBUF_VERSION "3.19.4" CACHE STRING "Version of Google Protocol Buffer library")
set(PROTOBUF_ROOT_DIR "$ENV{HOME}/protobuf-install/${GENOMICSDB_PROTOBUF_VERSION}" CACHE PATH "Path to installed protobuf")
set(PROTOBUF_URL "https://github.com/protocolbuffers/protobuf/archive/v${GENOMICSDB_PROTOBUF_VERSION}.zip" CACHE STRING "URL to protobuf github release")
set(PROTOBUF_REGENERATE True CACHE BOOL "Regenerate protocol buffers C++ files")
set(GENERATE_PROTOBUF_FILES_IN_BUILD_DIR True CACHE BOOL "Generate all protobuf files in build directory")

set(USE_LIBCSV False CACHE BOOL "Disable library components that import data from csv files")
set(LIBCSV_DIR "" CACHE PATH "Path to libcsv header and library")

set(DO_PROFILING False CACHE BOOL "Collect some stats during execution - doesn't add much overhead")
set(COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW False CACHE BOOL "Counting #cells traversed in columnar iterator")
set(PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL False CACHE BOOL "Counting #cells traversed in columnar iterator in sliding window")
set(DO_MEMORY_PROFILING False CACHE BOOL "Collect memory consumption in parts of the combine gVCF program - high overhead")

#See https://gperftools.github.io/gperftools/cpuprofile.html
set(USE_GPERFTOOLS False CACHE BOOL "Collect profiling information using gperftools/profiler for analysis with pprof")
#See https://gperftools.github.io/gperftools/heapprofile.html
set(USE_GPERFTOOLS_HEAP False CACHE BOOL "Collect heap profiling information using gperftools/tcmalloc for analysis with pprof")

set(GENOMICSDB_MAVEN_BUILD_DIR ${CMAKE_BINARY_DIR}/target CACHE PATH "Path to maven build directory")
set(MAVEN_QUIET False CACHE BOOL "Do not print mvn messages")
set(GENOMICSDB_MAVEN_PROFILE "" CACHE STRING "Profile to compile with Spark 2 or Java 11") # used in Maven builds
set(GPG_PASSPHRASE "" CACHE STRING "Gpg passphrase for deploying to maven")

set(IPPROOT "" CACHE PATH "Path to IPP libraries - used when intel optimized zlib is required")

set(LIBDBI_DIR "" CACHE PATH "Path to libdbi install directory")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(USE_HDFS True CACHE BOOL "Enables HDFS support")
set(BUILD_FOR_ARCH "" CACHE STRING "Arch for which to build - values passed to -march= option for gcc")

#See https://cmake.org/Wiki/CMake_RPATH_handling#Common_questions
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)
if(APPLE)
  set(CMAKE_MACOSX_RPATH True CACHE BOOL "Set rpath on OSX")
  set(CMAKE_FIND_FRAMEWORK LAST CACHE STRING "Try to find frameworks on Mac OS X last")
endif()

enable_testing()

#Platform check modules
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)
   
#Optimization flag
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND BUILD_FOR_ARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${BUILD_FOR_ARCH}")
endif()

if(BUILD_DISTRIBUTABLE_LIBRARY)
    set(DISABLE_MPI True)
    set(DISABLE_OPENMP True)
    #For the GNU compiler, link in static gcc libraries
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    endif()
    if(APPLE)  
        set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13")
    endif()
endif()

#Check C++ 2011 support
include(CheckAndSetStdCXX2011Flag)
CHECK_AND_SET_STD_CXX_2011_FLAG(test_cpp_2011)
if(NOT test_cpp_2011)
    message(FATAL_ERROR "Your compiler does not support C++ 2011, exiting")
endif()
#Check for stack-protector-strong
include(CheckAndSetStackProtectorStrong)
CHECK_AND_SET_STACK_PROTECTOR_STRONG_FLAG(test_stack_protector_strong)
if(NOT test_stack_protector_strong)
    message(STATUS "Using -fstack-protector instead of -fstack-protector-strong")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstack-protector")
endif()

#Compiler Warning Options
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -z noexecstack -z relro -z now")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z noexecstack -z relro -z now")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unused-variable -Wno-pessimizing-move -Wno-deprecated-declarations")
endif()

#MPI library/compiler etc
if(NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
else()
    add_definitions(-DDISABLE_MPI=1)
endif()

if(NOT DISABLE_OPENMP)
    find_package(OpenMPv4)
    if(NOT OPENMPV4_FOUND)
        set(DISABLE_OPENMP True)
    else()
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()
if(DISABLE_OPENMP)
    add_definitions(-DDISABLE_OPENMP=1)
endif()

add_definitions(-D_FILE_OFFSET_BITS=64)  #large file support
add_definitions(-DDUPLICATE_CELL_AT_END=1) #mandatory
add_definitions(-DGENOMICSDB_VERSION=\"${GENOMICSDB_VERSION}\")

if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_definitions(-D_FORTIFY_SOURCE=2)
endif()

find_package(StringView)
if(NOT STRING_VIEW_FOUND AND NOT STRING_VIEW_EXPERIMENTAL_FOUND)
  message(FATAL_ERROR "No string_view found. Use a compiler that supports std::string_view for now")
endif()

find_package(libuuid REQUIRED)
include_directories(${LIBUUID_INCLUDE_DIR})

#OpenSSL required by TileDB and CURL
if(BUILD_DISTRIBUTABLE_LIBRARY)
  set(OPENSSL_USE_STATIC_LIBS TRUE)
endif()
find_package(OpenSSL REQUIRED)
find_package(CURL REQUIRED)

#RapidJSON
find_package(RapidJSON REQUIRED)
include_directories(${RAPIDJSON_INCLUDE_DIR})

#SpdLog
find_package(Spdlog REQUIRED)
include_directories(${SPDLOG_INCLUDE_DIR})

#htslib
find_package(htslib REQUIRED)
include_directories(BEFORE ${HTSLIB_INCLUDE_DIR})

#TileDB
find_package(TileDB REQUIRED)
include_directories(${TILEDB_INCLUDE_DIR})
#include_directories(${TILEDB_INCLUDE_DIR})

#librt
find_library(LIBRT_LIBRARY rt)

#fmt library
find_package(fmt REQUIRED)
include_directories(${FMT_INCLUDE_DIR})
if(USING_FMT_BUNDLED_WITH_SPDLOG)
  add_definitions(-DUSING_FMT_BUNDLED_WITH_SPDLOG)
endif()

#pgsql driver and dbi libs
#find_package(libdbi)
#if(LIBDBI_FOUND)
    #include_directories(${LIBDBI_INCLUDE_DIR})
    #add_definitions(-DLIBDBI)
#endif()

#Protobuf library
find_package(ProtobufWrapper REQUIRED)
if(GENERATE_PROTOBUF_FILES_IN_BUILD_DIR)
  set(PROTOBUF_JAVA_OUTPUT_DIR "${CMAKE_BINARY_DIR}/src/main/java")
else()
  set(PROTOBUF_JAVA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/src/main/java")
endif()
include_directories(${PROTOBUF_INCLUDE_DIRS})

#JNI
if(BUILD_JAVA)
    if(CMAKE_CROSSCOMPILING)
	set(JAVA_HOME ${TARGET_JAVA_HOME})
	set(ENV{JAVA_HOME} ${TARGET_JAVA_HOME})
    endif()
    find_package(JNI REQUIRED)
    if(CMAKE_CROSSCOMPILING)
	unset(JAVA_HOME)
	unset(ENV{JAVA_HOME})
    endif()
    find_package(Java 1.8 REQUIRED)
    include_directories(${JNI_INCLUDE_DIRS})
    include_directories(src/main/jni/include)
    include(UseJava)
endif()

#libcsv
if(NOT BUILD_DISTRIBUTABLE_LIBRARY)
    if(USE_LIBCSV)
        find_package(libcsv REQUIRED)
    else()
        find_package(libcsv)
    endif()
endif()
if(LIBCSV_FOUND)
    include_directories(${LIBCSV_INCLUDE_DIR})
    set(USE_LIBCSV True)
    add_definitions(-DUSE_LIBCSV)
endif()

#Collect stats
if(DO_PROFILING)
    message(STATUS "Enabling stats collection")
    add_definitions(-DDO_PROFILING=1)
    #Collect some more stats
    if(COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW)
        message(STATUS "Enabling COUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW")
        add_definitions(-DCOUNT_NUM_CELLS_BETWEEN_TWO_CELLS_FROM_THE_SAME_ROW=1)
    endif()
    if(PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL)
        message(STATUS "Enabling PROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL")
        add_definitions(-DPROFILE_NUM_CELLS_TO_TRAVERSE_AT_EVERY_QUERY_INTERVAL=1)
    endif()
endif()

#Collect stats to use with gperf/pprof
#See https://gperftools.github.io/gperftools/cpuprofile.html
#See https://gperftools.github.io/gperftools/heapprofile.html
if(USE_GPERFTOOLS OR USE_GPERFTOOLS_HEAP)
    find_package(GPerftools)
    if(GPERFTOOLS_FOUND OR GPERFTOOLS_HEAP_FOUND)
        message(STATUS "Enabling profiling using GPerftools")
        include_directories(${GPERFTOOLS_INCLUDE_DIR})
        if(USE_GPERFTOOLS)
            add_definitions(-DUSE_GPERFTOOLS)
        endif()
        if(USE_GPERFTOOLS_HEAP)
            add_definitions(-DUSE_GPERFTOOLS_HEAP)
        endif()
    else()
        message(WARNING "GPerftools headers/library not found")
    endif()
endif()

#Collect memory consumption stats while producing combined VCF records
if(DO_MEMORY_PROFILING)
    message(STATUS "Enabling memory consumption profiling while producing combined VCF records")
    add_definitions(-DDO_MEMORY_PROFILING=1)
endif()

#Produce verbose output
if(VERBOSE)
    add_definitions(-DVERBOSE=${VERBOSE})
    message(STATUS "Programs will produce output with verbosity=${VERBOSE}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g3 -gdwarf-3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_COVERAGE "-DDEBUG -g3 -gdwarf-3 --coverage")

#Clean all
add_custom_target(clean-all ${CMAKE_MAKE_PROGRAM} clean)
if(TILEDB_SOURCE_DIR)
    add_custom_target(clean-TileDB ${CMAKE_COMMAND} -E remove_directory ${TILEDB_DIR_IN_BUILD_DIR})
    add_dependencies(clean-all clean-TileDB)
endif()
if(HTSLIB_SOURCE_DIR)
    add_custom_target(clean-htslib ${CMAKE_COMMAND} -E remove_directory ${HTSLIB_DIR_IN_BUILD_DIR})
    add_dependencies(clean-all clean-htslib)
endif()

#Uninstall GenomicsDB
add_custom_target(
   uninstall-genomicsdb
   COMMAND echo "-- Uninstalling GenomicsDB from ${CMAKE_INSTALL_PREFIX}..."
   COMMAND < install_manifest.txt xargs -I % sh -c 'echo -- Uninstalling % && rm -f %' %
   COMMAND echo "-- GenomicsDB is uninstalled!"
)

find_package(safestringlib)
if(SAFESTRINGLIB_FOUND)
    add_definitions(-DSAFESTRINGLIB_FOUND=1)
    include_directories(${SAFESTRINGLIB_INCLUDE_DIR})
endif()

include_directories (
    src/main/cpp/include/genomicsdb
    src/main/cpp/include/loader
    src/main/cpp/include/query_operations
    src/main/cpp/include/utils
    src/main/cpp/include/vcf
    src/main/cpp/include/config
    src/main/cpp/include/api
    src/test/cpp/include
)

set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${PROTOBUF_LIBRARIES})
if(MPI_FOUND)
  set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()
if(USE_HDFS)
    set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES} ${JAVA_JVM_LIBRARY})
endif()
if(OPENSSL_FOUND)
    set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES} ${OPENSSL_LIBRARIES})
endif()
if(CURL_FOUND)
    set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES} ${CURL_LIBRARIES})
endif()
if(SAFESTRINGLIB_FOUND)
    set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES} ${SAFESTRINGLIB_LIBRARY})
endif()
set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES}
  ${OPENSSL_LIBRARIES} ${ZLIB_LIBRARIES} ${LIBUUID_LIBRARY} ${CMAKE_DL_LIBS})

if(IPPROOT AND ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  find_package(IPP REQUIRED)
  set(GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES}
    ${LIBIPPDC_LIBRARY} ${LIBIPPS_LIBRARY} ${LIBIPPCORE_LIBRARY})
endif()

function(build_GenomicsDB_executable_common target)
  build_GenomicsDB_links(${target})
  install(TARGETS ${target} RUNTIME DESTINATION bin)
endfunction()

function(build_GenomicsDB_links target)
    set_property(TARGET ${target} PROPERTY POSITION_INDEPENDENT_CODE ON)
    #add_executable(${target} src/${target}.cc)
    target_link_libraries(${target} genomicsdb)
    if(TILEDB_SOURCE_DIR)
        target_link_libraries(${target} tiledb_static)
    else()
        target_link_libraries(${target} ${TILEDB_LIBRARY})
    endif()
    if(HTSLIB_SOURCE_DIR)
        add_dependencies(${target} htslib)
    endif()
    target_link_libraries(${target} ${HTSLIB_LIBRARY})
    if(LIBCSV_FOUND)
        target_link_libraries(${target} ${LIBCSV_LIBRARY})
    endif()
    if(LIBDBI_FOUND)
        target_link_libraries(${target} ${LIBPGSQL_DRIVER_LIBRARY} ${LIBDBI_DEV_LIBRARY})
    endif()
    if(LIBRT_LIBRARY)
        target_link_libraries(${target} ${LIBRT_LIBRARY})
    endif()
    if(USE_GPERFTOOLS AND GPERFTOOLS_FOUND)
        target_link_libraries(${target} ${GPERFTOOLS_PROFILER_LIBRARY})
    endif()
    if(USE_GPERFTOOLS_HEAP AND GPERFTOOLS_FOUND)
        target_link_libraries(${target} ${GPERFTOOLS_HEAP_PROFILER_LIBRARY})
    endif()
    target_link_libraries(${target} ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES})
endfunction()

function(build_GenomicsDB_executable target)
    add_executable(${target} src/${target}.cc)
    build_GenomicsDB_executable_common(${target})
endfunction()

add_subdirectory(src)
include_directories(${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS})
add_subdirectory(tools)
add_subdirectory(example)

if(BUILD_JAVA)
    set(JAVA_SCALA_SOURCES
        ./src/main/java/org/genomicsdb/exception/GenomicsDBException.java
        ./src/main/java/org/genomicsdb/importer/extensions/CallSetMapExtensions.java
        ./src/main/java/org/genomicsdb/importer/extensions/JsonFileExtensions.java
        ./src/main/java/org/genomicsdb/importer/extensions/VidMapExtensions.java
        ./src/main/java/org/genomicsdb/importer/model/ChromosomeInterval.java
        ./src/main/java/org/genomicsdb/importer/model/SampleInfo.java
        ./src/main/java/org/genomicsdb/importer/Constants.java
        ./src/main/java/org/genomicsdb/importer/GenomicsDBImporter.java
        ./src/main/java/org/genomicsdb/importer/GenomicsDBImporterJni.java
        ./src/main/java/org/genomicsdb/importer/GenomicsDBImporterStreamWrapper.java
        ./src/main/java/org/genomicsdb/importer/MultiChromosomeIterator.java
        ./src/main/java/org/genomicsdb/importer/SilentByteBufferStream.java
        ./src/main/java/org/genomicsdb/model/CommandLineImportConfig.java
        ./src/main/java/org/genomicsdb/model/ImportConfig.java
        ./src/main/java/org/genomicsdb/model/BatchCompletionCallbackFunctionArgument.java
        ./src/main/java/org/genomicsdb/reader/GenomicsDBFeatureIterator.java
        ./src/main/java/org/genomicsdb/reader/GenomicsDBFeatureReader.java
        ./src/main/java/org/genomicsdb/reader/GenomicsDBQueryStream.java
        ./src/main/java/org/genomicsdb/reader/GenomicsDBTimer.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBConfiguration.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBInputFormat.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBInputSplit.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBJavaSparkFactory.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBPartitionInfo.java
	./src/main/java/org/genomicsdb/spark/GenomicsDBQueryInfo.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBSchemaFactory.java
        ./src/main/java/org/genomicsdb/spark/GenomicsDBInput.java
        ./src/main/java/org/genomicsdb/Constants.java
        ./src/main/java/org/genomicsdb/GenomicsDBLibLoader.java
        ./src/main/scala/org/genomicsdb/GenomicsDBContext.scala
        ./src/main/scala/org/genomicsdb/GenomicsDBPartition.scala
        ./src/main/scala/org/genomicsdb/GenomicsDBRDD.scala
        ./src/main/scala/org/genomicsdb/GenomicsDBScalaSparkFactory.scala
        ./src/test/java/org/genomicsdb/importer/GenomicsDBImporterSpec.java
        ./src/test/java/org/genomicsdb/reader/ChrArrayFolderComparatorSpec.java
        ./src/test/java/org/genomicsdb/spark/GenomicsDBInputFormatTest.java
        ./src/test/java/org/genomicsdb/model/CommandLineImportConfigSpec.java
        ./src/test/java/org/genomicsdb/model/GenomicsDBCallsetsMapProtoSpec.java
        ./src/test/java/org/genomicsdb/model/GenomicsDBExportConfigurationSpec.java
        ./src/test/java/org/genomicsdb/model/GenomicsDBImportConfigurationSpec.java
        ./src/test/java/org/genomicsdb/model/GenomicsDBVidMappingProtoSpec.java
        ./src/test/java/org/genomicsdb/model/ImportConfigSpec.java
        ./src/test/java/org/genomicsdb/GenomicsDBTestUtils.java
        )
    
      if(GENOMICSDB_MAVEN_PROFILE STREQUAL "spark2")
      set(JAVA_SCALA_SOURCES 
          ${JAVA_SCALA_SOURCES}
          ./src/main/java/org/genomicsdb/spark/v2/GenomicsDBDataSourceReader.java
          ./src/main/java/org/genomicsdb/spark/v2/GenomicsDBDataSourceV2.java
          ./src/main/java/org/genomicsdb/spark/v2/GenomicsDBRecordReader.java
          ./src/main/java/org/genomicsdb/spark/v2/GenomicsDBInputPartition.java
	  ./src/main/java/org/genomicsdb/spark/v2/GenomicsDBInputPartitionReader.java
        )
    else()
      set(JAVA_SCALA_SOURCES
          ${JAVA_SCALA_SOURCES}
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBRecordReader.java
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBSource.java
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBTable.java
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBScan.java
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBScanBuilder.java
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBBatch.java
          ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBInputPartition.java
	  ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBInputPartitionReader.java
	  ./src/main/java/org/genomicsdb/spark/sources/GenomicsDBPartitionReaderFactory.java
        )
    endif()
    
    if(PROTOBUF_REGENERATE)
        #Must be set here - see https://cmake.org/cmake/help/v3.3/command/set_source_files_properties.html
        set_source_files_properties(${PROTOBUF_GENERATED_JAVA_SRCS} PROPERTIES GENERATED True)
        #Protobuf files are generated by Maven
        set(MAVEN_PROTOBUF_REGENERATE_ARGS "-Dprotoc.filepath=${PROTOBUF_PROTOC_EXECUTABLE}"
      "-Dprotobuf.output.directory=${PROTOBUF_JAVA_OUTPUT_DIR}")
    else()
        set(JAVA_SCALA_SOURCES ${JAVA_SCALA_SOURCES} ${PROTOBUF_GENERATED_JAVA_SRCS})
        set(MAVEN_PROTOBUF_REGENERATE_ARGS "")
    endif()
    set(MAVEN_QUIET_ARGS "")
    if(MAVEN_QUIET)
        set(MAVEN_QUIET_ARGS "--quiet")
    endif()
    set(MAVEN_ARGS
        ${MAVEN_QUIET_ARGS}
        -Dgenomicsdb.version=${GENOMICSDB_RELEASE_VERSION}
        -Dgenomicsdb_source_directory=${CMAKE_SOURCE_DIR}
        -Dgenomicsdb_build_directory=${GENOMICSDB_MAVEN_BUILD_DIR}
        -Dgenomicsdb_lib_directory=$<TARGET_FILE_DIR:tiledbgenomicsdb>
        -Dprotobuf.version=${GENOMICSDB_PROTOBUF_VERSION}
        ${MAVEN_PROTOBUF_REGENERATE_ARGS})
    
      if (GPG_PASSPHRASE) 
        set(MAVEN_ARGS 
          ${MAVEN_ARGS}
          -Dgpg.passphrase=${GPG_PASSPHRASE})
        message(STATUS "Setting gpg passphrase for maven")

      endif()

      if(GENOMICSDB_MAVEN_PROFILE)
        message(STATUS "Setting maven profile to ${GENOMICSDB_MAVEN_PROFILE}" )
        set(MAVEN_PROFILE 
          -P${GENOMICSDB_MAVEN_PROFILE})
        set(MAVEN_ARGS ${MAVEN_ARGS} ${MAVEN_PROFILE})
    else()
        set(MAVEN_PROFILE "")
    endif()

    #Maven build - depends on dynamic library
    add_custom_command(
        OUTPUT ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}.jar ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}-allinone.jar
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/pom.xml ${CMAKE_BINARY_DIR}/pom.xml
        COMMAND mvn versions:set ${MAVEN_QUIET_ARGS} -DnewVersion=${GENOMICSDB_RELEASE_VERSION} ${MAVEN_PROFILE}
        COMMAND mvn package -DskipTests ${MAVEN_ARGS}
        DEPENDS tiledbgenomicsdb ${JAVA_SCALA_SOURCES} pom.xml
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    install(FILES ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}.jar ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}-allinone.jar DESTINATION bin)

    execute_process(
      COMMAND ln -sf ${CMAKE_SOURCE_DIR}/tests ${CMAKE_BINARY_DIR})

    add_test(NAME javatests
        COMMAND mvn test
        ${MAVEN_ARGS}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    set(GENOMICSDB_EXAMPLE_SOURCES
        ${CMAKE_SOURCE_DIR}/example/java/TestGenomicsDB.java
        ${CMAKE_SOURCE_DIR}/example/java/TestBufferStreamGenomicsDBImporter.java
        ${CMAKE_SOURCE_DIR}/example/java/TestGenomicsDBSparkHDFS.java
        ${CMAKE_SOURCE_DIR}/example/java/TestGenomicsDBSource.java
        ${CMAKE_SOURCE_DIR}/example/java/TestGenomicsDBImporterWithMergedVCFHeader.java)

    add_jar(genomicsdb-${GENOMICSDB_RELEASE_VERSION}-examples
        SOURCES ${GENOMICSDB_EXAMPLE_SOURCES}
                log4j.properties
        INCLUDE_JARS ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}-allinone.jar
        OUTPUT_DIR ${GENOMICSDB_MAVEN_BUILD_DIR})

    #Deploy to Maven central
    #Cannot specify the jars in the DEPENDS clause. CMake manual specifies, quoting:
    #"Do not list the output in more than one independent target that may build in parallel or
    #the two instances of the rule may conflict"
    #Jars are already listed as dependencies in the add_jar command (-examples)
    #So, use the examples target as a dependency for this target
    add_custom_target(mvn_central_deploy
      COMMAND mvn ${MAVEN_ARGS} deploy
        COMMENT "Deploy to Maven central"
    )
    add_dependencies(mvn_central_deploy genomicsdb-${GENOMICSDB_RELEASE_VERSION}-examples)

    add_custom_target(mvn_local_deploy
            COMMAND mvn install:install-file -Dfile=target/genomicsdb-${GENOMICSDB_RELEASE_VERSION}.jar -DpomFile=pom.xml)
    add_dependencies(mvn_local_deploy genomicsdb-${GENOMICSDB_RELEASE_VERSION}-examples)

    if(NOT DISABLE_MPI)
        add_test(NAME CI_tests
            COMMAND python ${CMAKE_SOURCE_DIR}/tests/run.py
            ${CMAKE_BINARY_DIR}
            ${CMAKE_INSTALL_PREFIX}
            ${CMAKE_BUILD_TYPE})
    endif()
endif()
