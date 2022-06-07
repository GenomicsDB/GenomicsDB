#
# CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2021-2022 Omics Data Automation, Inc.
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

set(PROTOBUF_PROTO_FILES
  genomicsdb_coordinates.proto
  genomicsdb_callsets_mapping.proto
  genomicsdb_export_config.proto
  genomicsdb_import_config.proto
  genomicsdb_vid_mapping.proto
  )

if(PROTOBUF_REGENERATE)
  set(PROTO_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR})
else()
  set(PROTO_SRC_DIR ${CMAKE_SOURCE_DIR}/src/main/protobuf-generated/cpp/src)
endif()

foreach(PROTO_FILE ${PROTOBUF_PROTO_FILES})
  get_filename_component(CURR_FILENAME ${PROTO_FILE} NAME_WE)
  list(APPEND PROTO_GENERATED_CXX_SRCS ${PROTO_SRC_DIR}/${CURR_FILENAME}.pb.cc)
  list(APPEND PROTO_GENERATED_CXX_HDRS ${PROTO_SRC_DIR}/${CURR_FILENAME}.pb.h)
  if(BUILD_FOR_PYTHON)
    list(APPEND PROTO_GENERATED_PYTHON_SRCS ${PROTO_SRC_DIR}/python/${CURR_FILENAME}_pb2.py)
  endif()
endforeach()

if(PROTOBUF_REGENERATE)
  set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_command(OUTPUT ${PROTO_GENERATED_CXX_SRCS} ${PROTO_GENERATED_CXX_HDRS}
    COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} --cpp_out=${CMAKE_CURRENT_BINARY_DIR} ${PROTOBUF_PROTO_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    BYPRODUCTS ${PROTO_GENERATED_CXX_SRCS} ${PROTO_GENERATED_CXX_HDRS}
    COMMENT "Generating source/header files from .proto files"
    DEPENDS protobuf_ep ${PROTOBUF_PROTO_FILES})
  if(BUILD_FOR_PYTHON)
    add_custom_command(OUTPUT ${PROTO_GENERATED_PYTHON_SRCS}
      COMMAND echo "*** Generating python protobuf" &&
	      mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/python &&
	      ${PROTOBUF_PROTOC_EXECUTABLE} --python_out=${CMAKE_CURRENT_BINARY_DIR}/python ${PROTOBUF_PROTO_FILES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      BYPRODUCTS ${PROTO_GENERATED_PYTHON_SRCS}
      COMMENT "Generating GenomicsDB protobuf python bindings"
      DEPENDS protobuf_ep  ${PROTOBUF_PROTO_FILES})
  endif()
  add_custom_target(PROTOBUF_GENERATED_CXX_TARGET DEPENDS
    ${PROTO_GENERATED_CXX_SRCS} ${PROTO_GENERATED_CXX_HDRS} ${PROTO_GENERATED_PYTHON_SRCS})
  set(PROTOBUF_GENERATED_JAVA_SRCS
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBColumn.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBCallsetsMapProto.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBExportConfiguration.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBImportConfiguration.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBVidMapProto.java
    )
else()
  set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/src/main/protobuf-generated/cpp/include")
  set(PROTOBUF_GENERATED_JAVA_SRCS
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBCallsetsMapProto.java
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBExportConfiguration.java
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBImportConfiguration.java
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBVidMapProto.java
    )
endif()

set(PROTOBUF_GENERATED_CXX_SRCS ${PROTO_GENERATED_CXX_SRCS} CACHE INTERNAL "Protocol buffers generated C++ sources")
set(PROTOBUF_GENERATED_CXX_HDRS ${PROTO_GENERATED_CXX_HDRS} CACHE INTERNAL "Protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS ${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS} CACHE INTERNAL "Directory containing Protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_PYTHON_SRCS ${PROTO_GENERATED_PYTHON_SRCS} CACHE INTERNAL "Generated Python protobuf sources")
set(PROTOBUF_GENERATED_JAVA_SRCS ${PROTOBUF_GENERATED_JAVA_SRCS} CACHE INTERNAL "Protocol buffers generated Java sources")
