set(PROTOBUF_PROTO_FILES
    genomicsdb_coordinates.proto
    genomicsdb_callsets_mapping.proto
    genomicsdb_export_config.proto
    genomicsdb_import_config.proto
    genomicsdb_vid_mapping.proto
    )

#the following 2 variable values are useful iff PROTOBUF_REGENERATE is false 
set(CXX_SRCS_DIR "${CMAKE_SOURCE_DIR}/protobuf-generated/cpp/src")
set(CXX_HDRS_DIR "${CMAKE_SOURCE_DIR}/protobuf-generated/cpp/include")

GET_PROTOBUF_GENERATED_FILES(${PROTOBUF_REGENERATE} CXX_SRCS CXX_HDRS
  CXX_SRCS_DIR CXX_HDRS_DIR
  ${PROTOBUF_PROTO_FILES})

add_custom_target(PROTOBUF_GENERATED_CXX_TARGET DEPENDS ${CXX_SRCS} ${CXX_HDRS})

add_subdirectory(ga4gh)
list(APPEND CXX_SRCS ${PROTOBUF_GENERATED_GA4GH_CXX_SRCS})
list(APPEND CXX_HDRS ${PROTOBUF_GENERATED_GA4GH_CXX_HDRS})
list(APPEND CXX_HDRS_DIR ${PROTOBUF_GENERATED_GA4GH_CXX_HDRS_DIR})
add_dependencies(PROTOBUF_GENERATED_CXX_TARGET PROTOBUF_GENERATED_GA4GH_CXX_TARGET)

if(PROTOBUF_REGENERATE)
  set(PROTOBUF_GENERATED_JAVA_SRCS
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBColumn.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBCallsetsMapProto.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBExportConfiguration.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBImportConfiguration.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/GenomicsDBVidMapProto.java
    ${PROTOBUF_GENERATED_GA4GH_JAVA_SRCS}
    )
else()
  set(PROTOBUF_GENERATED_JAVA_SRCS
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBCallsetsMapProto.java
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBExportConfiguration.java
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBImportConfiguration.java
    src/main/protobuf-generated/java/org/genomicsdb/GenomicsDBVidMapProto.java
    ${PROTOBUF_GENERATED_GA4GH_JAVA_SRCS}
    )
endif()

set(PROTOBUF_GENERATED_CXX_SRCS ${CXX_SRCS} CACHE INTERNAL "Protocol buffers generated C++ sources")
set(PROTOBUF_GENERATED_CXX_HDRS ${CXX_HDRS} CACHE INTERNAL "Protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS ${CXX_HDRS_DIR} CACHE INTERNAL "Directory containing Protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_JAVA_SRCS ${PROTOBUF_GENERATED_JAVA_SRCS} CACHE INTERNAL "Protocol buffers generated Java sources")
