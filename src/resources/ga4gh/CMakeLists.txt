set(PROTOBUF_PROTO_FILES
  common.proto
  variants.proto
    )

#the following 2 variable values are useful iff PROTOBUF_REGENERATE is false 
set(CXX_SRCS_DIR "${CMAKE_SOURCE_DIR}/protobuf-generated/cpp/src")
set(CXX_HDRS_DIR "${CMAKE_SOURCE_DIR}/protobuf-generated/cpp/include")

GET_PROTOBUF_GENERATED_FILES(${PROTOBUF_REGENERATE} CXX_SRCS CXX_HDRS
  CXX_SRCS_DIR CXX_HDRS_DIR
  ${PROTOBUF_PROTO_FILES})

if(PROTOBUF_REGENERATE)
  set(PROTOBUF_GENERATED_GA4GH_JAVA_SRCS
    #${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/ga4gh/Common.java
    ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/org/genomicsdb/ga4gh/Variants.java
    PARENT_SCOPE
    )
else()
  set(PROTOBUF_GENERATED_GA4GH_JAVA_SRCS
    #src/main/protobuf-generated/java/org/genomicsdb/ga4gh/Common.java
    src/main/protobuf-generated/java/org/genomicsdb/ga4gh/Variants.java
    PARENT_SCOPE
    )
endif()

add_custom_target(PROTOBUF_GENERATED_GA4GH_CXX_TARGET DEPENDS ${CXX_SRCS} ${CXX_HDRS})

set(PROTOBUF_GENERATED_GA4GH_CXX_SRCS ${CXX_SRCS} PARENT_SCOPE)
set(PROTOBUF_GENERATED_GA4GH_CXX_HDRS ${CXX_HDRS} PARENT_SCOPE)
set(PROTOBUF_GENERATED_GA4GH_CXX_HDRS_DIR ${CXX_HDRS_DIR} PARENT_SCOPE)
