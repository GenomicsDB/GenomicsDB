# Following the standard CMake FindProtobuf module
# Determine compiler flags for protobuf
# Once done this will define
# PROTOBUF_LIBRARY_FOUND - protobuf found

find_path(PROTOBUF_INCLUDE_DIRS NAMES google/protobuf/service.h PATHS "${PROTOBUF_INCLUDE_DIR}" "${PROTOBUF_LIBRARY}/include")
if(PROTOBUF_STATIC_LINKING OR BUILD_DISTRIBUTABLE_LIBRARY)
    find_library(PROTOBUF_LIBRARIES NAMES libprotobuf.a protobuf PATHS "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
else()
    find_library(PROTOBUF_LIBRARIES NAMES protobuf PATHS "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Protobuf "Could not find Protobuf headers and/or libraries\n${DEFAULT_MSG}" PROTOBUF_INCLUDE_DIRS PROTOBUF_LIBRARIES)
find_program(PROTOBUF_PROTOC_EXECUTABLE NAMES protoc PATHS "${PROTOBUF_PROTOC_EXECUTABLE}" "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/bin")
if(PROTOBUF_FOUND)
    set(PROTOBUF_LIBRARY ${PROTOBUF_LIBRARIES})
    set(PROTOBUF_INCLUDE_DIR ${PROTOBUF_INCLUDE_DIRS})
endif()
include(FindProtobuf)

function(GET_PROTOBUF_GENERATED_FILES PROTOBUF_REGENERATE CXX_SRCS CXX_HDRS CXX_SRCS_DIR CXX_HDRS_DIR)
  set(LOCAL_CXX_SRCS) #local to this function stack
  set(LOCAL_CXX_HDRS)
  if(${PROTOBUF_REGENERATE})
    PROTOBUF_GENERATE_CPP(LOCAL_CXX_SRCS LOCAL_CXX_HDRS ${ARGN})
    set(${CXX_HDRS_DIR} ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)
    set(${CXX_SRCS_DIR} ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)
  else()
    foreach(PROTO_FILE ${ARGN})
      get_filename_component(CURR_FILENAME ${PROTO_FILE} NAME_WE)
      list(APPEND LOCAL_CXX_SRCS "${${CXX_SRCS_DIR}}/${CURR_FILENAME}.pb.cc")
      list(APPEND LOCAL_CXX_HDRS "${${CXX_HDRS_DIR}}/${CURR_FILENAME}.pb.h")
    endforeach()
  endif()
  set(${CXX_SRCS} ${LOCAL_CXX_SRCS} PARENT_SCOPE)
  set(${CXX_HDRS} ${LOCAL_CXX_HDRS} PARENT_SCOPE)
endfunction()
