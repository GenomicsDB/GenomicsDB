# Following the standard CMake FindProtobuf module
# Determine compiler flags for protobuf
# Once done this will define
# PROTOBUF_LIBRARY_FOUND - protobuf found

find_path(PROTOBUF_INCLUDE_DIRS NAMES google/protobuf/service.h PATHS "${PROTOBUF_INCLUDE_DIR}" "${PROTOBUF_LIBRARY}/include")
if(PROTOBUF_STATIC_LINKING OR BUILD_DISTRIBUTABLE_LIBRARY)
    if (CMAKE_VERSION VERSION_GREATER 3.10.0)
        set(Protobuf_USE_STATIC_LIBS ON)
    endif()
    find_library(PROTOBUF_LIBRARIES NAMES libprotobuf.a protobuf PATHS "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
else()
    find_library(PROTOBUF_LIBRARIES NAMES protobuf PATHS "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ProtobufWrapper "Could not find Protobuf headers and/or libraries\n${DEFAULT_MSG}" PROTOBUF_INCLUDE_DIRS PROTOBUF_LIBRARIES)
find_program(PROTOBUF_PROTOC_EXECUTABLE NAMES protoc PATHS "${PROTOBUF_PROTOC_EXECUTABLE}" "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/bin")
if(PROTOBUF_FOUND)
    set(PROTOBUF_LIBRARY ${PROTOBUF_LIBRARIES})
    set(PROTOBUF_INCLUDE_DIR ${PROTOBUF_INCLUDE_DIRS})
endif()
find_package(Protobuf)

Message(STATUS "Protobuf_LIBRARY=" ${Protobuf_LIBRARY})

include(CheckCXXSourceCompiles)
function(CHECK_IF_USING_PROTOBUF_V_3_0_0_BETA_1 FLAG_VAR_NAME)
  set(PB_test_source
    "
    #include <google/protobuf/util/json_util.h>
    int main() {
      google::protobuf::util::JsonParseOptions parse_opt;
      parse_opt.ignore_unknown_fields = true;
      return 0;
    }
    "
    )
  set(CMAKE_REQUIRED_INCLUDES ${PROTOBUF_INCLUDE_DIRS})
  check_cxx_source_compiles("${PB_test_source}" PROTOBUF_V3_STABLE_FOUND)
  if(PROTOBUF_V3_STABLE_FOUND)
    set(${FLAG_VAR_NAME} False PARENT_SCOPE)
  else()
    set(${FLAG_VAR_NAME} True PARENT_SCOPE)
  endif()
endfunction()
