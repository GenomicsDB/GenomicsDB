#
# CMakeLists.txt
#
# The MIT License
#
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

if(PROTOBUF_STATIC_LINKING OR BUILD_DISTRIBUTABLE_LIBRARY)
    if(CMAKE_VERSION VERSION_GREATER 3.10.0)
        set(Protobuf_USE_STATIC_LIBS ON)
    endif()
    set(PROTOBUF_WRAPPER_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
    set(PROTOBUF_WRAPPER_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX}) 
endif()

if (PROTOBUF_LIBRARY)
    find_path(Protobuf_INCLUDE_DIR google/protobuf/service.h REQUIRED NO_DEFAULT_PATH PATHS "${PROTOBUF_LIBRARY}/include")
    find_library(Protobuf_LIBRARY libprotobuf${PROTOBUF_WRAPPER_LIBRARY_SUFFIX} REQUIRED NO_DEFAULT_PATH PATHS "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
    find_library(Protobuf_PROTOC_LIBRARY libprotoc${PROTOBUF_WRAPPER_LIBRARY_SUFFIX} REQUIRED NO_DEFAULT_PATH PATHS "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
    find_program(Protobuf_PROTOC_EXECUTABLE protoc REQUIRED NO_DEFAULT_PATH PATHS "${PROTOBUF_LIBRARY}/bin")
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(ProtobufWrapper "Check PROTOBUF_LIBRARY option, could not find Protobuf headers and/or binaries at ${PROTOBUF_LIBRARY}${DEFAULT_MSG}" Protobuf_INCLUDE_DIR Protobuf_LIBRARY Protobuf_PROTOC_LIBRARY Protobuf_PROTOC_EXECUTABLE)
endif()

find_package(Protobuf REQUIRED)

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
