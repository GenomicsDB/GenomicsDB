#
# FindRapidJSON.cmake
#
# The MIT License
#
# Copyright (c) 2020 Omics Data Automation, Inc.
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
# Determine compiler flags for RapidJSON
# Once done this will define
# RAPIDJSON_FOUND - RapidJSON found
# RapidJSON_C_FLAGS
# RapidJSON_CXX_FLAGS

# Update git submodule if necessary
if(NOT EXISTS "${RAPIDJSON_INCLUDE_DIR}")
  MESSAGE(STATUS "Installing submodule RapidJSON at ${CMAKE_SOURCE_DIR}/dependencies")
  execute_process(
    COMMAND git submodule update --init ${CMAKE_SOURCE_DIR}/dependencies/RapidJSON
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE submodule_update_exit_code)
  if(NOT(submodule_update_exit_code EQUAL 0))
    message(FATAL_ERROR "Failure to update git submodule RapidJSON")
  endif()
  set(RAPIDJSON_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/dependencies/RapidJSON/include" CACHE PATH "Path to RapidJSON include directory" FORCE)
endif()

find_path(RAPIDJSON_INCLUDE_DIR
  NAMES rapidjson/rapidjson.h
  PATHS "${RAPIDJSON_INCLUDE_DIR}"
  CMAKE_FIND_ROOT_PATH_BOTH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RapidJSON "Could not find RapidJSON header: rapidjson/rapidjson.h - specify the variable RAPIDJSON_INCLUDE_DIR to point to the directory <RapidJSON>/include"
  RAPIDJSON_INCLUDE_DIR)
