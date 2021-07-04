# Determine compiler flags for fmt library
# Once done this will define
# FMT_FOUND - fmt found

if(USING_FMT_BUNDLED_WITH_SPDLOG)
  #fmt headers are in dependencies/spdlog/include/spdlog/fmt/bundled - source files will include "fmt/bundled/format.h"
  find_path(FMT_INCLUDE_DIR NAMES spdlog/fmt/bundled/core.h HINTS "${FMT_ROOT_DIR}/include" "${FMT_ROOT_DIR}")
else()
  find_path(FMT_INCLUDE_DIR NAMES fmt/core.h HINTS "${FMT_ROOT_DIR}/include" "${FMT_ROOT_DIR}")
endif()

if(BUILD_DISTRIBUTABLE_LIBRARY)
  find_library(FMT_LIBRARY NAMES libfmt.a fmt HINTS "${FMT_ROOT_DIR}/lib64" "${FMT_ROOT_DIR}/lib" "${FMT_ROOT_DIR}")
else()
  find_library(FMT_LIBRARY NAMES fmt HINTS "${FMT_ROOT_DIR}/lib64" "${FMT_ROOT_DIR}/lib" "${FMT_ROOT_DIR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(fmt "Could not find fmt library headers ${DEFAULT_MSG}" FMT_INCLUDE_DIR)
