# Determine compiler flags for fmt library
# Once done this will define
# FMT_FOUND - fmt found

find_path(FMT_INCLUDE_DIR NAMES fmt/core.h HINTS "${FMT_PREFIX_DIR}/include" "${FMT_PREFIX_DIR}")

if(BUILD_DISTRIBUTABLE_LIBRARY)
  find_library(FMT_LIBRARY NAMES libfmt.a fmt HINTS "${FMT_PREFIX_DIR}/lib64" "${FMT_PREFIX_DIR}/lib" "${FMT_PREFIX_DIR}")
else()
  find_library(FMT_LIBRARY NAMES fmt HINTS "${FMT_PREFIX_DIR}/lib64" "${FMT_PREFIX_DIR}/lib" "${FMT_PREFIX_DIR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(fmt "Could not find libuuid headers and/or libraries ${DEFAULT_MSG}" FMT_INCLUDE_DIR FMT_LIBRARY)
