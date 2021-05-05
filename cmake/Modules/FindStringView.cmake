# Once done this will define
# STRING_VIEW_FOUND if string_view found

include(CheckCXXSourceCompiles)
set(STRING_VIEW_SRC
  "
  #include <string_view>
  int main() {
    std::string_view(\"qwerty\");
    return 0;
  }
  "
  )
set(STRING_VIEW_EXPERIMENTAL_SRC
  "
  #include <experimental/string_view>
  int main() {
    std::experimental::string_view(\"qwerty\");
    return 0;
  }
  "
  )
check_cxx_source_compiles("${STRING_VIEW_SRC}" STRING_VIEW_FOUND)
if(NOT STRING_VIEW_FOUND)
  check_cxx_source_compiles("${STRING_VIEW_EXPERIMENTAL_SRC}" STRING_VIEW_EXPERIMENTAL_FOUND)
  if(STRING_VIEW_EXPERIMENTAL_FOUND)
    message(STATUS "experimental/string_view found")
    add_definitions(-DSTRING_VIEW_EXPERIMENTAL_FOUND=1)
  endif()
else()
  message(STATUS "string_view found")
  add_definitions(-DSTRING_VIEW_FOUND=1)
endif()
