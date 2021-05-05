include(CheckCXXCompilerFlag)
#Check C++ 2011 support
macro (CHECK_AND_SET_STD_CXX_2011_FLAG CXX_STD_2011_FOUND)
  set(BACKUP_CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  #string_view requires C++ 2014 support - most compilers that
  #support string_view also support 2014 (or 201y)
  set(CMAKE_CXX_FLAGS "${BACKUP_CMAKE_CXX_FLAGS} -std=c++14")
  CHECK_CXX_COMPILER_FLAG(-std=c++14 ${CXX_STD_2011_FOUND})
  if(NOT ${CXX_STD_2011_FOUND})
    set(CMAKE_CXX_FLAGS "${BACKUP_CMAKE_CXX_FLAGS} -std=c++1y")
    CHECK_CXX_COMPILER_FLAG(-std=c++1y ${CXX_STD_2011_FOUND})
    if(NOT ${CXX_STD_2011_FOUND})
      set(CMAKE_CXX_FLAGS "${BACKUP_CMAKE_CXX_FLAGS} -std=c++11")
      CHECK_CXX_COMPILER_FLAG(-std=c++11 ${CXX_STD_2011_FOUND})
      if(NOT ${CXX_STD_2011_FOUND})
	set(CMAKE_CXX_FLAGS "${BACKUP_CMAKE_CXX_FLAGS}")
      endif()
    endif()
  endif()
endmacro()
