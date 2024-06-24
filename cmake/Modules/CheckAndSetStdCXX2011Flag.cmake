#
# CheckAndSetStdCXX2011Flag.cmake
#
# The MIT License
#
# Copyright (c) 2024 dātma, inc™
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

include(CheckCXXCompilerFlag)
#Check C++ Compiler Support
macro (CHECK_AND_SET_STD_CXX_2011_FLAG CXX_STD_2011_FOUND)
  set(BACKUP_CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  if(BUILD_NANOARROW)
    #std::binary_semaphore requires C++ 2020
    CHECK_CXX_COMPILER_FLAG(-std=c++20 ${CXX_STD_2011_FOUND})
    if(${CXX_STD_2011_FOUND})
      if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 11.3)
        set(CMAKE_CXX_STANDARD 20)
      else()
        set(CMAKE_CXX_FLAGS "${BACKUP_CMAKE_CXX_FLAGS} -std=c++2a")
      endif()
    else()
      message(STATUS_ERROR, "BUILD_NANOARROW uses std::semaphore that requires C++20 compatible compilers")
    endif()
  else()
    #string_view requires C++ 2017 support - most compilers that
    #support string_view also support 2017 (or 201y)
    #set(CMAKE_CXX_FLAGS "${BACKUP_CMAKE_CXX_FLAGS} -std=c++17")
    CHECK_CXX_COMPILER_FLAG(-std=c++17 ${CXX_STD_2011_FOUND})
    if(${CXX_STD_2011_FOUND})
      set(CMAKE_CXX_STANDARD 17)
    endif()
  endif()
  
  if(NOT ${CXX_STD_2011_FOUND})
    message(WARNING "Need at least C++17 for compiling source as is, but it can be modified easily to make it C++11 compliant")
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
  
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endmacro()
