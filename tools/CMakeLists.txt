include_directories("include")

# CMake/GCC doesn't add "-pie" by default for executables (CMake issue #14983)
if(NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pie")
endif()

if(NOT DISABLE_MPI)
  set(deplist src/create_genomicsdb_workspace.cc src/gt_mpi_gather.cc src/vcf2genomicsdb_init.cc src/vcf2genomicsdb.cc src/vcfdiff.cc src/vcf_histogram.cc src/consolidate_genomicsdb_array.cc)
  build_GenomicsDB_executable_dependencies(cli "${deplist}")

  add_test(NAME tools_tests COMMAND ${CMAKE_SOURCE_DIR}/tests/test_tools.sh ${CMAKE_SOURCE_DIR}/tests/inputs/vcfs ${CMAKE_INSTALL_PREFIX})
endif()
