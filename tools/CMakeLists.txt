include_directories("include")

# CMake/GCC doesn't add "-pie" by default for executables (CMake issue #14983)
if(NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pie")
endif()

if(NOT DISABLE_MPI)
  build_GenomicsDB_executable(create_genomicsdb_workspace)
  build_GenomicsDB_executable(gt_mpi_gather)
  build_GenomicsDB_executable(vcf2genomicsdb_init)
  build_GenomicsDB_executable(vcf2genomicsdb)
  build_GenomicsDB_executable(vcfdiff)
  build_GenomicsDB_executable(vcf_histogram)
  build_GenomicsDB_executable(consolidate_genomicsdb_array)

  add_test(NAME tools_tests COMMAND ${CMAKE_SOURCE_DIR}/tests/test_tools.sh ${CMAKE_SOURCE_DIR}/tests/inputs/vcfs ${CMAKE_INSTALL_PREFIX})
endif()
