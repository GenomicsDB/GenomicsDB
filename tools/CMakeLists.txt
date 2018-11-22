include_directories("include")

# CMake doesn't add "-pie" by default for executables (CMake issue #14983)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pie")

if(NOT DISABLE_MPI)
    build_GenomicsDB_executable(create_tiledb_workspace)
    build_GenomicsDB_executable(gt_mpi_gather)
    build_GenomicsDB_executable(vcf2tiledb)
    build_GenomicsDB_executable(vcfdiff)
    build_GenomicsDB_executable(vcf_histogram)
    build_GenomicsDB_executable(consolidate_tiledb_array)
endif()
