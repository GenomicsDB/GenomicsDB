#include <catch2/catch.hpp>
#include <iostream>
#include "variant_query_config.h"
#include "genomicsdb.h"
#include <mpi.h>
#include <stdlib.h>

TEST_CASE("test bgen", "[bgen]") {
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  int procs = 1;
  int my_world_mpi_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);

  std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

  // test 1 (t0_1_2_combined)
  VariantQueryConfig query_config;
  query_config.read_from_file(ctests_input_dir + "/bgen/query_1.json", my_world_mpi_rank);

  std::string array_name = query_config.get_array_name(my_world_mpi_rank);

  GenomicsDB gdb(query_config.get_workspace(my_world_mpi_rank),
                 query_config.get_callset_mapping_file(),
                 query_config.get_vid_mapping_file(),
                 query_config.get_reference_genome());
  gdb.generate_plink(array_name, &query_config, 7, 1);

  std::string cmd = "diff output.bgen " + ctests_input_dir + "/bgen/output_1.bgen > /dev/null";
  int status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("first BGEN test output did not match reference output");
  }
  cmd = "diff output.bed " + ctests_input_dir + "/bgen/output_1.bed > /dev/null"; // compare beds
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("first BED test output did not match reference output");
  }
  cmd = "diff output.bim " + ctests_input_dir + "/bgen/output_1.bim > /dev/null"; // compare bim
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("first BIM test output did not match reference output");
  }
  cmd = "diff output.fam " + ctests_input_dir + "/bgen/output_1.fam > /dev/null"; // compare fams
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("first FAM test output did not match reference output");
  }
  cmd = "diff output.tped " + ctests_input_dir + "/bgen/output_1.tped > /dev/null"; // compare beds
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("first tped test output did not match reference output");
  }

  // test 2 (min_PL_spanning_deletion)
  
  VariantQueryConfig query_config2;
  query_config2.read_from_file(ctests_input_dir + "/bgen/query_2.json", my_world_mpi_rank);

  array_name = query_config2.get_array_name(my_world_mpi_rank);

  GenomicsDB gdb2(query_config2.get_workspace(my_world_mpi_rank),
                 query_config2.get_callset_mapping_file(),
                 query_config2.get_vid_mapping_file(),
                 query_config2.get_reference_genome());
  gdb2.generate_plink(array_name, &query_config2, 1, 1);

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/output_2.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("second BGEN test output did not match reference output");
  }

  // test 3 (merged tcga vcfs)

  VariantQueryConfig query_config3;
  query_config3.read_from_file(ctests_input_dir + "/bgen/query_3.json", my_world_mpi_rank);

  array_name = query_config3.get_array_name(my_world_mpi_rank);

  GenomicsDB gdb3(query_config3.get_workspace(my_world_mpi_rank),
                 query_config3.get_callset_mapping_file(),
                 query_config3.get_vid_mapping_file(),
                 query_config3.get_reference_genome());
  gdb3.generate_plink(array_name, &query_config3, 1, 1);

  //std::vector<std::string> types = {"bgen", "bed", "tped", "bim", "fam"};

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/output_3.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("third BGEN test output did not match reference output");
  }
 
  system("rm -f output.bgen");
}
