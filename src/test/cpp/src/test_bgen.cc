#include <catch2/catch.hpp>

#if(INCLUDE_TEST_BGEN)

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

  // basic - direct
  std::string workspace = ctests_input_dir + "/bgen/workspace_1";
  GenomicsDB gdb_basic(workspace, workspace+"/callset.json", workspace+"/vidmap.json", {"REF", "ALT", "DP", "GT", "ID"});
  gdb_basic.generate_plink("1$1$249250621", SCAN_FULL, {{0,20}}, 7, 1);
  std::string cmd_basic = "diff output.bgen " + ctests_input_dir + "/bgen/output_1.bgen > /dev/null";
  int status = system(cmd_basic.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("Direct - first BGEN test output did not match reference output");
  }

  // test 1 (t0_1_2_combined)
  GenomicsDB gdb(ctests_input_dir + "/bgen/query_1.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb.generate_plink(7, 1);

  std::string cmd = "diff output.bgen " + ctests_input_dir + "/bgen/output_1.bgen > /dev/null";
  status = system(cmd.c_str());
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
  GenomicsDB gdb2(ctests_input_dir + "/bgen/query_2.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb2.generate_plink(1, 1);

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/output_2.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("second BGEN test output did not match reference output");
  }

  // test 3 (merged tcga vcfs)
  GenomicsDB gdb3(ctests_input_dir + "/bgen/query_3.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb3.generate_plink(1, 1);

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/output_3.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("third BGEN test output did not match reference output");
  }
 
  system("rm -f output.bgen");

  // test 4 (PL/GL with one sample empty, standard two pass generation)
  GenomicsDB gdb4(ctests_input_dir + "/bgen/bgen_test_query_1.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb4.generate_plink(1, 1); 

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/bgen_test_output_1_2p.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("fourth BGEN test output did not match reference output");
  }
 
  // test 5 (PL/GL with one sample empty, one pass generation)
  GenomicsDB gdb5(ctests_input_dir + "/bgen/bgen_test_query_1.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb5.generate_plink(1, 1, true);

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/bgen_test_output_1_1p.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("fifth BGEN test output did not match reference output");
  }

  // test 6 (PL/GL, standard two pass generation)
  GenomicsDB gdb6(ctests_input_dir + "/bgen/bgen_test_query_2.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb6.generate_plink(1, 1);

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/bgen_test_output_2_2p.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("sixth BGEN test output did not match reference output");
  }

  // test 7 (PL/GL, one pass generation, verbose, progress interval)
  GenomicsDB gdb7(ctests_input_dir + "/bgen/bgen_test_query_2.json", GenomicsDB::JSON_FILE, "", my_world_mpi_rank);
  gdb7.generate_plink(1, 1, true, true, 1000);

  cmd = "diff output.bgen " + ctests_input_dir + "/bgen/bgen_test_output_2_1p.bgen > /dev/null";
  status = system(cmd.c_str());
  if(WEXITSTATUS(status)) {
    throw std::runtime_error("seventh BGEN test output did not match reference output");
  }
}

#endif //INCLUDE_TEST_BGEN
