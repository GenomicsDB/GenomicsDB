#include <catch2/catch.hpp>
#include <iostream>
#include <mpi.h>
#include "tiledb_loader_rc.h"

TEST_CASE("rc constructor", "[rc_constr]") {
  std::cout << "test" << std::endl;

  std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

  //int my_world_mpi_rank = 0;
  //MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
  int idx = 0;

  auto loader_json_config_file = std::string(GENOMICSDB_CTESTS_DIR) + "/rc/bed_workspace/loader.json";

  SinglePosition2TileDBLoader loader(loader_json_config_file, idx);
  loader.read_all();
}
