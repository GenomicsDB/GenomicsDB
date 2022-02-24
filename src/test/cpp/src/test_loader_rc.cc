#include <catch2/catch.hpp>
#include <iostream>
#include "tiledb_loader_rc.h"

TEST_CASE("rc constructor", "[rc_constr]") {
  std::cout << "test loader_rc" << std::endl;

  std::string ctests_input_dir(GENOMICSDB_CTESTS_DIR);

  //int my_world_mpi_rank = 0;
  //MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
  int idx = 0;

  auto workspace = std::string(GENOMICSDB_CTESTS_DIR) + "rc/bed_workspace";

  //auto loader_json_config_file = std::string(GENOMICSDB_CTESTS_DIR) + "rc/bed_workspace/loader.json";
  auto loader_json = workspace + "/loader.json";
  auto gtf_file = std::string(GENOMICSDB_CTESTS_DIR) + "rc/small.gtf";

  SinglePosition2TileDBLoader loader(loader_json, idx, gtf_file);
  loader.read_all();

  std::string array_name = "1$1$249250621";

  GenomicsDBTranscriptomics gdb(workspace); // TODO pass callset/vid mapper etc
  gdb.query_variant_calls(array_name);

  GenomicsDBVariantCallProcessor processor;
  processor.initialize(gdb.get_genomic_field_types());

  auto retval = gdb.query_variant_calls(processor, array_name);
  std::string test_output = GenomicsDBTranscriptomics::print_transcriptomics_cells(retval);

  std::ifstream gold_file(std::string(GENOMICSDB_CTESTS_DIR) + "rc/gold_output");
  std::stringstream ss;
  ss << gold_file.rdbuf();
  std::string gold_output = ss.str();

  CHECK(test_output == gold_output);
}
