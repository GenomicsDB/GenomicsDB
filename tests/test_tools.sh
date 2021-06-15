#!/bin/bash
#
# The MIT License (MIT)
# Copyright (c) 2021 Omics Data Automation Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# $1 contains the test vcfs - t0.vcf.gz, t1.vcf.gz and t2.vcf.gz
if [[ $# -ne 2 ]]; then
  echo "Usage: ./test_tools.sh <vcfs_dir> <install_dir>"
  echo "<vcf_dir> can be a local or remote directory containing the test vcfs -"
  echo "            t0.vcf.gz, t1.vcf.gz and t2.vcf.gz"
  exit 1
fi

if [[ -n $2 ]]; then
  PATH=$2/bin:$2:$PATH
fi

VCFS_DIR=$(cd $1 && pwd)
TEMP_DIR=$(mktemp -d -t test_tools-XXXXXXXXXX)

WORKSPACE=$TEMP_DIR/ws
LOADER_JSON=$TEMP_DIR/ws/loader.json
CALLSET_MAPPING_JSON=$TEMP_DIR/ws/callset_mapping.json
TEMPLATE_HEADER=$TESTS_DIR/ws/vcf_header.vcf

cleanup() {
  rm -fr $TEMP_DIR
}

die() {
  cleanup
  if [[ $# -eq 1 ]]; then
    echo $1
  fi
  exit 1
}

create_sample_list() {
  SAMPLE_DIR=$TEMP_DIR/samples_$RANDOM
  rm -fr $SAMPLE_DIR
  mkdir $SAMPLE_DIR
  SAMPLE_LIST=$TEMP_DIR/sample_list_$RANDOM
  rm -f $SAMPLE_LIST
  for sample in $@; do
    cp $VCFS_DIR/$sample* $SAMPLE_DIR
    echo $SAMPLE_DIR/$sample >> $SAMPLE_LIST
  done
}

create_interval_list() {
  INTERVAL_LIST=$TEMP_DIR/interval_list_$RANDOM
  rm -f $INTERVAL_LIST
  for interval in $@; do
    echo $interval >> $INTERVAL_LIST
  done
}

create_template_loader_json() {
  TEMPLATE=$TEMP_DIR/template_loader_json_$RANDOM
  cat > $TEMPLATE  << EOF
{
    "treat_deletions_as_intervals": true,
    "compress_tiledb_array": true,
    "produce_tiledb_array": true,
    "size_per_column_partition": 700,
    "delete_and_create_tiledb_array": true,
    "num_parallel_vcf_files": 1,
    "discard_vcf_index": true,
    "num_cells_per_tile": 3,
    "offload_vcf_output_processing": false,
    "row_based_partitioning": false,
    "segment_size": 400,
    "do_ping_pong_buffering": false
}
EOF
}

# run_command
#    $1 : command to be executed
#    $2 : optional - 0(default) if command should return successfully
#                    any other value if the command should return a failure
run_command() {
  declare -i EXPECTED_RC
  EXPECTED_RC=0
  if [[ $# -eq 2 && $2 -ne 0 ]]; then
    EXPECTED_RC=255
  fi
  # Execute the command redirecting all output to $TEMP_DIR/output
  $($1 &> $TEMP_DIR/output)
  if [[ $? -ne $EXPECTED_RC ]]; then
    cat $TEMP_DIR/output
    die "command '`echo $1`' returned $EXPECTED_RC unexpectedly"
  fi
}

# Sanity Checks
run_command "vcf2genomicsdb_init" 1
run_command "vcf2genomicsdb_init --help"
run_command "vcf2genomicsdb_init --version"
run_command "vcf2genomicsdb_init --notanargument" 1

WORKSPACE=$TEMP_DIR/ws_$RANDOM
run_command "vcf2genomicsdb_init -w $WORKSPACE" 1

#    $1 actual
#    $2 expected
#    $3 is error message
STATUS=0
assert_true() {
  if [[ $1 -ne $2 ]]; then
    echo "Assertion Failed : $3, actual=$1 expected=$2"
    $STATUS=1
  fi
}

#    $1 - command
#    $2 - number of samples in callsets.json
#    $3 - number of partitions in loader.json
#    $4 - number of fields in vidmap.json
#    $5 - number of contigs in vidmap.json
#    $6 - test #
run_command_and_check_results() {
  run_command "$1" 0
  declare -i n_samples
  declare -i n_partitions
  declare -i n_fields
  declare -i n_contigs
  n_samples=$(grep sample_name $WORKSPACE/callset.json | wc -l)
  n_partitions=$(grep array_name $WORKSPACE/loader.json | wc -l)
  n_fields=$(grep vcf_field_class $WORKSPACE/vidmap.json | wc -l)
  n_contigs=$(grep tiledb_column_offset $WORKSPACE/vidmap.json | wc -l)
  assert_true $n_samples $2 "Test $6 Number of samples in callset.json"
  assert_true $n_partitions $3 "Test $6 Number of partitions in loader.json"
  assert_true $n_fields $4 "Test $6 Number of fields in vidmap.json"
  assert_true $n_contigs $5 "Test $6 Number of contigs in vidmap.json"
}

ERR=1
      
# Basic Tests
create_sample_list t0.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST" 1 85 24 85 "#1"
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST" ERR
run_command "vcf2genomicsdb_init -w $WORKSPACE -s non-existent-file" ERR
run_command "vcf2genomicsdb_init -w $WORKSPACE -S non-existent-dir" ERR
NO_SAMPLE_LIST=$TEMP_DIR/no_samples_$RANDOM
touch $NO_SAMPLE_LIST
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $NO_SAMPLE_LIST" ERR
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 1 85 24 85 "#2"
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o" 1 85 24 85 "#3"

# Partition Tests
# single partition -n 0
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -n 0"  1 1 24 85 "#4"
# at least 10 partitions -n 10, will still create one partition per chromosome/contig - 85 contigs in all
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -n 10"  1 85 24 85 "#4"
# at least 100 partitions -n 100
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -n 100"  1 168 24 85 "#5"
# size of partitions -z 1000
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -z 10000000"  1 376 24 85 "#6"
# -z 1000 overrides -n 100
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -n 100 -z 10000000"  1 376 24 85 "#7"
# merge small contigs
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -m"  1 25 24 85 "#8"
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -m -n 100"  1 108 24 85 "#9"
# with interval list
create_interval_list 1:1-100000 1:100001-1000000
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -i $INTERVAL_LIST"  1 2 24 85 "#10"
# -n 0 overrides -i <interval_list>
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -i $INTERVAL_LIST -n 0"  1 1 24 85 "#11"

# Combined vcfs Test
create_sample_list t0_1_2_combined.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 3 85 18 85 "#12"

# Append/Incremental update Test
create_sample_list t0.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 1 85 24 85 "#13"
create_sample_list t1.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -a" 2 85 24 85 "#14"
create_sample_list t0.vcf.gz t1.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -a" 2 85 24 85 "#15"

# Fields test
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -o -S $SAMPLE_DIR -f GT,DP" 2 85 2 85 "#16"

# Template loader json
create_template_loader_json
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE" 2 85 24 85 "#17"
assert_true $(grep '"segment_size": 400' $WORKSPACE/loader.json | wc -l) 1 "Test #16 segment_size from template loader json was not applied"

# Validate by running vcf2genomicsdb with the generated loader json
vcf2genomicsdb -r 1 $WORKSPACE/loader.json

cleanup
exit $STATUS
