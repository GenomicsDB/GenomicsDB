#!/bin/bash
#
# The MIT License (MIT)
# Copyright (c) 2021-2023 Omics Data Automation Inc.
# Copyright (c) 2023 da̅tma, inc.
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

# Sanity test command line tools -
#    vcf2genomicsdb_init
#    vcf2genomicsdb
#    gt_mpi_gather
# $1 contains the test vcfs - t0.vcf.gz, t1.vcf.gz and t2.vcf.gz
if [[ $# -lt 2 ]]; then
  echo "Usage: ./test_tools.sh <vcfs_dir> <install_dir> <disable_mpi>"
  echo "<vcf_dir> can be a local or remote directory containing the test vcfs -"
  echo "            t0.vcf.gz, t1.vcf.gz and t2.vcf.gz"
  echo "optional <disable_mpi>, if set gt_mpi_gather_tests are not run"
  exit 1
fi

if [[ -n $2 ]]; then
  PATH=$2/bin:$2:$PATH
fi

if [[ $# -ge 3 ]]; then
  DISABLE_MPI=$3
else
  DISABLE_MPI=0
fi

VCFS_DIR=$(cd $1 && pwd)
TEMP_DIR=$(mktemp -d -t test_tools-XXXXXXXXXX)

WORKSPACE=$TEMP_DIR/ws

REFERENCE_GENOME=$VCFS_DIR/../chr1_10MB.fasta.gz

EMPTY=""

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

sanity_check() {
  CMD=$1
  run_command "$CMD" ERR
  run_command "$CMD --help"
  run_command "$CMD --version"
  run_command "$CMD --notanargument" ERR
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

# create_template_loader_json
#    (Optional) $1 tiledb_compression_type
#    (Optional) $2 tiledb_compression_level
create_template_loader_json() {
  TEMPLATE=$TEMP_DIR/template_loader_json_$RANDOM
  cat > $TEMPLATE  << EOF
{
    "treat_deletions_as_intervals": true,
    "compress_tiledb_array": true,
    "produce_combined_vcf": false,
    "produce_tiledb_array": true,
    "size_per_column_partition": 700,
    "delete_and_create_tiledb_array": true,
    "num_parallel_vcf_files": 1,
    "discard_vcf_index": true,
    "num_cells_per_tile": 3,
    "offload_vcf_output_processing": false,
    "row_based_partitioning": false,
    "segment_size": 400,
EOF
  if [[ $# -ge 1 ]]; then
    cat >> $TEMPLATE << EOF
    "tiledb_compression_type": $1,
EOF
  elif [[ $# -ge 2 ]]; then
    cat >> $TEMPLATE << EOF
    "tiledb_compression_level": $2,
EOF
  fi
  cat >> $TEMPLATE << EOF
    "do_ping_pong_buffering": false
}
EOF
}

# create_query_json
#    (Optional) $1 additional json entries
create_query_json() {
  QUERY_JSON=$TEMP_DIR/query_json_$RANDOM
  cat > $QUERY_JSON  << EOF
{
    "workspace": "$WORKSPACE",
    "generate_array_name_from_partition_bounds": true,
    "query_contig_intervals": [{
        "contig": "1",
        "begin": 1
    }],
    "query_row_ranges": [{
        "range_list": [{
            "low" : 0,
            "high": 3
        }]
     }],
EOF
  if [[ $# -ge 1 ]]; then
    cat >> $QUERY_JSON  << EOF
    $1,
EOF
  fi
  cat >> $QUERY_JSON << EOF
    "attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT" ]
}
EOF
}

# create_query_json with intervals
#    $1 is the query interval specification
#    (Optional) $2 additional json entries
create_query_json_with_intervals() {
  QUERY_JSON=$TEMP_DIR/query_json_$RANDOM
  cat > $QUERY_JSON  << EOF
{
    "workspace": "$WORKSPACE",
    "generate_array_name_from_partition_bounds": true,
    $1,
    "query_row_ranges": [{
        "range_list": [{
            "low" : 0,
            "high": 3
        }]
     }],
EOF
  if [[ $# -ge 2 ]]; then
    cat >> $QUERY_JSON  << EOF
    $2,
EOF
  fi
  cat >> $QUERY_JSON << EOF
    "attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PGT" ]
}
EOF
}

# run_command
#    $1 : command to be executed
#    $2 : optional - 0(default) if command should return successfully
#                    any other value if the command should return a failure
run_command() {
  echo $EMPTY > $TEMP_DIR/output
  declare -i EXPECT_NZ
  declare -i GOT_NZ  
  EXPECT_NZ=0
  GOT_NZ=0
  if [[ $# -eq 2 && $2 -ne 0 ]]; then
    EXPECT_NZ=1
  fi
  # Execute the command redirecting all output to $TEMP_DIR/output
  $($1 &> $TEMP_DIR/output)
  retval=$?

  if [[ $retval -ne 0 ]]; then
    GOT_NZ=1
  fi

  if [[ $(($GOT_NZ ^ $EXPECT_NZ)) -ne 0 ]]; then
    cat $TEMP_DIR/output
    die "command '`echo $1`' returned $retval unexpectedly"
  fi
}

time_command() {
  if [[ $(uname) == "Darwin" ]]; then
    DT="date -u +%s"
  else
    DT="date -u +%s%N"
  fi
  START=$($DT)
  run_command "$@"
  END=$($DT)
  TIME=$(($END-$START))
}

OK=0
ERR=1

# Sanity Checks for Tools
sanity_check "vcf2genomicsdb_init"
sanity_check "vcf2genomicsdb"
sanity_check "consolidate_genomicsdb_array"
if [[ $DISABLE_MPI == 0 ]]; then
  sanity_check "gt_mpi_gather"
fi

WORKSPACE=$TEMP_DIR/ws_$RANDOM
run_command "vcf2genomicsdb_init -w $WORKSPACE" ERR

#    $1 actual
#    $2 expected
#    $3 is error message
STATUS=0
assert_true() {
  if [[ $1 -ne $2 ]]; then
    echo "Assertion Failed : $3, actual=$1 expected=$2"
    STATUS=1
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
  # Validate by running vcf2genomicsdb with the generated loader json
  run_command "vcf2genomicsdb --progress $WORKSPACE/loader.json"
}
      
# Basic Tests for vcf2genomicsdb_init
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

# Check with explicit progress update set
export VCF2GENOMICSDB_INIT_PROGRESS_UPDATE_SAMPLE_SIZE=1
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 2 85 24 85 "#16"

# Fields test
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -o -S $SAMPLE_DIR -f GT,DP" 2 85 2 85 "#17"

# Template loader json
create_template_loader_json
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE" 2 85 24 85 "#18"
assert_true $(grep '"segment_size": 400' $WORKSPACE/loader.json | wc -l) 1 "Test #18 segment_size from template loader json was not applied"

# Fail if same field in INFO and FORMAT have different types
create_sample_list inconsistent_DP_t0.vcf.gz
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" ERR

# Try supported compression types/levels,
# see https://github.com/OmicsDataAutomation/TileDB/blob/b338ac9f84f5afde3b083a148d74019f37495fec/core/include/c_api/tiledb_constants.h#L146
TILEDB_COMPRESSION_ZLIB=1
TILEDB_COMPRESSION_ZSTD=2
TILEDB_COMPRESSION_LZ4=3

create_sample_list t0.vcf.gz t1.vcf.gz
create_template_loader_json $TILEDB_COMPRESSION_ZLIB -1
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE" 2 85 24 85 "$20"
create_template_loader_json $TILEDB_COMPRESSION_ZSTD -1
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE" 2 85 24 85 "$21"
create_template_loader_json $TILEDB_COMPRESSION_LZ4 -1
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE" 2 85 24 85 "$22"
# No compression level specified
create_template_loader_json $TILEDB_COMPRESSION_ZLIB
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE" 2 85 24 85 "$23"

# Fail if generate_array_name_from_partition_bounds(not supported yet!) is set
create_template_loader_json
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o -t $TEMPLATE"
sed -i -e 's/"array_name": .*,/"generate_array_name_from_partition_bounds": true,/g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" ERR

# Try running with produce_combined_vcf
create_template_loader_json
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o -t $TEMPLATE"
sed -i -e 's/"produce_combined_vcf": false/"produce_combined_vcf": true/g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" 134
sed -i -e 's|"produce_tiledb_array": true|"reference_genome": "'${REFERENCE_GENOME}'"|g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" OK
sed -i -e 's|"row_based_partitioning": false|"vcf_header_filename": "'${WORKSPACE}'/vcfheader.vcf"|g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" OK
sed -i -e 's|"workspace": "'${WORKSPACE}'"|"workspace": "'${WORKSPACE}'", "vcf_header_filename": "'${WORKSPACE}'/vcfheader.vcf", "vcf_output_filename": "'${WORKSPACE}'/out"|g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" OK

# Test with too small size_per_column_partition
sed -i -e 's/"size_per_column_partition": 700/"size_per_column_partition": 70/g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" OK
sed -i -e 's/"size_per_column_partition": 70/"size_per_column_partition": 0/g' $WORKSPACE/loader.json
run_command "vcf2genomicsdb $WORKSPACE/loader.json" ERR

# Test --progress switch with an interval
create_template_loader_json
run_command "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE"
run_command "vcf2genomicsdb --progress=2 $WORKSPACE/loader.json"
run_command "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE"
run_command "vcf2genomicsdb --progress=10.5s $WORKSPACE/loader.json"
run_command "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE"
run_command "vcf2genomicsdb --progress=.1m $WORKSPACE/loader.json"
run_command "vcf2genomicsdb_init -w $WORKSPACE -S $SAMPLE_DIR -o -t $TEMPLATE"
run_command "vcf2genomicsdb --progress=.001h $WORKSPACE/loader.json"

# Fail with unsupported compression levels
create_template_loader_json -5 -5
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o -t $TEMPLATE"
run_command "vcf2genomicsdb $WORKSPACE/loader.json" ERR

if [[ $DISABLE_MPI == 1 ]]; then
  echo "Build with DISABLE_MPI option set, not running gt_mpi_gather tests"
  exit 0
fi

# Sanity test gt_mpi_gather
create_sample_list t0.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 1 85 24 85 "#30"

create_query_json
run_command "gt_mpi_gather -w "non-existent-workspace" $WORKSPACE/loader.json -j $QUERY_JSON --print-AC" ERR
run_command "gt_mpi_gather -A "non-existent-array" $WORKSPACE/loader.json -j $QUERY_JSON --print-AC" ERR
run_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --print-AC"
run_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --print-calls"

# Test gt_mpi_gather with a partition not imported
run_command "gt_mpi_gather --rank 1 -l $WORKSPACE/loader.json -j $QUERY_JSON" ERR

# Diff results from two gt_mpi_gather invocations
# $1 Results from 1st invocation
# $2 Results from 2nd invocation
# $3 Expected Result(OK or ERR)
# $4 Error Message
diff_gt_mpi_gather_output() {
  sed "1 d" $1 > $1.cut
  sed -i -e '/MD5/d' $1.cut
  sed "1 d" $2 > $2.cut
  sed -i -e '/MD5/d' $2.cut
  DIFF_RESULT=$(diff $1.cut $2.cut)
  if [[ ! -z  $DIFF_RESULT && $3 == OK ]]; then
    echo "gt_mpi_gather results differ : $4"
    echo "First output:"
    cat $1.cut
    echo "Second output:"
    cat $2.cut
    STATUS=1
  elif [[ -z $DIFF_RESULT && $3 == ERR ]]; then
    echo "gt_mpi_gather results should differ but were identical : $4"
    STATUS=1
  fi
}

# $1 : intervals
test_bypass_intersecting_intervals_phase() {
  create_query_json_with_intervals "$1"
  time_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --print-calls"
  DIFF=$TIME
  mv -f $TEMP_DIR/output $TEMP_DIR/output.without_bypass
  create_query_json_with_intervals "$1" "\"bypass_intersecting_intervals_phase\": true"
  time_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --print-calls"
  if [ $DIFF -lt $TIME ]; then
    echo "Bypass intersecting intervals phase took longer than the normal two pass iterator"
    STATUS=1
  fi
}

# Test bypass of intersecting intervals phase in the genomicsdb iterators with no output
INTERVALS="\"query_contig_intervals\": [{\"contig\": \"1\", \"begin\": 1, \"end\":100 }]"
test_bypass_intersecting_intervals_phase "$INTERVALS"
diff_gt_mpi_gather_output $TEMP_DIR/output.without_bypass $TEMP_DIR/output OK "Setting bypass intersecting intervals should not affect gt_mpi_gather for $INTERVALS for t0.vcf.gz"

# Test bypass of intersecting intervals phase in the genomicsdb iterators
INTERVALS="\"query_contig_intervals\": [{\"contig\": \"1\", \"begin\": 1, \"end\":12142 }]"
test_bypass_intersecting_intervals_phase "$INTERVALS"
diff_gt_mpi_gather_output $TEMP_DIR/output.without_bypass $TEMP_DIR/output OK "Setting bypass intersecting intervals should not affect gt_mpi_gather for $INTERVALS for t0.vcf.gz"

# Test bypass of intersecting intervals phase in the genomicsds iterators with two intervals
INTERVALS="\"query_contig_intervals\": [{\"contig\": \"1\", \"begin\": 1, \"end\":12142 }, {\"contig\": \"1\", \"begin\": 17385 }]"
test_bypass_intersecting_intervals_phase "$INTERVALS"
diff_gt_mpi_gather_output $TEMP_DIR/output.without_bypass $TEMP_DIR/output OK "Setting bypass intersecting intervals should not affect gt_mpi_gather for $INTERVALS for t0.vcf.gz"

# Test bypass of intersecting intervals phase with indels
create_sample_list t6_asa.vcf.gz
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o"
run_command "vcf2genomicsdb $WORKSPACE/loader.json"

# Test bypass of intersecting intervals phase with indels with query interval encompassing entire indel
INTERVALS="\"query_contig_intervals\": [{\"contig\": \"1\", \"begin\": 8029500, \"end\": 8029510 }]"
test_bypass_intersecting_intervals_phase "$INTERVALS"
diff_gt_mpi_gather_output $TEMP_DIR/output.without_bypass $TEMP_DIR/output OK "Setting bypass intersecting intervals should not affect gt_mpi_gather results to be different for $INTERVALS for t6_asa.vcf.gz as the entire indel is included in the query interval"

# Test bypass of intersecting intervals phase with indels with query interval overlapping one position in the indel
INTERVALS="\"query_contig_intervals\": [{\"contig\": \"1\", \"begin\": 8029501 }]"
test_bypass_intersecting_intervals_phase "$INTERVALS"
diff_gt_mpi_gather_output $TEMP_DIR/output.without_bypass $TEMP_DIR/output ERR "Setting bypass intersecting intervals should cause gt_mpi_gather to be different for $INTERVALS for t6_asa.vcf.gz as the start of the indel is not included in the query interval"

# Test bypass of intersecting intervals phase with indels with query interval overlapping start of indel
INTERVALS="\"query_contig_intervals\": [{\"contig\": \"1\", \"begin\": 8029500, \"end\": 8029501 }]"
test_bypass_intersecting_intervals_phase "$INTERVALS"
diff_gt_mpi_gather_output $TEMP_DIR/output.without_bypass $TEMP_DIR/output OK "Setting bypass intersecting intervals should not affect gt_mpi_gather results to be different for $INTERVALS for t6_asa.vcf.gz as the start of the indel is included in the query interval"

# Test exception thrown with vcfs having overlapping variants, where the call is neither a deletion, mnv or a reference block
create_sample_list t9_unhandled_overlapping_variants.vcf.gz
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o"
run_command "vcf2genomicsdb -r 12 $WORKSPACE/loader.json"
create_query_json_with_intervals "\"query_contig_intervals\": [{\"contig\": \"13\", \"begin\": 1, \"end\":115169878 }]" "\"reference_genome\" : \"$REFERENCE_GENOME\""
sed -i -e 's/"RAW_MQ",//g' $QUERY_JSON
run_command "gt_mpi_gather --rank 12 -l $WORKSPACE/loader.json -j $QUERY_JSON --produce-Broad-GVCF" ERR

# Fail if the reference genome cannot be parsed by htslib for any reason
create_sample_list t0.vcf.gz
run_command "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o"
run_command "vcf2genomicsdb $WORKSPACE/loader.json"

# Fail when the reference file is not a reference file
create_query_json "\"reference_genome\" : $WORKSPACE/loader.json"
run_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --produce-Broad-GVCF" ERR

# Run with valid reference genome successfully
create_query_json "\"reference_genome\" : \"$REFERENCE_GENOME\""
run_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --produce-Broad-GVCF"

# Test consolidate_genomicsdb_array and gt_mpi_gather after consolidation
run_consolidation_and_check_results() {
  run_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --print-AC"
  mv -f $TEMP_DIR/output $TEMP_DIR/output.no.consolidation
  run_command "$1"
  run_command "gt_mpi_gather -l $WORKSPACE/loader.json -j $QUERY_JSON --print-AC"
  diff_gt_mpi_gather_output  $TEMP_DIR/output.no.consolidation $TEMP_DIR/output $2 OK "After consolidation"
}

ARRAY="1\$1\$249250621"

# Workspace and array specified cannot be null
run_command "consolidate_genomicsdb_array --workspace $EMPTY --array $ARRAY" ERR
run_command "consolidate_genomicsdb_array -w $WORKSPACE -a $EMPTY" ERR

create_sample_list t0.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 1 85 24 85 "#40"
run_consolidation_and_check_results  "consolidate_genomicsdb_array --workspace $WORKSPACE --array $ARRAY" "#41"

create_sample_list t1.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -a" 2 85 24 85 "#42"
run_consolidation_and_check_results  "consolidate_genomicsdb_array -w $WORKSPACE -a $ARRAY -p" "#43"

# Test consolidate_genomicsdb_array with buffer-size
create_sample_list t2.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -a" 3 85 24 85 "#44"
run_consolidation_and_check_results "consolidate_genomicsdb_array -w $WORKSPACE -a $ARRAY -z 1024 -p" "#45"

# Test consolidate_genomicsdb_array with batch-size
create_sample_list t0.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -o" 1 85 24 85 "#46"
create_sample_list t1.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -a" 2 85 24 85 "#47"
create_sample_list t2.vcf.gz
run_command_and_check_results "vcf2genomicsdb_init -w $WORKSPACE -s $SAMPLE_LIST -a" 3 85 24 85 "#48"
run_consolidation_and_check_results "consolidate_genomicsdb_array -w $WORKSPACE -a $ARRAY -z 1024 -b 2 -p" "#49"

cleanup
exit $STATUS
