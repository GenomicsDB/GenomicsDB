#!/bin/bash
#
# The MIT License (MIT)
# Copyright (c) 2020 Omics Data Automation Inc.
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
  echo "Usage: ./test_hfile_plugin.sh <vcfs_dir> <install_dir>"
  echo "<vcf_dir> can be a local or remote directory containing the test vcfs -"
  echo "            t0.vcf.gz, t1.vcf.gz and t2.vcf.gz"
  exit 1
fi

if [[ -n $2 ]]; then
  PATH=$2/bin:$PATH
fi

TESTS_DIR=`dirname $0`
if [[ $(uname) == "Darwin" ]]; then
  TEMP_DIR=$(mktemp -d -t test_genomicsdb_hfile_plugin)
else
  TEMP_DIR=$(mktemp -d -t test_genomicsdb_hfile_plugin-XXXXXXXXXX)
fi

LOADER_JSON=$TEMP_DIR/loader.json
CALLSET_MAPPING_JSON=$TEMP_DIR/callset_mapping.json
WORKSPACE=$TEMP_DIR/ws
TEMPLATE_HEADER=$TESTS_DIR/inputs/template_vcf_header.vcf

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

check_rc() {
  if [[ $# -eq 1 ]]; then
    if [[ $1 -ne 0 ]]; then
      die "command returned $1."
    fi
  fi
}

# create callset.json
TEMP_CALLSET_MAPPING_JSON=$TEMP_DIR/temp_callset_mapping.json
cp $TESTS_DIR/inputs/callsets/t0_1_2.json $TEMP_CALLSET_MAPPING_JSON
if [[ -n $1 ]]; then
  if [[ $1 == */ ]]; then
    sed -e "s?inputs/vcfs/?$1?g" $TEMP_CALLSET_MAPPING_JSON > $CALLSET_MAPPING_JSON
  else
    sed -e "s?inputs/vcfs/?$1/?g" $TEMP_CALLSET_MAPPING_JSON > $CALLSET_MAPPING_JSON
  fi
else
  cp $TEMP_CALLSET_MAPPING_JSON $CALLSET_MAPPING_JSON
fi

# create loader.json
cat > $LOADER_JSON  << EOF
{
    "treat_deletions_as_intervals": true,
    "callset_mapping_file": "$CALLSET_MAPPING_JSON",
    "compress_tiledb_array": true,
    "produce_tiledb_array": true,
    "produce_combined_vcf": true,
    "reference_genome": "$TESTS_DIR/inputs/chr1_10MB.fasta.gz",
    "size_per_column_partition": 700,
    "delete_and_create_tiledb_array": true,
    "num_parallel_vcf_files": 1,
    "discard_vcf_index": true,
    "num_cells_per_tile": 3,
    "offload_vcf_output_processing": false,
    "column_partitions": [
        {
            "begin": 0,
            "array_name": "t0_1_2",
            "workspace": "$WORKSPACE"
        }
    ],
    "vcf_header_filename": "$TEMPLATE_HEADER",
    "row_based_partitioning": false,
    "segment_size": 40,
    "do_ping_pong_buffering": false,
    "vid_mapping_file": "$TESTS_DIR/inputs/vid.json"
}
EOF

COMBINED_VCF="$TEMP_DIR/combined.vcf"
check_rc `create_genomicsdb_workspace $WORKSPACE`
check_rc `vcf2genomicsdb $LOADER_JSON > $COMBINED_VCF`

if [[ ! -f $COMBINED_VCF ]]; then
  die "vcf2genomicsdb does not seem to have generated any output"
fi

#TODO: compare the combined vcf with golden_outputs/t0_1_2_combined
ACTUAL_FILESIZE=`wc -c $COMBINED_VCF | awk '{print $1}'`
HEADER_FILESIZE=`wc -c $TEMPLATE_HEADER | awk '{print $1}'`
if [[ $ACTUAL_FILESIZE -le $HEADER_FILESIZE ]]; then
  die "vcf2genomicsdb did not generate the expected vcf"
fi

# Test out vcf2genomicsdb_init/vcf2genomicsdb combo

TEMPLATE_LOADER_JSON=$TEMP_DIR/template_loader.json
cat > $TEMPLATE_LOADER_JSON  << EOF
{
    "treat_deletions_as_intervals": true,
    "compress_tiledb_array": true,
    "produce_tiledb_array": true,
    "produce_combined_vcf": true,
    "reference_genome": "$TESTS_DIR/inputs/chr1_10MB.fasta.gz",
    "size_per_column_partition": 700,
    "delete_and_create_tiledb_array": true,
    "num_parallel_vcf_files": 1,
    "discard_vcf_index": true,
    "num_cells_per_tile": 3,
    "offload_vcf_output_processing": false,
    "row_based_partitioning": false,
    "segment_size": 40,
    "vcf_header_filename": "$TEMPLATE_HEADER",
    "do_ping_pong_buffering": true
}
EOF

SAMPLE_LIST=$TEMP_DIR/sample.list
echo $1/t0.vcf.gz > $SAMPLE_LIST
echo $1/t1.vcf.gz >> $SAMPLE_LIST
echo $1/t2.vcf.gz >> $SAMPLE_LIST
echo "cat sample.list"
cat $SAMPLE_LIST

if [[ $1 == *//* ]]; then
  echo "Processing from samples dir $1"
  check_rc `vcf2genomicsdb_init -w $WORKSPACE -o -S $1 -t $TEMPLATE_LOADER_JSON`
else
  echo "Processing from sample list"
  check_rc `vcf2genomicsdb_init -w $WORKSPACE -o -s $SAMPLE_LIST -t $TEMPLATE_LOADER_JSON`
fi

check_rc `vcf2genomicsdb $WORKSPACE/loader.json > ${COMBINED_VCF}_new`
ACTUAL_FILESIZE_NEW=`wc -c ${COMBINED_VCF}_new | awk '{print $1}'`

# Approximate check rounded to the tenth place as the file sizes are going to be different
if [[ $ACTUAL_FILESIZE_NEW/10 -ne $ACTUAL_FILESIZE/10 ]]; then
  echo "Filesizes : $ACTUAL_FILESIZE_NEW $ACTUAL_FILESIZE"
  die "vcf2genomicsdb with vcf2genomicsdb_init did not generate the expected vcf"
fi

cleanup
