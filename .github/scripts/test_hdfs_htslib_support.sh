#
# The MIT License (MIT)
# Copyright (c) 2020 Omics Data Automation, Inc.
# Copyright (c) 2023 dātma, inc™
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

#!/bin/bash

set -e

# VCFs are copied into hdfs path
# Name Node is $1 - optional
# hdfs path is $2 - optional

NAME_NODE=${1:-"hdfs://localhost:9000"}
HDFS_PATH=${2:-"/home/hadoop/input/vcfs"}

VCF_DIRS=$GITHUB_WORKSPACE/tests/inputs/vcfs

hdfs dfs -mkdir -p $NAME_NODE$HDFS_PATH
hdfs dfs -copyFromLocal $VCF_DIRS/t0.vcf.gz $NAME_NODE$HDFS_PATH
hdfs dfs -copyFromLocal $VCF_DIRS/t0.vcf.gz.tbi $NAME_NODE$HDFS_PATH
hdfs dfs -copyFromLocal $VCF_DIRS/t1.vcf.gz $NAME_NODE$HDFS_PATH
hdfs dfs -copyFromLocal $VCF_DIRS/t1.vcf.gz.tbi $NAME_NODE$HDFS_PATH
hdfs dfs -copyFromLocal $VCF_DIRS/t2.vcf.gz $NAME_NODE$HDFS_PATH
hdfs dfs -copyFromLocal $VCF_DIRS/t2.vcf.gz.tbi $NAME_NODE$HDFS_PATH

echo "$GITHUB_WORKSPACE/tests/test_hfile_plugin.sh $NAME_NODE$HDFS_PATH $CMAKE_INSTALL_PREFIX"
$GITHUB_WORKSPACE/tests/test_hfile_plugin.sh $NAME_NODE$HDFS_PATH $CMAKE_INSTALL_PREFIX

