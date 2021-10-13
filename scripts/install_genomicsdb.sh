#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2020-2021 Omics Data Automation, Inc.
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
#

set -e

if [[ -z $DOCKER_BUILD ]]; then
	echo "This script is meant to be used with Docker. Exiting."
	exit 1
fi

GENOMICSDB_USER=${1:-genomicsdb}
GENOMICSDB_BRANCH=${2:-develop}
GENOMICSDB_INSTALL_DIR=${3:-/usr/local}
BUILD_DISTRIBUTABLE_LIBRARY=${4:-false}
ENABLE_BINDINGS=$5

GENOMICSDB_USER_DIR=`eval echo ~$GENOMICSDB_USER`
GENOMICSDB_DIR=$GENOMICSDB_USER_DIR/GenomicsDB
echo GENOMICSDB_DIR=$GENOMICSDB_DIR

CMAKE=`which cmake3`
if [[ -z $CMAKE ]]; then
  echo "cmake3 not found. Cannot continue"
  exit 1
fi

if [[ $ENABLE_BINDINGS == *java* ||  $BUILD_DISTRIBUTABLE_LIBRARY == true ]]; then
	BUILD_JAVA=true
else
	BUILD_JAVA=false
fi

# Autoconf version 2.73 is the highest version supported by yum install on centos 6,
# and 2.73 does not support m4_esyscmd_s in configure.ac.
# obviously htslib is not supporting centos 6! so hardcoding something for now.
repair_htslib() {
	if [[ -f /etc/centos-release ]]; then
		if grep -q "release 6" /etc/centos-release; then
			pushd dependencies/htslib &&
				VERSION=`./version.sh` &&
				sed -i s/m4_esyscmd_s\(.*version.*\)/[$VERSION]/g configure.ac &&
				popd
		fi
	fi
}

build_genomicsdb() {
	. /etc/profile &&
	git clone https://github.com/GenomicsDB/GenomicsDB -b ${GENOMICSDB_BRANCH} $GENOMICSDB_DIR &&
	pushd $GENOMICSDB_DIR &&
	git submodule update --recursive --init &&
	repair_htslib &&
	echo "Building GenomicsDB" &&
	mkdir build &&
	pushd build &&
	echo "	$CMAKE .. -DCMAKE_INSTALL_PREFIX=$GENOMICSDB_INSTALL_DIR -DBUILD_DISTRIBUTABLE_LIBRARY=$BUILD_DISTRIBUTABLE_LIBRARY -DBUILD_JAVA=$BUILD_JAVA" &&
	$CMAKE .. -DCMAKE_INSTALL_PREFIX=$GENOMICSDB_INSTALL_DIR -DBUILD_DISTRIBUTABLE_LIBRARY=$BUILD_DISTRIBUTABLE_LIBRARY -DBUILD_JAVA=$BUILD_JAVA && make && make install &&
	popd &&
	echo "Building GenomicsDB DONE" &&
	popd
}

setup_genomicsdb_env() {
	GENOMICSDB_ENV=/etc/profile.d/genomicsdb_env.sh
	echo "export GENOMICSDB_HOME=$GENOMICSDB_INSTALL_DIR" > $GENOMICSDB_ENV
	if [[ $GENOMICSDB_INSTALL_DIR != "/" ]] &&  [[ $GENOMICSDB_INSTALL_DIR != "/usr" ]]; then
		echo "export PATH=\$GENOMICSDB_HOME/bin:\$PATH" >> $GENOMICSDB_ENV
		echo "export LD_LIBRARY_PATH=\$GENOMICSDB_HOME/lib:\$LD_IBRARY_PATH" >> $GENOMICSDB_ENV
		echo "export C_INCLUDE_PATH=\$GENOMICSDB_HOME/include" >> $GENOMICSDB_ENV
	fi
	if [[ $ENABLE_BINDINGS == *java* ]]; then
			echo "export CLASSPATH=\$GENOMICSDB_HOME/bin:\$CLASSPATH" >> $GENOMICSDB_ENV
	fi
	. /etc/profile
}

build_genomicsdb &&
setup_genomicsdb_env &&
echo "$GENOMICSDB_USER:$GENOMICSDB_USER" | chpasswd
chown -R $GENOMICSDB_USER:genomicsdb $GENOMICSDB_USER_DIR
