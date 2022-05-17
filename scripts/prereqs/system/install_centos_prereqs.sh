#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2019-2021 Omics Data Automation, Inc.
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

install_devtoolset() {
  echo "Installing devtoolset"
  if [[ $CENTOS_VERSION -ge 8 ]]; then
    yum -y makecache &&
      yum -y grouplist &&
      yum -y groupinstall "Development Tools"
  else
    yum install -y centos-release-scl &&
      yum install -y devtoolset-7 &&
      echo "source /opt/rh/devtoolset-7/enable" >> $PREREQS_ENV
  fi
}

install_openjdk() {
  echo "Installing openjdk" &&
  yum install -y -q java-1.8.0-openjdk-devel &&
  echo "export JRE_HOME=/usr/lib/jvm/jre" >> $PREREQS_ENV
}

install_mpi() {
  yum install -y -q mpich-devel &&
  echo "export PATH=/usr/lib64/mpich/bin:\$PATH" >> $PREREQS_ENV &&
  echo "export LD_LIBRARY_PATH=/usr/lib64:/usr/lib:\$LD_LIBRARY_PATH" >> $PREREQS_ENV
}

install_csv() {
  if [[ $CENTOS_VERSION -lt 8 ]]; then 
    yum install -y -q libcsv libcsv-devel &&
      yum install -y -q csv
  fi
}

install_system_prerequisites() {
  yum install -y -q deltarpm
  yum update -y -q &&
    yum install -y -q sudo &&
    yum install -y -q which wget git make &&
    install_mpi &&
    install_devtoolset &&
    install_openjdk &&
    yum install -y -q autoconf automake libtool unzip &&
    yum update -y -q autoconf &&
    yum install -y -q epel-release &&
    yum install -y zlib-devel &&
    yum install -y -q openssl-devel &&
    yum install -y -q cmake3 &&
    yum install -y -q libuuid libuuid-devel &&
    yum install -y -q patch &&
    yum install -y -q zstd &&
    yum install -y -q curl libcurl-devel &&
    install_csv
}

install_nobuild_prerequisites() {
  yum install -y -q deltarpm
  yum update -y -q &&
    install_mpi &&
    yum install -y -q epel-release &&
    yum install -y zlib-devel &&
    yum install -y -q openssl-devel &&
    yum install -y -q cmake3 &&
    yum install -y -q patch &&
    yum install -y -q libuuid libuuid-devel &&
    yum install -y -q zstd &&
    yum install -y -q curl libcurl-devel
}
