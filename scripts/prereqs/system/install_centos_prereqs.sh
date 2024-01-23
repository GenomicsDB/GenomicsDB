#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2019-2021 Omics Data Automation, Inc.
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
#

install_devtoolset() {
  echo "Installing devtoolset"
  if [[ $CENTOS_VERSION -ge 8 ]]; then
    yum -y makecache &&
      yum -y grouplist &&
      yum -y groupinstall "Development Tools"
  else
    yum install -y centos-release-scl &&
      if [[ $CENTOS_VERSION -ge 6 ]]; then
        yum install -y devtoolset-11 &&
        echo "source /opt/rh/devtoolset-11/enable" >> $PREREQS_ENV
      else
        yum install -y devtoolset-7 &&
        echo "source /opt/rh/devtoolset-7/enable" >> $PREREQS_ENV
      fi
  fi
}

install_cmake3() {
  CMAKE=`which cmake`
  if [[ ! -z $CMAKE ]]; then
    CMAKE_VERSION=$($CMAKE -version | awk '/version/{print $3}')
  fi
  if [[ -z $CMAKE_VERSION || CMAKE_VERSION < "3.6" ]]; then
    echo "Installing cmake..."
    wget -nv https://github.com/Kitware/CMake/releases/download/v3.19.1/cmake-3.19.1-Linux-x86_64.sh -P /tmp &&
      chmod +x /tmp/cmake-3.19.1-Linux-x86_64.sh &&
      /tmp/cmake-3.19.1-Linux-x86_64.sh --prefix=/usr/local --skip-license
    if [ ! -f /usr/local/bin/cmake3 ]; then
      ln -s /usr/local/bin/cmake /usr/local/bin/cmake3
    fi
  fi
}

install_openjdk() {
  if [[ $CENTOS_VERSION -ge 8 ]]; then
    echo "Installing openjdk" &&
      yum install -y -q java-17-openjdk-devel &&
      echo "export JRE_HOME=/usr/lib/jvm/jre" >> $PREREQS_ENV
  else
    curl -O https://download.java.net/java/GA/jdk17.0.2/dfd4a8d0985749f896bed50d7138ee7f/8/GPL/openjdk-17.0.2_linux-x64_bin.tar.gz &&
    tar zxvf openjdk-17.0.2_linux-x64_bin.tar.gz &&
    mv jdk-17.0.2 /usr/local/ &&
    echo "export JAVA_HOME=/usr/local/jdk-17.0.2" >> $PREREQS_ENV
  fi
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

# Sufficient to build/install the genomicsdb libraries
install_minimum_prerequisites() {
   yum install -y -q deltarpm
   yum update -y -q &&
     yum install -y -q epel-release &&
     yum install -y -q which wget git make &&
     install_devtoolset &&
     yum install -y -q autoconf automake libtool unzip &&
     yum install -y zlib-devel &&
     install_cmake3 &&
     yum install -y -q patch &&
     yum install -y -q perl perl-IPC-Cmd
   if [[ $BUILD_DISTRIBUTABLE_LIBRARY == false ]]; then
     yum install -y -q libuuid libuuid-devel &&
       yum install -y -q curl libcurl-devel
   fi
}

install_openssl() {
  if [[ ! -f $OPENSSL_PREFIX/lib/libcrypto.a ]]; then
    echo "Installing OpenSSL"
    pushd /tmp
    wget $WGET_NO_CERTIFICATE https://www.openssl.org/source/openssl-$OPENSSL_VERSION.tar.gz &&
    tar -xvzf openssl-$OPENSSL_VERSION.tar.gz &&
    cd openssl-$OPENSSL_VERSION &&
    CFLAGS=-fPIC ./config -fPIC no-shared --prefix=$OPENSSL_PREFIX --openssldir=$OPENSSL_PREFIX &&
    make && make install && echo "Installing OpenSSL DONE"
    rm -fr /tmp/openssl*
    popd
  fi
  add_to_env OPENSSL_ROOT_DIR $OPENSSL_PREFIX
  add_to_env LD_LIBRARY_PATH $OPENSSL_PREFIX/lib
}

install_system_prerequisites() {
  yum install -y -q deltarpm
  yum update -y -q &&
    yum install -y -q sudo &&
    install_mpi &&
    install_openjdk &&
    yum install -y -q openssl-devel &&
    yum install -y -q zstd &&
    yum install -y -q libuuid libuuid-devel &&
    yum install -y -q curl libcurl-devel
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
