#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2019-2020 Omics Data Automation, Inc.
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

set -e

install_jdk17() {
  if type -p javac; then
    JAVAC=java
  elif [[ -n $JAVA_HOME ]] && [[ -x $JAVA_HOME/bin/javac ]];  then
    JAVAC="$JAVA_HOME/bin/java"
  fi
  if [[ ! -z $JAVAC ]]; then
    JDK_VERSION=$($JAVAC -version 2>&1 | awk '/version/{print $2}')
  fi
  if [[ -z $JDK_VERSION || $JDK_VERSION < "1.17" ]]; then
    apt-get -y install openjdk-17-jdk
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

install_gcc() {
  apt-get -y install build-essential software-properties-common
}

install_package() {
  if [[ $# -ne 2 ]]; then
    dpkg -s $1 &> /dev/null
  else
    echo "Checking if executable exists"
    which $1
  fi
  if [ $? -ne 0 ]; then
    echo "Installing $1..."
    apt-get install $1
    echo "Installing DONE"
  fi
}

install_R() {
  echo "NYI: install R functionality"
  return 0
  apt-get -y install gnupg &&
  apt-get install -y software-properties-common &&
  apt-get update -q &&
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &&
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' &&
  apt-get -y install libxml2-dev &&
  apt-get update -q &&
  apt-get -y install r-base
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

# Sufficient to build/install the genomicsdb libraries
install_minimum_prerequisites() {
  apt-get update -q &&
    apt install -y tzdata \
        zlib1g-dev \
        libssl-dev \
        rsync \
        libidn11-dev \
        wget \
        autoconf \
        automake \
        libtool \
        zip \
        unzip \
        curl \
        git \
        zstd \
        sudo &&
    if [[ $BUILD_DISTRIBUTABLE_LIBRARY == false ]]; then
      apt install -y uuid-dev \
          libcurl4-openssl-dev
    fi &&       
    apt-get update -q &&
    install_gcc &&
    install_cmake3 &&
    apt-get clean &&
    apt-get purge -y &&
    rm -rf /var/lib/apt/lists*
}

install_system_prerequisites() {
  apt-get update -q &&
    apt install -y tzdata \
       lcov \
       mpich \
       zlib1g-dev \
       libssl-dev \
       rsync \
       libidn11-dev \
       uuid-dev \
       libcurl4-openssl-dev \
       wget \
       autoconf \
       automake \
       libtool \
       zip \
       unzip \
       curl \
       git \
       libcsv-dev \
       zstd \
       sudo &&
    apt-get update -q &&
    install_gcc &&
    install_jdk17 &&
    install_cmake3 &&
    install_R &&
    apt-get clean &&
    apt-get purge -y &&
    rm -rf /var/lib/apt/lists*
}

install_nobuild_prerequisites() {
    apt-get update -q &&
    mkdir -p /usr/share/man/man7 && mkdir -p /usr/share/man/man1 &&
    apt-get -y install --no-install-recommends \
                                    zlib1g \
                                    libbz2-1.0 \
                                    libssl3 \
                                    libgomp1 \
                                    mpich \
                                    zstd \
                                    libcsv3 \
                                    libcurl4 &&
    apt-get clean &&
    apt-get purge -y &&
    rm -rf /var/lib/apt/lists*
}
