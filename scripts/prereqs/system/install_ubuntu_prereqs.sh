#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2019-2020 Omics Data Automation, Inc.
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

install_jdk8() {
  if type -p javac; then
    JAVAC=java
  elif [[ -n $JAVA_HOME ]] && [[ -x $JAVA_HOME/bin/javac ]];  then
    JAVAC="$JAVA_HOME/bin/java"
  fi
  if [[ ! -z $JAVAC ]]; then
    JDK_VERSION=$($JAVAC -version 2>&1 | awk '/version/{print $2}')
  fi
  if [[ -z $JDK_VERSION || $JDK_VERSION < "1.8" ]]; then
    apt-get -y install openjdk-8-jdk icedtea-plugin &&
      apt-get update -q &&
      pushd $HOME &&
      git clone https://github.com/michaelklishin/jdk_switcher.git &&
      source jdk_switcher/jdk_switcher.sh && jdk_switcher use openjdk8 &&
      echo "export JAVA_HOME=$JAVA_HOME" >> $PREREQS_ENV &&
      popd
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
  fi
}

install_gcc() {
  GCC=`which gcc`
  GPLUSPLUS=`which g++`
  if [[ -z $GCC || -z $GPLUSPLUS ]]; then
    apt-get -y install build-essential software-properties-common
  fi
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

install_system_prerequisites() {
  apt-get update -q &&
    install_package tzdata &&
    install_package install 1 &&
    install_package lcov 1 &&
    install_package mpich &&
    install_package zlib1g-dev &&
    install_package libssl-dev &&
    install_package rsync 1 &&
    install_package libidn11 &&
    install_package uuid-dev &&
    install_package libcurl4-openssl-dev &&
    install_package wget 1 &&
    install_package autoconf 1 &&
    install_package automake 1 &&
    install_package libtool 1 &&
    install_package zip 1 &&
    install_package unzip 1 &&
    install_package curl 1 &&
    install_package git &&
    install_package libcsv-dev &&
    install_package sudo 1 &&
    apt-get update -q &&
    install_gcc &&
    install_jdk8 &&
    install_cmake3 &&
    install_R &&
    apt-get clean &&
    apt-get purge -y &&
    rm -rf /var/lib/apt/lists*
}
