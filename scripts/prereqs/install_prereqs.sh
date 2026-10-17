#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2019-2022 Omics Data Automation, Inc.
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

# Arguments to the script
#  $1 - 'full' if build prerequisites should be installed, 'release' if only runtime prerequisites should be installed

OPENSSL_VERSION=${OPENSSL_VERSION:-3.0.12}
MAVEN_VERSION=${MAVEN_VERSION:-3.9.2}
BUILD_TYPE=${1:-"full"}

INSTALL_OPENSSL=${INSTALL_OPENSSL:-false}

# Check for the following overriding env variables
#    $INSTALL_PREFIX allows for dependencies maven/protobuf/etc. that are built to be installed to $INSTALL_PREFIX for user installs
#    $PREREQS_ENV will set up file that can be sourced to set up the ENV for building GenomicsDB
if [[ `uname` == "Darwin" || `id -u` -ne 0 ]]; then
  INSTALL_PREFIX=${INSTALL_PREFIX:-$HOME/genomicsdb_prereqs}
  MAVEN_INSTALL_PREFIX=$INSTALL_PREFIX
  PREREQS_ENV=${PREREQS_ENV:-$HOME/genomicsdb_prereqs.sh}
  echo "GenomicsDB dependencies(e.g. maven, protobuf, etc. that are built from source will be installed to \$INSTALL_PREFIX=$INSTALL_PREFIX"
  echo "GenomicsDB environment will be persisted to \$PREREQS_ENV=$PREREQS_ENV"
  echo "Sleeping for 10 seconds. Ctrl-C if you want to change or setup \$INSTALL_PREFIX or \$PREREQS_ENV"
  sleep 10
else
  INSTALL_PREFIX=${INSTALL_PREFIX:-/usr/local}
  MAVEN_INSTALL_PREFIX=/opt
  PREREQS_ENV=${PREREQS_ENV:-/etc/profile.d/genomicsdb_prereqs.sh}
fi
if [ -f $PREREQS_ENV ]; then
  echo "Found PREREQS_ENV=$PREREQS_ENV from a previous run"
  source $PREREQS_ENV
fi
touch $PREREQS_ENV

################################# Should not have to change anything below ############################

CENTOS_VERSION=0
WGET_NO_CERTIFICATE=""
PARENT_DIR="$(dirname $0)"

# $1 - path variable name
# $2 - path variable value
add_to_env() {
  SEP=":"
  if [[ $1 != *PATH ]]; then
    echo "export $1=$2" >> $PREREQS_ENV
  elif grep -q -w $1 $PREREQS_ENV; then
    VALUE=`grep -w $1 $PREREQS_ENV`
    if [[ $VALUE == *$2* ]]; then
      return 0
    fi
    NEW_VALUE="export $1=$2$SEP${VALUE/export $1=/}"
    sed -i "s|$VALUE\$||g" $PREREQS_ENV
    echo $NEW_VALUE >> $PREREQS_ENV
  elif [[ -z $1 ]]; then
    echo "export $1=$2" >> $PREREQS_ENV
  else
    echo "export $1=$2$SEP\$$1" >> $PREREQS_ENV
  fi
}

install_mvn() {
  MVN=`which mvn`
  if [ -z $MVN ]; then
    if [ ! -d $MAVEN_INSTALL_PREFIX/apache-maven-$MAVEN_VERSION ]; then
      echo "Installing Maven"
      MAVEN_MAJOR_VERSION=$(cut -d '.' -f 1 <<< $MAVEN_VERSION)
      wget $WGET_NO_CERTIFICATE -nv https://archive.apache.org/dist/maven/maven-$MAVEN_MAJOR_VERSION/$MAVEN_VERSION/binaries/apache-maven-$MAVEN_VERSION-bin.tar.gz -P /tmp &&
        tar xf /tmp/apache-maven-*.tar.gz -C $MAVEN_INSTALL_PREFIX &&
        rm /tmp/apache-maven-*.tar.gz
      echo "Installing Maven DONE"
    fi
    rm -f $MAVEN_INSTALL_PREFIX/maven
    ln -s $MAVEN_INSTALL_PREFIX/apache-maven-$MAVEN_VERSION $MAVEN_INSTALL_PREFIX/maven &&
      MVN=$MAVEN_INSTALL_PREFIX/maven/bin/mvn &&
      test -f $MVN
  fi
  MVN_BIN_DIR=`dirname ${MVN}`
  M2_HOME=`dirname ${MVN_BIN_DIR}`
  test -d $M2_HOME  &&
    add_to_env M2_HOME $M2_HOME &&
    add_to_env PATH "\${M2_HOME}/bin"
}

centos_version() {
  CENTOS_RELEASE_FILE=/etc/centos-release
  if [[ ! -f $CENTOS_RELEASE_FILE ]]; then
    CENTOS_RELEASE_FILE=/etc/redhat-release
    if [[ ! -f $CENTOS_RELEASE_FILE ]]; then
      return 0
    fi
  fi
  if grep -q "release 6" $CENTOS_RELEASE_FILE; then
    CENTOS_VERSION=6
  elif grep -q "release 7" $CENTOS_RELEASE_FILE; then
    CENTOS_VERSION=7
  elif grep -q "release 8" $CENTOS_RELEASE_FILE; then
    CENTOS_VERSION=8
  fi
}

install_os_prerequisites() {
  case `uname` in
    Linux )
      if apt-get --version >/dev/null 2>&1; then
        export DEBIAN_FRONTEND=noninteractive
        source $PARENT_DIR/system/install_ubuntu_prereqs.sh
      else
        source $PARENT_DIR/system/install_centos_prereqs.sh
      fi
      if [[ $1 == "release" ]]; then
        install_nobuild_prerequisites
      else 
        install_system_prerequisites
      fi
      ;;
    Darwin )
      source $PARENT_DIR/system/install_macos_prereqs.sh
      install_minimum_prerequisites
      install_system_prerequisites
      ;;
    * )
      echo "OS=`uname` not supported"
      exit 1
  esac
}

install_prerequisites() {
  PREREQS_ENV=$PREREQS_ENV install_os_prerequisites $1 && echo "Install OS prerequistes successful" &&
    source $PREREQS_ENV
  if [[ $1 == "full" ]]; then
    install_mvn
  fi
}

finalize() {
  echo
  echo "--"
  if [ $1 -ne 0 ]; then
    echo "Error(s) encountered! Please check stderr and correct the errors before proceeding."
  else
    echo "Installing prerequisites DONE. The GenomicsDB env is written out at $PREREQS_ENV for subsequent usage - "
    echo "    Invoke 'source $PREREQS_ENV' to setup environment for building GenomicsDB."
  fi
  return $1
}

centos_version

RC=1
export -f add_to_env
install_prerequisites $BUILD_TYPE &&
  if [[ $INSTALL_OPENSSL == true ]]; then
    echo "Installing openssl $OPENSSL_VERSION"
    OPENSSL_PREFIX=$INSTALL_PREFIX/ssl install_openssl
  fi &&
  RC=0
finalize $RC
