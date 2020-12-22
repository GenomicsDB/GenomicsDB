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

# Check for the following overriding env variables
#    $INSTALL_PREFIX allows for dependencies maven/protobuf/etc. that are built to be installed to $INSTALL_PREFIX for user installs
#    $PREREQS_ENV will set up file that can be sourced to set up the ENV for building GenomicsDB
if [[ `uname` == "Darwin" || `id -u` -ne 0 ]]; then
  INSTALL_PREFIX=${INSTALL_PREFIX:-$HOME/genomicsdb}
  MAVEN_INSTALL_PREFIX=INSTALL_PREFIX
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
  rm -f $PREREQS_ENV
fi
touch $PREREQS_ENV

BUILD_DISTRIBUTABLE_LIBRARY=${1:-false}

OPENSSL_VERSION=1.0.2o
MAVEN_VERSION=3.6.3

################################# Should not have to change anything below ############################

CENTOS_VERSION=0
PARENT_DIR="$(dirname $(python -c "import os; import sys; print(os.path.realpath(sys.argv[1]))" $0))"

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
      wget -nv https://www-us.apache.org/dist/maven/maven-3/$MAVEN_VERSION/binaries/apache-maven-$MAVEN_VERSION-bin.tar.gz -P /tmp &&
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

PROTOBUF_PREFIX=$INSTALL_PREFIX
install_protobuf() {
  if [ ! -f $PROTOBUF_PREFIX/bin/protoc ]; then 
    echo "Installing Protobuf"
    pushd /tmp
    if [[ $BUILD_DISTRIBUTABLE_LIBRARY == true ]]; then
      wget -nv https://github.com/protocolbuffers/protobuf/releases/download/v3.0.0-beta-1/protobuf-cpp-3.0.0-beta-1.zip &&
        unzip protobuf-cpp-3.0.0-beta-1.zip &&
        cp $PARENT_DIR/protobuf-v3.0.0-beta-1.autogen.sh.patch protobuf-3.0.0-beta-1/autogen.sh &&
        mv protobuf-3.0.0-beta-1 protobuf
    else
      git clone -b 3.0.x https://github.com/google/protobuf.git
    fi
    pushd protobuf &&
      ./autogen.sh &&
      ./configure --prefix=$INSTALL_PREFIX --with-pic &&
      make -j4 && make install &&
        echo "Installing Protobuf DONE"
    popd
    rm -fr /tmp/protobuf*
    popd
  fi
  add_to_env PATH $PROTOBUF_PREFIX/bin &&
    add_to_env LD_LIBRARY_PATH  $PROTOBUF_PREFIX/lib
}

OPENSSL_PREFIX=$INSTALL_PREFIX/ssl
install_openssl() {
  if [[ ! -d $OPENSSL_PREFIX ]]; then
    echo "Installing OpenSSL"
    pushd /tmp
    wget https://www.openssl.org/source/openssl-$OPENSSL_VERSION.tar.gz &&
      tar -xvzf openssl-$OPENSSL_VERSION.tar.gz &&
      cd openssl-$OPENSSL_VERSION &&
      if [[ `uname` == "Linux" ]]; then CFLAGS=-fPIC ./config -fPIC -shared --prefix=$OPENSSL_PREFIX; else    ./Configure darwin64-x86_64-cc shared -fPIC --prefix=$OPENSSL_PREFIX; fi &&
      make && make install && echo "Installing OpenSSL DONE"
    rm -fr /tmp/openssl*
    popd
  fi
  add_to_env OPENSSL_ROOT_DIR $OPENSSL_PREFIX
  add_to_env LD_LIBRARY_PATH $OPENSSL_PREFIX/lib
}

CURL_PREFIX=$INSTALL_PREFIX
install_curl() {
  if [[ `uname` == "Darwin" ]]; then
    # curl is supported natively in macOS
    return 0
  fi
  if [[ ! -f $CURL_PREFIX/libcurl.so ]]; then
    echo "Installing CURL"
    pushd /tmp
    git clone https://github.com/curl/curl.git &&
      cd curl &&
      autoreconf -i &&
      ./configure --enable-lib-only --with-pic  --with-ssl=$OPENSSL_PREFIX --prefix $CURL_PREFIX &&
      make && make install && echo "Installing CURL DONE"
    rm -fr /tmp/curl
    popd
  fi
  add_to_env LD_LIBRARY_PATH $CURL_PREFIX/lib
}

UUID_PREFIX=$INSTALL_PREFIX
install_uuid() {
  if [[ `uname` == "Darwin" ]]; then
    # libuuid comes via brew install ossp-uuid; it might even be supported natively in macOS now
    return 0
  fi
  if [[ ! -f $UUID_PREFIX/libuuid.a ]]; then
    echo "Installing libuuid"
    pushd /tmp
    wget https://sourceforge.net/projects/libuuid/files/libuuid-1.0.3.tar.gz &&
      tar -xvzf libuuid-1.0.3.tar.gz &&
      cd libuuid-1.0.3 &&
      sed -i s/2.69/2.63/ configure.ac &&
      aclocal &&
      automake --add-missing &&
      ./configure --with-pic CFLAGS="-I/usr/include/x86_64-linux-gnu" --disable-shared --enable-static --prefix $UUID_PREFIX &&
      autoreconf -i -f &&
      make && make install && echo "Installing libuuid DONE"
    rm -fr /tmp/libuuid*
    popd
  fi
  add_to_env LD_LIBRARY_PATH $UUID_PREFIX/lib
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
        source system/install_ubuntu_prereqs.sh
      else
        source system/install_centos_prereqs.sh
      fi
      install_system_prerequisites
      ;;
    Darwin )
      system/install_macos_prereqs.sh
      ;;
    * )
      echo "OS=`uname` not supported"
      exit 1
  esac
}

install_prerequisites() {
  echo "1 PREREQS_ENV=$PREREQS_ENV"
  install_os_prerequisites &&
    source $PREREQS_ENV &&
    install_mvn &&
    install_protobuf
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
if [[ $BUILD_DISTRIBUTABLE_LIBRARY == false && $CENTOS_VERSION -eq 6 ]]; then
  echo "Centos 6 is supported only when build-arg distributable_jar=true"
  exit 1
fi

RC=1
install_prerequisites $2 &&
  if [[ $BUILD_DISTRIBUTABLE_LIBRARY == true ]]; then
    echo "Installing static libraries"
    install_openssl &&
      install_curl &&
      install_uuid
  fi && RC=0
finalize $RC
