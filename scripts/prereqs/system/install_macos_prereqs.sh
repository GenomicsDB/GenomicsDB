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

PREREQS_ENV=${PREREQS_ENV:-$HOME/genomicsdb_prereqs.sh}

# Install GenomicsDB Dependencies
HOMEBREW_NO_AUTO_UPDATE=1
HOMEBREW_NO_INSTALL_CLEANUP=1

build_openssl() {
  pushd openssl-$OPENSSL_VERSION &&
    ./Configure darwin64-$ARCH-cc no-tests -fPIC --prefix=$OPENSSL_PREFIX/$ARCH &&
    make && make install_sw && make clean &&
    popd
}

install_universal_openssl() {
  mkdir -p $OPENSSL_PREFIX/lib &&
    mkdir -p $OPENSSL_PREFIX/include &&
    lipo -create $OPENSSL_PREFIX/x86_64/lib/libcrypto.dylib $OPENSSL_PREFIX/arm64/lib/libcrypto.dylib -output $OPENSSL_PREFIX/lib/libcrypto.dylib &&
    lipo -create $OPENSSL_PREFIX/x86_64/lib/libcrypto.a $OPENSSL_PREFIX/arm64/lib/libcrypto.a -output $OPENSSL_PREFIX/lib/libcrypto.a &&
    lipo -create $OPENSSL_PREFIX/x86_64/lib/libssl.dylib $OPENSSL_PREFIX/arm64/lib/libssl.dylib -output $OPENSSL_PREFIX/lib/libssl.dylib &&
    lipo -create $OPENSSL_PREFIX/x86_64/lib/libssl.a $OPENSSL_PREFIX/arm64/lib/libssl.a -output $OPENSSL_PREFIX/lib/libssl.a &&
    mv $OPENSSL_PREFIX/x86_64/include/openssl $OPENSSL_PREFIX/include
}

install_openssl() {
  if [[ ! -f $OPENSSL_PREFIX/lib/libcrypto.a ]]; then
    echo "Installing OpenSSL"
    pushd /tmp
    wget -nv $WGET_NO_CERTIFICATE https://www.openssl.org/source/openssl-$OPENSSL_VERSION.tar.gz &&
      tar -xzf openssl-$OPENSSL_VERSION.tar.gz &&
      ARCH=x86_64 build_openssl &&
      ARCH=arm64 build_openssl &&
      install_universal_openssl &&
      echo "Installing Universal OpenSSL DONE"
    rm -fr /tmp/openssl*
    popd
  fi
  add_to_env OPENSSL_ROOT_DIR $OPENSSL_PREFIX
}

# Sufficient to build/install the genomicsdb libraries
install_minimum_prerequisites() {
  echo "Installing minimum prerequisites"
  brew list cmake &>/dev/null || brew install cmake
  brew list automake &> /dev/null || brew install automake
  brew list pkg-config &> /dev/null || brew install pkg-config
  echo "Installing minimum prerequisites DONE"
}

install_system_prerequisites() {
  echo "Installing system prerequisites..."
  brew list mpich &>/dev/null || brew install mpich
  brew list openssl@3 &> /dev/null || brew install openssl@3
  add_to_env OPENSSL_ROOT_DIR $(brew --prefix openssl@3)

  # brew has started installing lcov 2.0 and some GenomicsDB sources are erroring out while running lcov
  # For example -
  # geninfo: ERROR: "/Users/runner/work/GenomicsDB/GenomicsDB/src/main/cpp/include/query_operations/variant_operations.h":50: function _ZN23RemappedDataWrapperBaseC2Ev end line 37 less than start line
  # The errors can be suppressed, but installing the older version 1.16 explicitly for now
  wget -nv https://github.com/Homebrew/homebrew-core/raw/e92d2ae54954ebf485b484d8522104700b144fee/Formula/lcov.rb
  brew list lcov &> /dev/null && brew uninstall lcov
  brew install -s lcov.rb

  #brew list openjdk@17 &> /dev/null || brew install openjdk@17
  echo "Installing system prerequisites DONE"
}
