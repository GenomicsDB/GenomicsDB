#!/bin/bash

install_protobuf() {
  if [ ! -d $PROTOBUF_HOME ]; then
    mkdir -p $PROTOBUF_HOME
  fi
  cd $TRAVIS_BUILD_DIR/dependencies
  if [ ! -f $CACHE_DIR/protobuf.tar ]; then
    echo "Installing Protobuf"
    PROTOBUF_DIR=$TRAVIS_BUILD_DIR/dependencies/protobuf-$PROTOBUF_VERSION
    wget -nv https://github.com/protocolbuffers/protobuf/releases/download/v$PROTOBUF_VERSION/protobuf-cpp-$PROTOBUF_VERSION.zip &&
    unzip protobuf-cpp-$PROTOBUF_VERSION.zip &&
    cd protobuf-$PROTOBUF_VERSION &&
    cp $TRAVIS_BUILD_DIR/.travis/scripts/protobuf-v$PROTOBUF_VERSION.autogen.sh.patch autogen.sh &&
    ./autogen.sh &&
    ./configure --prefix=$PROTOBUF_HOME --with-pic &&
    make -j4 && make install &&
    cd $PROTOBUF_HOME/
    tar -cf  $CACHE_DIR/protobuf.tar bin lib include
  else
	echo "Extracting protobuf"
    cd $PROTOBUF_HOME
    tar -xf $CACHE_DIR/protobuf.tar
  fi
}

echo "Installing OSX Dependencies"
install_protobuf &&
echo "Installing OSX Dependencies DONE"
