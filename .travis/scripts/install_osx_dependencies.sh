#!/bin/bash

install_gtest() {
  GTEST_DIR=$TRAVIS_BUILD_DIR/dependencies/googletest-release-1.8.1
  if [ ! -f $CACHE_DIR/googletest.tar ]; then
     echo "Installing google test"
     cd $TRAVIS_BUILD_DIR/dependencies
     wget -nv https://github.com/google/googletest/archive/release-1.8.1.tar.gz &&
     tar -xz release-1.8.1.tar.gz &&
     cd googletest-release-1.8.1/googletest && cmake . -DBUILD_SHARED_LIBS=1 && make &&
     cd $GTEST_DIR && tar cf $CACHE_DIR/googletest.tar googletest/libgtest* googletest/include/gtest
  else
     echo "Extracting gtest"
     if [ ! -d $GTEST_DIR ]; then
       mkdir -p $GTEST_DIR
     fi
     cd $GTEST_DIR
     tar -xf $CACHE_DIR/googletest.tar
  fi
  cd $GTEST_DIR/googletest &&
  cp libgtest* /usr/local/lib && 
  cp -r include/gtest /usr/local/include &&
  echo "Installing google test SUCCESSFUL"
}

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
install_gtest &&
install_protobuf &&
echo "Installing OSX Dependencies DONE"
