#!/bin/bash

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
	yum install -y java-1.8.0-openjdk-devel &&
	echo "export JRE_HOME=/usr/lib/jvm/jre" >> $PREREQS_ENV
}

install_mpi() {
	yum install -y mpich-devel &&
	echo "export PATH=/usr/lib64/mpich/bin:\$PATH" >> $PREREQS_ENV &&
	echo "export LD_LIBRARY_PATH=/usr/lib64:/usr/lib:\$LD_LIBRARY_PATH" >> $PREREQS_ENV
}

install_curl() {
  yum install -y curl libcurl-devel
  return 0
}

install_csv() {
  if [[ $CENTOS_VERSION -lt 8 ]]; then 
    yum install -y libcsv libcsv-devel &&
      yum install -y csv
  fi
}

install_test_prerequisites() {
  if [[ $CENTOS_VERSION -ge 8 ]]; then
    yum install -y python3 &&
      pip3 install jsondiff
  else
	  yum install -y python-pip &&
	    pip install virtualenv &&
	    pip install jsondiff &&
	    yum install -y lcov
  fi
}

install_prerequisites_centos() {
	yum update -y -q &&
	  yum install -y sudo &&
	  yum install -y -q which wget git make &&
	  install_mpi &&
	  install_devtoolset &&
	  install_openjdk &&
	  yum install -y autoconf automake libtool unzip &&
	  yum update -y autoconf &&
	  yum install -y epel-release &&
	  yum install -y zlib-devel &&
	  yum install -y openssl-devel &&
    yum install -y cmake3 &&
	  yum install -y libuuid libuuid-devel &&
    install_csv &&
    install_curl &&
    install_test_prerequisites
}
