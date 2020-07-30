#!/bin/bash

set -e

install_openjdk8() {
	apt-get -y install openjdk-8-jdk icedtea-plugin &&
		apt-get update -q &&
		pushd $HOME &&
		git clone https://github.com/michaelklishin/jdk_switcher.git &&
		source jdk_switcher/jdk_switcher.sh && jdk_switcher use openjdk8 &&
		echo "export JAVA_HOME=$JAVA_HOME" >> $PREREQS_ENV &&
		popd
}

install_cmake3() {
	wget -nv https://github.com/Kitware/CMake/releases/download/v3.4.0/cmake-3.4.0-Linux-x86_64.sh -P /tmp &&
		chmod +x /tmp/cmake-3.4.0-Linux-x86_64.sh &&
		/tmp/cmake-3.4.0-Linux-x86_64.sh --prefix=/usr/local --skip-license
}

install_R() {
  apt-get -y install gnupg &&
  apt-get install -y software-properties-common &&
  apt-get update -q &&
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &&
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' &&
  apt-get -y install libxml2-dev &&
	apt-get update -q &&
	apt-get -y install r-base
}

install_prerequisites_ubuntu() {
	apt-get update -q &&
    apt-get install -y tzdata &&
		apt-get -y install \
						lcov \
						mpich \
						zlib1g-dev \
						libssl-dev \
						rsync \
						libidn11 \
						uuid-dev \
						libcurl4-openssl-dev \
						build-essential \
						software-properties-common \
						wget \
						git \
						autoconf \
						automake \
						libtool \
						zip \
						unzip \
						curl \
						sudo &&
		apt-get update -q &&
		apt-get update -q &&
		install_openjdk8 &&
		install_R &&
		apt-get clean &&
		apt-get purge -y &&
		install_cmake3 &&
		rm -rf /var/lib/apt/lists*
}
