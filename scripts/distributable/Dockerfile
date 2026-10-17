# The MIT License (MIT)
# Copyright (c) 2026 Tim Fennell
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
# Toolchain for GenomicsDB's distributable Linux library (scripts/distributable/build.sh): GCC 15, clang 21 and
# CMake 4 on AlmaLinux 8. The library needs a glibc at least as new as the image's, 2.28, so the base image sets the
# oldest Linux distributions the library loads on.
FROM almalinux:8
RUN dnf -y install gcc-toolset-15-gcc gcc-toolset-15-gcc-c++ gcc-toolset-15-binutils llvm-toolset \
      git make which patch diffutils file zlib-devel python3 java-17-openjdk-devel \
      autoconf automake libtool perl-IPC-Cmd perl-Time-Piece \
      && dnf -y module install maven:3.8 && dnf clean all
ARG CMAKE_VERSION=4.4.4
RUN curl -sfL https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-$(uname -m).tar.gz \
      | tar xz -C /opt && ln -s /opt/cmake-${CMAKE_VERSION}-linux-$(uname -m)/bin/* /usr/local/bin/
# The source tree is mounted from the host, owned by someone else, and the build runs git in it
RUN git config --system safe.directory '*'
# clang uses GCC 15's libstdc++ rather than the system GCC 8's
RUN echo "--gcc-install-dir=/opt/rh/gcc-toolset-15/root/usr/lib/gcc/$(uname -m)-redhat-linux/15" \
      | tee /usr/bin/clang.cfg > /usr/bin/clang++.cfg
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk
