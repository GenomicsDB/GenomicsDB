#!/bin/bash

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
# Builds GenomicsDB's distributable native library and jars the way releases are built, runs the C++ tests on the
# library, checks what it depends on (tests/check_distributable_library.sh) and loads it from the jar. The build goes
# in build/ and the dependencies built from source in deps/.
#
# On Linux, run it from the repository root in a container of the image this directory's Dockerfile builds, with CC
# and CXX naming the compiler (GCC by default):
#   docker build -t genomicsdb-toolchain scripts/distributable
#   docker run --rm -v "$PWD":/src -w /src genomicsdb-toolchain scripts/distributable/build.sh <version> [options]
# On macOS, set MACOSX_DEPLOYMENT_TARGET to the oldest macOS the library should load on, e.g. 14.0.
#
# usage: scripts/distributable/build.sh <version> [cmake options...]

set -euo pipefail

VERSION=$1
shift
if [ "$(uname)" = Darwin ]; then
  : "${MACOSX_DEPLOYMENT_TARGET:?must name the oldest macOS the library should load on, e.g. 14.0}"
  JOBS=$(sysctl -n hw.ncpu)
else
  source /opt/rh/gcc-toolset-15/enable
  export CC=${CC:-gcc} CXX=${CXX:-g++}
  JOBS=$(nproc)
fi
export MAVEN_ARGS="-Dgenomicsdb.version=$VERSION"

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_DISTRIBUTABLE_LIBRARY=1 -DBUILD_JAVA=1 \
  -DGENOMICSDB_RELEASE_VERSION="$VERSION" "$@"
cmake --build build --target "genomicsdb-$VERSION-examples" ctests api_tests -j "$JOBS" 2>&1 | tee build/build.log
# The macOS linker only warns when something linked in targets a newer macOS than the library
if grep 'built for newer' build/build.log; then
  exit 1
fi
(cd build && ctest --output-on-failure -R '^(distributable_library|ctests|api_tests)$')
java -cp "build/target/genomicsdb-$VERSION-allinone-spark.jar" scripts/distributable/SmokeLoadNativeLib.java
