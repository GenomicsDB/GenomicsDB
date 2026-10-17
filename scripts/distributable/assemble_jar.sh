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
# Assembles the multi-platform release from the Linux x86-64 build's jars, sources, javadoc and pom plus the Linux
# aarch64 and macOS arm64 native libraries: genomicsdb-<version>.jar and genomicsdb-<version>-allinone-spark.jar with
# all three libraries, genomicsdb-<version>-sources.jar, -javadoc.jar and genomicsdb-<version>.pom.
#
# usage: assemble_jar.sh <version> <x86-64 target dir> <x86-64 pom.xml> <aarch64 libtiledbgenomicsdb.so>
#                        <macOS libtiledbgenomicsdb.dylib> <output dir>
set -euo pipefail
VERSION=$1
TARGET=$2
POM=$3
AARCH64_LIB=$4
MACOS_LIB=$5
OUT=$6
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

mkdir -p "$OUT" "$STAGE/linux-aarch64" "$STAGE/macos-aarch64"
cp "$AARCH64_LIB" "$STAGE/linux-aarch64/libtiledbgenomicsdb.so"
cp "$MACOS_LIB" "$STAGE/macos-aarch64/libtiledbgenomicsdb.dylib"
OUT=$(cd "$OUT" && pwd)
for jar in "genomicsdb-$VERSION.jar" "genomicsdb-$VERSION-allinone-spark.jar"; do
  cp "$TARGET/$jar" "$OUT/$jar"
  (cd "$STAGE" && zip -q "$OUT/$jar" linux-aarch64/libtiledbgenomicsdb.so macos-aarch64/libtiledbgenomicsdb.dylib)
  libraries=$(unzip -l "$OUT/$jar" | grep -cE '(linux-x86_64|linux-aarch64|macos-aarch64)/libtiledbgenomicsdb')
  if [ "$libraries" != 3 ]; then
    echo "$jar has $libraries native libraries instead of 3" >&2
    exit 1
  fi
done
cp "$TARGET/genomicsdb-$VERSION-sources.jar" "$TARGET/genomicsdb-$VERSION-javadoc.jar" "$OUT/"
cp "$POM" "$OUT/genomicsdb-$VERSION.pom"
ls -l "$OUT"
