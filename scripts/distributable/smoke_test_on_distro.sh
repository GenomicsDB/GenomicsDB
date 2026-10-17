#!/bin/sh

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
# Runs in a plain Linux distribution's image, with a JDK from the host mounted, and loads the jar's native library for
# this machine in a JVM, which fails if a library it needs is missing.
#
# usage: smoke_test_on_distro.sh <JDK directory> <genomicsdb-<version>-allinone-spark.jar>
set -e
JDK=$1
JAR=$2
. /etc/os-release
echo "$PRETTY_NAME, $(getconf GNU_LIBC_VERSION)"
"$JDK/bin/java" -cp "$JAR" "$(dirname "$0")/SmokeLoadNativeLib.java"
