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

# Install GenomicsDB Dependencies
# Workaround for FindMPI bug in latest released cmake, this is apparently fixed in cmake 3.20
# See issue filed with cmake - https://gitlab.kitware.com/cmake/cmake/-/issues/21955
# brew list cmake &>/dev/null || brew install cmake
# Start workaround ----
brew list cmake &>/dev/null && brew uninstall cmake
wget https://github.com/Homebrew/homebrew-core/raw/f3a486150931102cf4ce7b5943b6d209c34c05fe/Formula/cmake.rb
brew install --HEAD -s ./cmake.rb
# End workaround ------

brew list mpich &>/dev/null || brew install mpich
brew list ossp-uuid &>/dev/null || brew install ossp-uuid
brew list libcsdb &>/dev/null || brew install libcsv
brew list automake &> /dev/null || brew install automake
brew list openssl &> /dev/null || brew install openssl
