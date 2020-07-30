#!/bin/bash

set -e 

# Install GenomicsDB Dependencies
brew list cmake &>/dev/null || brew install cmake
brew list mpich &>/dev/null || brew install mpich
brew list ossp-uuid &>/dev/null || brew install ossp-uuid
brew list libcsdb &>/dev/null || brew install libcsv
brew list automake &> /dev/null || brew install automake
