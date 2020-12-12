#!/bin/bash

## Build script for genomicsdbproto python package

SRC=../src/resources/
BUILD=./proto/genomicsdbproto

mkdir -p proto/genomicsdbproto

## Step1: Generate python code from proto
protoc --proto_path=$SRC --python_out=$BUILD $SRC/*.proto

cd proto

## Step2: Setup package

cp ../../LICENSE .

cat << EOF > README.md
# GenomicsDB Protocol buffers for Python

This offers protocol buffers for GenomicsDB [GenomicsDB](https://github.com/GenomicsDB/GenomicsDB/)
EOF

touch genomicsdbproto/__init__.py

cat <<EOF > setup.py 
import setuptools
from distutils.core import setup

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="genomicsdbproto",
  version="0.0.1",
  author="Melvin Lathara",
  author_email="melvin@omicsautomation.com",
  description="Protocol buffers for GenomicsDB",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/GenomicsDB/GenomicsDB",
  packages=setuptools.find_packages(),
  install_requires=['protobuf','wheel'],
  python_requires='>=3.0',
  classifiers=[
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
  ],
)
EOF

## Step3: Create the wheel
python3 -m venv venv
source venv/bin/activate

pip3 install --upgrade pip setuptools wheel

python3 setup.py sdist bdist_wheel
deactivate
