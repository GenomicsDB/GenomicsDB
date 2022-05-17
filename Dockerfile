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
# Description: Docker file for building and installing GenomicsDB

# OS currently tested are ubuntu and centos
#ARG os=ubuntu:trusty
ARG os=centos:7
FROM $os as full

ARG user=genomicsdb
ARG install_dir=/usr/local
ARG distributable_jar=false
# Options to enable_bindings are only java for now
# e.g. enable_bindings="java"
ARG enable_bindings=""

COPY . /build/GenomicsDB/
ENV DOCKER_BUILD=true
WORKDIR /build/GenomicsDB
RUN ./scripts/prereqs/install_prereqs.sh ${distributable_jar} full ${enable_bindings} &&\
    useradd -r -U -m ${user} &&\
    ./scripts/install_genomicsdb.sh ${user} ${install_dir} ${distributable_jar} ${enable_bindings}

ARG os=centos:7

FROM $os as release

ARG user=genomicsdb
ARG install_dir=/usr/local
ARG distributable_jar=false
ARG enable_bindings=""
COPY ./scripts/prereqs /build
WORKDIR /build
RUN ./install_prereqs.sh ${distributable_jar} release ${enable_bindings} &&\
    useradd -r -U -m ${user}

COPY --from=full /usr/local/bin/*genomicsdb* /usr/local/bin/gt_mpi_gather /usr/local/bin/vcf* ${install_dir}/bin/

USER ${user}
WORKDIR /home/${user}
ENTRYPOINT ["/bin/bash", "--login"]
