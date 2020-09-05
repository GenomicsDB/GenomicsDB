# GenomicsDB - Build and Install
Experimental scripts to build and install GenomicsDB libraries and tools.

## JVM
GenomicsDB jars bundled with pre-built native libraries are regularly published to [![Maven Central](https://img.shields.io/maven-central/v/org.genomicsdb/genomicsdb.svg)](https://mvnrepository.com/artifact/org.genomicsdb). These artifacts can be included as dependencies appropriate to the build(maven/gradle/sbt) system in use. See [Examples](../example/java/) on java/spark usage.

## Bash
Scripts [install_prereqs.sh](prereqs/install_prereqs.sh) is meant to work with Docker, but should work on any Linux distribution from a bash shell. Most of the code is self-explanatory, feel free to modify the scripts as desired.  Note that [install_prereqs.sh](prereqs/install_prereqs.sh) requires sudoable access.

```bash
# Install the prerequisites to build GenomicsDB and create a genomicsdb_prereqs.sh in $HOME
sudo prereqs/install_prereqs.sh
source $HOME/genomicsdb_prereqs.sh
cmake -DCMAKE_PREFIX_INSTALL=$HOME /path/to/GenomicsDB # simplest
# Build and install the include/lib/bin files to $HOME
make && make install
```

## Docker
To build and install GenomicsDB using Docker, specify the following *optional* build arguments

  | Build Argument | Default |
  | --- | --- |
  | os=ubuntu:trusty\|centos:7\|\<linux_base:ver\> | centos:7 |
  | user=<user_name> | genomicdb |
  | branch=master\|develop\|<any_branch> | master |
  | install_dir=<my_install_dir> | /usr/local |
  | distributable_jar=true\|false | false |  
  | enable_bindings=java | none |
  
Examples:
```
docker build --build-arg os=ubuntu:trusty --build-arg branch=develop --build-arg install_dir=/home/$USER -t genomicsdb:build . 
```

To run and enter the bash shell:
```
# Use the -t argument value used with docker build ...
docker run -it genomicsdb:build
```

To build and copy all built artifacts from the docker image:
```
export docker_os=centos
export docker_repo=genomicsdb
export docker_tag=`date "+%Y-%m-%d-%H:%M:%S"`
docker build --build-arg os=$docker_os --build-arg install_dir=/tmp/artifacts -t $docker_repo:$docker_tag .
docker create -it --name genomicsdb $docker_repo:$docker_tag bash
docker cp genomicsdb:/tmp/artifacts $docker_os
docker rm -fv genomicsdb
```

