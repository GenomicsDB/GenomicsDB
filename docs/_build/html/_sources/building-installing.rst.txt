.. _Building + Installing: 

###############################
Building + Installing
###############################

This section refers to building and installing the core GenomicsDB executables and tools. For including published artifacts in your code, see the appropriate :ref:`API docs <API Docs>` and *Examples* sections. 
For more details on dependencies and manually compiling GenomicsDB, see the `old docs`_ for now. 

.. _old docs: https://github.com/GenomicsDB/GenomicsDB/wiki/Compiling-GenomicsDB#building

Bash
*******************************
The script *install_prereqs.sh* is meant to work with Docker, but should work on any Linux distribution from a bash shell. Note that *install_prereqs.sh* requires sudo-able access.

.. code-block:: bash

    # Clone GenomicsDB and change to scripts directory
    git clone https://github.com/GenomicsDB/GenomicsDB && cd GenomicsDB/scripts/

    # Install the prerequisites to build GenomicsDB and create a genomicsdb_prereqs.sh in $HOME
    sudo prereqs/install_prereqs.sh
    source $HOME/genomicsdb_prereqs.sh
    cmake -DCMAKE_PREFIX_INSTALL=$HOME /path/to/GenomicsDB # simplest

    # Build and install the include/lib/bin files to $HOME
    make && make install
    

Docker
*******************************
To build and install GenomicsDB using Docker, specify the following *optional* build arguments

.. list-table::
   :widths: 50 25 
   :header-rows: 1

   * - Build Argument
     - Default
   * - os=ubuntu:trusty|centos:7|<linux_base:ver>
     - centos:7
   * - user=<user_name>
     - genomicsdb
   * - branch=master|develop|<any_branch>
     - master
   * - install_dir=<my_install_dir>
     - /usr/local
   * - distributable_jar=true|false
     - false
   * - enable_bindings=java
     - none

  
Examples:

.. code-block:: bash

    docker build --build-arg os=ubuntu:trusty --build-arg branch=develop --build-arg install_dir=/home/$USER -t genomicsdb:build . 


To run and enter the bash shell:

.. code-block:: bash

    # Use the -t argument value used with docker build ...
    docker run -it genomicsdb:build


To build and copy all built artifacts from the docker image:

.. code-block:: bash

    export docker_os=centos
    export docker_repo=genomicsdb
    export docker_tag=`date "+%Y-%m-%d-%H:%M:%S"`
    docker build --build-arg os=$docker_os --build-arg install_dir=/tmp/artifacts -t $docker_repo:$docker_tag .
    docker create -it --name genomicsdb $docker_repo:$docker_tag bash
    docker cp genomicsdb:/tmp/artifacts $docker_os
    docker rm -fv genomicsdb



CLI Tools
*******************************
The CLI tools are compiled using *cmake* during the build process and placed in <build-dir>/tools.