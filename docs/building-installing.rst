.. _Building + Installing: 

###############################
Building + Installing
###############################

This section refers to building and installing the core GenomicsDB executables and tools. For including published artifacts in your code, see the appropriate :ref:`API docs <API Docs>` and :ref:`Examples <Examples>` sections. 

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
Docker images for GenomicsDB are published to GitHub Container Registry. `This page`_ shows the image versions available. Below we show some examples of using the docker image to interact with GenomicsDB executables.

.. _This page: https://github.com/GenomicsDB/GenomicsDB/pkgs/container/genomicsdb/versions

 
Examples:

To ingest vcfs from `/path/to/vcfs` to `/path/to/workspace` (both paths refer to paths on the host machine)

.. code-block:: bash

    host> docker run -v /path/to/vcfs:/opt/vcfs -v /path/to/where/workspace/should/be/created:/opt/workspace_parent -u $( id -u $USER ):$( id -g $USER ) -it ghcr.io/genomicsdb/genomicsdb:v1.4.4
    docker> vcf2genomicsdb_init -w /opt/workspace_parent/workspace -S /path/to/vcfs -n 0
    docker> vcf2genomicsdb /opt/workspace_parent/loader.json


To query a GenomicsDB workspace, and return the results as a VCF/gVCF:

.. code-block:: bash

    host> docker run -v /path/to/workspace:/opt/workspace -u $( id -u $USER ):$( id -g $USER ) -it ghcr.io/genomicsdb/genomicsdb:v1.4.4
    # run with no arguments for help
    docker> gt_mpi_gather
    # Create a query.json file based on help
    docker> gt_mpi_gather -j /path/to/query.json --produce-Broad-GVCF


CLI Tools
*******************************
The CLI tools are compiled during the build process and placed in <build-dir>/tools.