.. _CLI Tools:

###############################
CLI Tools
###############################
GenomicsDB contains several native command line tools. These tools are also available through GenomicsDB docker. 
Refer to :ref:`this link <Building + Installing>` for instructions on building and installing GenomicsDB, as well as how to work with artifacts in docker images.

.. _CLI Tools vcf2genomicsdb_init:

vcf2genomicsdb_init 
-------------------
  
  Tool for creating a GenomicsDB workspace, and initializing it with configuration json files. 
  
  Example usage:

  .. code-block:: bash

    vcf2genomicsdb_init -w /path/to/workspace -s <list-of-sample-vcfs>

.. ------------------------------

.. _CLI Tools vcf2genomicsdb:

vcf2genomicsdb
--------------
  
  Tool for importing VCF files into GenomicsDB. 
  
  Example usage:

  .. code-block:: bash

    vcf2genomicsdb <loader-json.json>

.. ------------------------------

.. _CLI Tools gt-mpi-gather: 

gt_mpi_gather
-------------
  
  Tool for querying GenomicsDB. Outputs results as either *VariantCalls* or *Variants*
  
  Example usage:

  .. code-block:: bash

    gt_mpi_gather -j <query.json>

.. ------------------------------

.. _CLI Tools create_genomicsdb_workspace:

create_genomicsdb_workspace
---------------------------

  Tool for creating a new GenomicsDB workspace.

  Example usage:

  .. code-block:: bash

    create_genomicsdb_workspace /path/to/workspace

.. ------------------------------

.. _CLI Tools consolidate_genomicsdb_array:

consolidate_genomicsdb_array
----------------------------

  Consolidate fragments from incremental imports into a single fragment.

  Example usage:

  .. code-block:: bash

    consolidate_genomicsdb_array -w /path/to/workspace -a <name-of-array-to-consolidate>

.. ------------------------------

.. _CLI Tools vcf_histogram:

vcf_histogram
-------------

  Create a histogram showing the number of variants per histogram bin. Can be useful in deciding how to partition the GenomicsDB workspace.

  Example usage:

  .. code-block:: bash

    vcf_histogram <json>

.. ------------------------------

.. _CLI Tools vcf_diff:

vcfdiff
-------

  Check whether VCFs are identical

  Example usage:

  .. code-block:: bash

    vcfdiff /path/to/vcf1 /path/to/vcf2