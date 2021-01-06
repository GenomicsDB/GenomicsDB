.. _CLI Tools:

###############################
CLI Tools
###############################
GenomicsDB contains several native command line tools. The artifacts are not published online. 
Users need to :ref:`build and install <Building + Installing>` GenomicsDB.

.. _CLI Tools vcf2genomicsdb:

* *vcf2genomicsdb*: 
  
  Tool for importing VCF files into GenomicsDB. 
  
  Example usage:

  .. code-block:: bash

    vcf2genomicsdb <loader-json.json>
.. ------------------------------

.. _CLI Tools gt-mpi-gather: 

* *gt_mpi_gather*: 
  
  Tool for querying GenomicsDB. Outputs results as either *VariantCalls* or *Variants*
  
  Example usage:

  .. code-block:: bash

    gt_mpi_gather <loader-json.json>

.. ------------------------------

* create_genomicsdb_workspace
* consolidate_genomicsdb_array
* histogram
* vcfdiff

