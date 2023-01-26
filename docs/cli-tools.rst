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

  Here is an example `query.json` file:

  .. code-block:: python

    {
        "workspace" : "/tmp/ws/",
        "array" : "t0_1_2",
        "query_column_ranges" : [ [ [0, 100 ], 500 ] ],
        "query_row_ranges" : [ [ [0, 2 ] ] ],
        "vid_mapping_file": "tests/inputs/vid.json",
        "callset_mapping_file": "tests/inputs/callset_mapping.json",
        "query_attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP" ]
    }

  Most of the fields are self-explanatory (take a look at the :ref:`terminology section <Concepts + Terminology>` section).
  The following fields are mandatory:
  
  * *workspace* (type: string or list of strings)
  
  * *array* (type: string or list of strings)
  
  * *query_column_ranges* : This field contains is a list of lists. Each member list contains the column ranges that are to be queried. Each element of the inner list can be a single integer, representing the single column position to be queried, or a list of size 2 representing the column range to be queried.
  
  In the above example, the process will query column range \[0-100\] (inclusive) and the position 500 and return all *VariantCalls* intersecting with these query ranges/positions.
  
  For a more detailed explanation as to why this field is a list of lists (and other fields are lists of strings), we 
  refer the reader to the [[wiki page explaining how we use MPI in the context of GenomicsDB|MPI-with-GenomicsDB]].
  
  * *query_attributes* : List of strings specifying attributes to be fetched
  
  Optional field(s):
  
  * *query_row_ranges* : Similar to *query_column_ranges* but for rows. If this field is omitted, all rows are included in the query.
  
  * *vid_mapping_file* and *callset_mapping_file* (type: string or list of strings): Paths to JSON files specifying the :ref:`vid mapping <Import / ETL>` and :ref:`callset mapping <Import / ETL>`.
  
    These two fields are optional because a user might specify them in the loader JSON while creating an array and pass the loader JSON to the query tool(s) (see below). This allows the user to write many different query JSON files without repeating the information.
  
    If the *vid_mapping_file* and/or *callset_mapping_file* are specified in both the loader and query JSON files and passed to the query tool(s), then the parameter value in the query JSON gets precedence.

  Query result format
  ^^^^^^^^^^^^^^^^^^^

  `gt_mpi_gather` returns data in a JSON format by default. Query results may be returned as VariantCalls or Variants depending on command line parameters passed by the user. Other formats (gVCF, PLINK) are also supported but these are well described elsewhere and so omitted here. 

  A VariantCall object contains information stored for a given sample/CallSet for a given location/genomic interval in GenomicsDB, i.e, it contains all the fields stored in the corresponding GenomicsDB array cell. A sample VariantCall formatted as a JSON is shown below:
  
  .. code-block:: python

    {
        "row": 2,
        "interval": [ 17384, 17386 ],
        "fields": {
            "REF": "GAT",
            "ALT": [ "GT","<NON_REF>" ],
            "PL": [ 1018,0,1116,1137,1224,2361 ],
            "BaseQRankSum": [ 1.046 ]
        }
    }

  The above example shows a VariantCall object for the sample/CallSet corresponding to row 2 in the GenomicsDB at column 17384. Since the VariantCall is a deletion, it spans multiple columns and ends at column 17386 (inclusive). The fields field is self explanatory.
  
  A Variant is a set of VariantCalls which meet the following properties:
  - All VariantCalls span the exact same column interval
  - All VariantCalls have the same value of REF
  - All VariantCalls have the same value of ALT
  A sample Variant formatted as a JSON is shown below:
  
  .. code-block:: python

    {
        "interval": [ 17384, 17384 ],
            "common_fields" : {
                "REF": "G",
                "ALT": [ "A","<NON_REF>" ]
            },
            "variant_calls": [
            {
                "row": 0,
                "interval": [ 17384, 17384 ],
                "fields": {
                    "REF": "G",
                    "ALT": [ "A","<NON_REF>" ],
                    "BaseQRankSum": [ -2.096000 ]
                }
            },
            {
                "row": 2,
                "interval": [ 17384, 17384 ],
                "fields": {
                    "REF": "G",
                    "ALT": [ "A","<NON_REF>" ],
                    "BaseQRankSum": [ 1.046000 ]
                }
            }
        ]
    }

  In the above example, the Variant object corresponds to interval [17384, 17384] with "G" as the reference allele and "A", "<NON_REF>" as the alternate alleles. It consists of two VariantCalls at rows 0 and 2.
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