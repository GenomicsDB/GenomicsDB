
.. _API Docs: 

###############################
APIs
###############################


Overview
*******************************
In addition to the native C++ interface, GenomicsDB has support/bindings for Python, Java, and R.
The data structure returned from a GenomicsDB query depends on the API or CLI tool being used and the method invoked. 
Similarly, the method of passing in query field information depends on the route chosen.

For example, the CLI tool gt_mpi_gather takes in a query in json format. 

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


The various APIs take in the same query fields in one form or another, so check the approptiate reference for your use case. 
If you are new to GenomicsDB, we recommend starting with the :ref:`gt_mpi_gather <CLI Tools gt-mpi-gather>` CLI tool to familiarize yourself with the concepts. 

All the APIs support taking configuration options as protobuf. Here's the protobuf documentation:

.. toctree::
   :maxdepth: 2

   protodoc/proto.md


C++
*******************************
* Native interface - See `cpp source code`_ (Reference coming soon)
* Querying CLI tools available - see :ref:`vcf2genomicsdb <CLI Tools vcf2genomicsdb>`

.. _cpp source code: https://github.com/GenomicsDB/GenomicsDB/tree/master/src/main/cpp

API Reference
===============================

.. toctree::
   :maxdepth: 3

   doxyrst/cplusplus

Python
*******************************
* Python bindings for GenomicsDB. Requires Python version *>= 3.6*
* Built using Cython.
* Currently only supports querying (no importing).
* Data can be returned in either a numpy matrix or a list of lists.

The GenomicsDB python bindings are `published on PyPI`_ (currently experimental).

.. _published on PyPI: https://pypi.org/project/genomicsdb/

Install ``genomicsdb`` with the following command:

.. code-block:: python

   pip install genomicsdb


API Reference
===============================
* Coming Soon. See the `bindings source code`_ for the time being. 

.. _bindings source code: https://github.com/GenomicsDB/GenomicsDB-Python

Java
*******************************
* Data can be returned in either a htsjdk `VariantContext`_ object or <other>.
* The GenomicsDB java bindings are published to `Maven Central`_ every release.

.. _VariantContext: https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
.. _Maven Central: https://mvnrepository.com/artifact/org.genomicsdb/genomicsdb


* Add GenomicsDB to your maven project with the following snippet (update the version, if necessary):

.. code-block:: java

   <!-- https://mvnrepository.com/artifact/org.genomicsdb/genomicsdb -->
   <dependency>
      <groupId>org.genomicsdb</groupId>
      <artifactId>genomicsdb</artifactId>
      <version>1.4.5</version>
   </dependency>


API Reference
===============================

.. toctree::
   :maxdepth: 3

   doxyrst/javaimporter
   doxyrst/javareader





