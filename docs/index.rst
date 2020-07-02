.. GenomicsDB documentation master file, created by
   sphinx-quickstart on Wed Jun 10 21:36:23 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. ######################################
   Welcome to GenomicsDB's documentation!
   ######################################
  
    
.. .. toctree::
   :maxdepth: 3
   :caption: Contents:

.. 
   Indices and tables
   ==================
.. 
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
 


###############################
Overview (DOCUMENTATION WORK-IN-PROGRESS)
###############################

.. 
   .. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
      :target: https://opensource.org/licenses/MIT



GenomicsDB is a highly performant scalable data storage written in C++ for importing, querying and transforming genomic
variant data. 


.. _Concepts + Terminology:

Concepts + Terminology
*******************************
This section outlines basic GenomicsDB concepts and terminology used throughout this documentation.

* *Cell*: Cell of an array (for example: A[i][j] in a 2-D array)

* *Workspace*: A directory on the file system under which multiple GenomicsDB arrays can be stored. The underlying array stores some metadata files under the workspace directory. All GenomicsDB import/create programs assume that if the specified workspace directory exists, then it's a valid workspace directory. If the user points to an existing directory which doesn't contain the workspace metadata, the create/import programs will exit with an exception. If the directory doesn't exist, then the programs will initialize a new workspace correctly.

* *Array*: Name of a GenomicsDB array

   Given a workspace and array name, the framework will store its data in the directory *<workspace>/<array>*.

* *Column-major* and *row-major* ordering: Denotes the order in which cells are stored on disk.

   * Column major storage implies that cells belonging to the same column are stored contiguously on disk (cells belonging to the same column but different rows are sorted by row id). When cells are stored in column major order queries such as "retrieve all cells belonging to genomic position *X*" are fast due to high spatial locality. However, queries such as "retrieve cells for sample *Y*" are relatively slow. For GenomicsDB, we expect the former type of queries to be more frequent and hence, by default all arrays are stored in column major order (even when partitioning by rows across GenomicsDB instances).

   * Row major ordering implies that cells belonging to the same row are stored contiguously on disk.

* *Bulk importing and incremental importing*: *Bulk importing* of data implies that all the data is loaded at once into GenomicsDB and the array is never modified later on. *Incremental importing* implies that the array can be modified by adding new samples/CallSet over time. GenomicsDB support both modes of operation - however, repeated incremental updates (for example, incrementally adding one sample/CallSet at a time, ~1000 times) may lead to reduced performance during querying. If you are interested in why this occurs, please read the section on writing arrays to understand the concept of fragments and how updates are performed. Ideally, users should perform incremental imports a few times (~10s of times) per array and each import should contain a large number of samples.




.. 
   Architecture 
   *******************************
   <diagram showing stack and APIs>
..
   description




.. ---------------------------------------------------------------------

.. _Building + Installing: 

###############################
Building + Installing
###############################
Coming Soon

.. <copy exisiting info>

.. 
   copy information from GenomicsDB-Install once it is merged
   example with docker
   outline of environment variables

.. 
   Advanced
   *******************************
   performance tuning


.. ---------------------------------------------------------------------

###############################
APIs
###############################


Overview
*******************************
GenomicsDB has support/bindings for the following langauges.
<explanation of query jsons>

C++
*******************************
* native interface
* Querying CLI tools available - see :ref:`vcf2genomicsdb <vcf2genomicsdb>`


Python
*******************************
* Python bindings for GenomicsDB. Requires Python version *>= 3.6*
* Built using Cython.
* Currently, there is only support for querying (no importing).
* Data can be returned in either a numpy matrix or a list of lists.

The GenomicsDB python bindings are `published on PyPI`_.

.. _published on PyPI: https://pypi.org/project/genomicsdb/

Install ``genomicsdb`` with the following command:

.. code-block:: python

   pip install genomicsdb


API Reference
===============================
* Coming Soon. See the `bindings source code`_ for the time being. 

.. _bindings source code: https://github.com/nalinigans/GenomicsDB-Python

Java
*******************************
* java - variant context <link>
* Data can be returned in either a VariantContext object or <other>.
* The GenomicsDB java bindings are published to `Maven Central`_ every release.

.. _Maven Central: https://mvnrepository.com/artifact/org.genomicsdb/genomicsdb


* Add GenomicsDB to your maven project with the following snippet:

.. code-block:: java

   <!-- https://mvnrepository.com/artifact/org.genomicsdb/genomicsdb -->
   <dependency>
      <groupId>org.genomicsdb</groupId>
      <artifactId>genomicsdb</artifactId>
      <version>1.3.0</version>
   </dependency>


API Reference
===============================
* Coming Soon. See the `source code`_ for the time being. 

.. _source code: https://github.com/GenomicsDB/GenomicsDB/tree/master/src/main/java/org/genomicsdb

.. specific detail of GenomicsDB Importer usages

R
*******************************
* Support Coming Soon!

.. 
Spark
*******************************
.. <spark java - variant context>
   <explanation of how to arrange arrays + executors>
   <pom snippet>

.. ---------------------------------------------------------------------

.. _CLIs Tools:

###############################
CLI Tools
###############################
GenomicsDB contains several command line tools. These executables are not published online, 
you must :ref:`build and install <Building + Installing>` GenomicsDB to use them. 

.. _vcf2genomicsdb:

* vcf2genomicsdb
* gt mpi gather
* create genomicsdb workspace
* consolidate genomicsdb array
* histogram
* vcfdiff

.. ---------------------------------------------------------------------

.. _Import / ETL:

###############################
Import / ETL
###############################
GenomicsDB supports importing genomics data in several common formats, and with a variety of methods.

GenomicsDB can ingest data in VCF, gVCF. and CSV formats. When importing VCFs or gVCFs, 
you need to run a few preprocessing steps before importing. 

The primary and suggested method of importing is to use the native vcf2genomicsdb importer tool (see CLI Tools).
There is also a Java API for importing (see here).


VCF
*******************************
preprosessing steps

.. <copy exisiting docs>

gVCF
*******************************
preprosessing steps

.. <copy exisiting docs>

CSV
*******************************
Importing CSVs into GenomicsDB

Preliminaries
===============================
*  You need libcsv while `compiling`_.
.. _compiling: https://github.com/GenomicsDB/GenomicsDB/wiki/Compiling-GenomicsDB
*  License terms: We use libcsv to parse CSV files. 
   `libcsv`_ is licensed under the GNU Library or Lesser General Public License version 2.0 (LGPLv2). 
   So, if you are re-distributing binaries or object files, they may be subject to LGPLv2 terms. 
   Please ensure that any binaries/object files you distribute are compliant with LGPLv2. 
   You can disable libcsv usage by not setting the USE_LIBCSV and LIBCSV_DIR flags during compilation. 
   However, your binaries/executables will not be able to import CSV files into GenomicsDB.
.. _libcsv: https://sourceforge.net/projects/libcsv/

CSV format
===============================
Given a variant call at a specific position in the genome for a particular CallSet/sample, 
the CSV file format supported by GenomicsDB contains one line describing the call. 
Essentially, each line in the CSV describes one cell in the GenomicsDB array.

The exact format of the CSV is shown below:

.. code-block:: python

   <row>,<begin_column>,<end_column>,<REF>,<concatenated_ALT>,<QUAL>,<FILTER_field_specification>[,<other_fields_specification>]


Fixed fields:

* row (mandatory, type: int64): Row number in GenomicsDB for this sample/CallSet

* begin_column (mandatory, type: int64): Column in GenomicsDB at which the variant call begins

* end_column (mandatory, type: int64): Column in GenomicsDB at which the variant call ends (inclusive).

* REF (mandatory, type:string): Reference allele at the given position.

* concatenated_ALT (mandatory, type: string): Concatenated alternate alleles separated by the character '|'. For example, if a given call has two alternate alleles TAG and TG, then the concatenated string would be TAG|TG.

* QUAL (optional, type: float): Represents the QUAL field in a VCF. It can be left empty implying that the value was missing for this call.

* FILTER_field_specification (mandatory, type: variable length list of integers): Each element of the FILTER field list is an integer representing the FILTERs that this call failed (similar to how a BCF represents FILTERs). The first element of this field is an integer displaying the number of elements in the list. A value of 0 indicates that the list is empty.

Additional fields can be optionally specified

* <other_fields_specification>: The format depends on the type of field:

   * String type fields (fixed or variable length strings): The field should contain the string - an empty token indicates that the field is missing.

   * Fixed length field (int or float): The field should contain exactly N elements where N is the length of the field (fixed constant). One or more tokens may be left empty to indicate that those elements are missing.

   * Variable length field (int or float): The first element of this field should be an integer denoting the number of elements in the field for this call. It should then be followed by the elements of this field. An empty or missing field can be specified by setting the first element (field length) to 0 - no other elements should follow an empty field.


Example
-------------------------------
The following line contains 2 fields in addition to the fixed fields:

* SB: Fixed length field of 4 integers

* PL: Variable length field of integers

   .. code-block:: python

      2,1857210,1857210,G,A|T,894.77,0,,,,,6,923,0,599,996,701,1697

The line specifies the variant call for row id 2, beginning at column 1857210 and ending at 1857210. 
The REF allele is 'G' and the call has 2 alternate alleles 'A' and 'T' (SNVs). 
The QUAL value is 894.77 and there are no FILTERs specified (hence FILTER field length = 0). 
The SB field is missing - denoted by the 4 empty tokens. 
The PL field consists of 6 integers - the length appears first (since PL is a variable length field) followed by the elements [923,0,599,996,701,1697].


Special fields
-------------------------------
* GT is represented in the CSV as a variable length list of integers - each element of the list refers to the allele index (0 for reference allele, 1 for the first alternate allele and so on). The length of the list represents the ploidy of the sample/CallSet and must be specified in the CSV line (since GT is treated as a variable length list).


Organizing your data
===============================
* All CSV files imported into a GenomicsDB array must respect the number and order of fields as defined in the `vid_mapping_file`_.
.. _vid_mapping_file: https://github.com/GenomicsDB/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB#fields-information

* The import program cannot handle CSV files where 2 lines have the same value of row and begin_column - this restriction is similar to that imposed on `loading VCFs`_.
.. _loading VCFs: https://github.com/GenomicsDB/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB#fields-information


Compression Schemas
*******************************
GenomicsDB supports several compression schemas, which can be configured during build time. 

.. list-table:: Compression Schemas
   :widths: 25 25 50
   :header-rows: 1

   * - Compression Alogirthm
     - Build Variable Value
     - Notes
   * - None
     - 0
     - Default value if not specified
   * - Gzip
     - 1
     - Description
   * - Zstd
     - 2
     - Description
   * - Lz4
     - 3
     - Description


Best Practices
*******************************
* best practices

.. ---------------------------------------------------------------------

###############################
Supported Storage
###############################
GenomicsDB has support for several filesystems. 

Local / Cluster
*******************************
* *POSIX*: Tested on Centos 6 + 7 

* *HDFS*: HDFS support

Cloud
*******************************
* *AWS*: Support for S3 and EMRFS
* *GCP*: Support for G3
* *Azure*: Support for Azure Blob storage


.. ---------------------------------------------------------------------

###############################
Examples / Tutorials
###############################

Python Notebook
*******************************
Example Coming Soon

.. <code snippet querying example>

.. link out to full jupyter notebook code in github

Scala Notebook
*******************************
Example Coming Soon

.. <code snippet querying example>

.. <link out to full scala notebook code in github>

.. 
   Importing Data w/Docker
   *******************************
   <documentation from Nalini's repo>
.. ---------------------------------------------------------------------

###############################
Advanced
###############################

Performance Tuning
*******************************
* multi sample vs single sample
* array fragmentation
* tile sizes

Detailed Explanations
*******************************
.. in depth details of GenomicsDB that don't fit well elsewhere

.. ---------------------------------------------------------------------

###############################
Using with GATK
###############################

GenomicsDB is packaged into
`gatk4 <https://software.broadinstitute.org/gatk/documentation/article?id=11091>`__
and benefits qualitatively from a large user base.

.. ---------------------------------------------------------------------

.. 
   ###############################
   FAQs
   ###############################
   <faqs>

.. ---------------------------------------------------------------------
.. 
   ###############################
   Common Issues
   ###############################
   <issues>

.. ---------------------------------------------------------------------

###############################
Documentation WORK-IN-PROGRESS
###############################

|License: MIT|

GenomicsDB is a highly performant scalable data storage written in C++ for importing, querying and transforming genomic
variant data. 


Supported platforms and filesystems: 
*******************************

* Linux and MacOS. 
* POSIX, HDFS, EMRFS(S3), GCS and Azure Blob.

Included are
=============

* JVM/Spark wrappers that allow for streaming VariantContext_ buffers to/from the C++ layer among other functions. 
GenomicsDB jars with native libraries and only zlib dependencies are regularly published on `Maven Central`_.
* Native tools for incremental ingestion of variants in the form of VCF/BCF/CSV into GenomicsDB for performance.
* MPI and Spark support for parallel querying of GenomicsDB.

GenomicsDB is packaged into
`gatk4 <https://software.broadinstitute.org/gatk/documentation/article?id=11091>`__
and benefits qualitatively from a large user base.

The GenomicsDB documentation for users is hosted as a Github wiki:
https://github.com/GenomicsDB/GenomicsDB/wiki

GenomicsDB is open source and all participation is welcome. Please read
the `guidelines <contrib/README.md>`__ to help with contributions.

.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
.. _VariantContext: https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
.. .. _Maven Central: https://repo1.maven.org/maven2/org/genomicsdb/genomicsdb

