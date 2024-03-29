
###############################
Overview
###############################

|License: MIT| 

GenomicsDB is a highly performant scalable data storage written in C++ for importing, querying and transforming genomic
variant data. 


.. _Concepts + Terminology:

Concepts + Terminology
*******************************
This section outlines basic GenomicsDB concepts and terminology used throughout this documentation.

* *Cell*: Cell of an array (for example: A[i][j] in a 2-D array)

* *Workspace*: A directory on the file system under which multiple GenomicsDB arrays can be stored. The underlying array stores some metadata files under the workspace directory. All GenomicsDB import/create programs assume that if the specified workspace directory exists, then it's a valid workspace directory. If the user points to an existing directory which doesn't contain the workspace metadata, the create/import programs will exit with an exception. If the directory doesn't exist, then the programs will initialize a new workspace correctly.

* *Array*: Name of a GenomicsDB array

* Given a workspace and array name, the framework will store its data in the directory *<workspace>/<array>*.

* *Column-major* and *row-major* ordering: Denotes the order in which cells are stored on disk.

   * Column major storage implies that cells belonging to the same column are stored contiguously on disk (cells belonging to the same column but different rows are sorted by row id). When cells are stored in column major order queries such as "retrieve all cells belonging to genomic position *X*" are fast due to high spatial locality. However, queries such as "retrieve cells for sample *Y*" are relatively slow. For GenomicsDB, we expect the former type of queries to be more frequent and hence, by default all arrays are stored in column major order (even when partitioning by rows across GenomicsDB instances).

   * Row major ordering implies that cells belonging to the same row are stored contiguously on disk.

* *Bulk importing and incremental importing*: *Bulk importing* of data implies that all the data is loaded at once into GenomicsDB and the array is never modified later on. *Incremental importing* implies that the array can be modified by adding new samples/CallSet over time. GenomicsDB support both modes of operation - however, repeated incremental updates (for example, incrementally adding one sample/CallSet at a time, ~1000 times) may lead to reduced performance during querying. If you are interested in why this occurs, please read the section on writing arrays to understand the concept of fragments and how updates are performed. Ideally, users should perform incremental imports a few times (~10s of times) per array and each import should contain a large number of samples.


.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT

