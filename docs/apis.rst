
.. _API Docs: 

###############################
APIs
###############################


Overview
*******************************
In addition to the native C++ interface, GenomicsDB has support/bindings for Python, Java, and R.
The data structure returned from a GenomicsDB query depends on the API or CLI tool being used and the method invoked. 
Similarly, the method of passing in query field information depends on the route chosen.

The various APIs take in the same query fields in one form or another, so check the appropriate reference for your use case. 
If you are new to GenomicsDB, we recommend starting with the :ref:`gt_mpi_gather <CLI Tools gt-mpi-gather>` CLI tool to familiarize yourself with the concepts. 

All the APIs support taking configuration options as protobuf. Here's the protobuf documentation:

.. toctree::
   :maxdepth: 2

   protodoc/proto.md


C++
*******************************

* Native interface - See `cpp source code`_
* Querying CLI tools available - see :ref:`vcf2genomicsdb <CLI Tools vcf2genomicsdb>`

.. _cpp source code: https://github.com/GenomicsDB/GenomicsDB/tree/master/src/main/cpp

API Reference
===============================

.. toctree::
   :maxdepth: 2

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

.. toctree::
   :maxdepth: 2

   python/genomicsdb

Java
*******************************

* Data can be returned in either a htsjdk `VariantContext`_ object or <other>.
* The GenomicsDB java bindings are published to `Maven Central`_ every release.
* Add GenomicsDB to your maven project with the following snippet (update the version, if necessary):

.. code-block:: java

   <!-- https://mvnrepository.com/artifact/org.genomicsdb/genomicsdb -->
   <dependency>
      <groupId>org.genomicsdb</groupId>
      <artifactId>genomicsdb</artifactId>
      <version>1.4.5</version>
   </dependency>

.. _VariantContext: https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
.. _Maven Central: https://mvnrepository.com/artifact/org.genomicsdb/genomicsdb


API Reference
===============================

.. toctree::
   :maxdepth: 2

   doxyrst/javaimporter
   doxyrst/javareader





