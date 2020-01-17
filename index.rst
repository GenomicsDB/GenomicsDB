Documentation WORK-IN-PROGRESS
*******************************

|License: MIT|

GenomicsDB, originally from `Intel Health and
Lifesciences <https://github.com/Intel-HLS/GenomicsDB>`__, is built on
top of a fork of `htslib <https://github.com/samtools/htslib>`__ and a
tile-based array storage system for importing, querying and transforming
variant data. Variant data is sparse by nature (sparse relative to the
whole genome) and using sparse array data stores is a perfect fit for
storing such data. GenomicsDB is a highly performant scalable data
storage written in C++ for importing, querying and transforming genomic
variant data. 

Supported platforms and filesystems: 
-------------------------------------

* Linux and MacOS. 
* POSIX, HDFS, EMRFS(S3), GCS and Azure Blob.

Included are
-------------

* JVM/Spark wrappers that allow for streaming
`VariantContext <https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html>`
buffers to/from the C++ layer among other functions. GenomicsDB jars
with native libraries and only zlib dependencies are regularly published
on `Maven
Central <https://repo1.maven.org/maven2/org/genomicsdb/genomicsdb>`.
* Native tools for incremental ingestion of variants in the form of
VCF/BCF/CSV into GenomicsDB for performance.
* MPI and Spark support
for parallel querying of GenomicsDB.

GenomicsDB is packaged into
`gatk4 <https://software.broadinstitute.org/gatk/documentation/article?id=11091>`__
and benefits qualitatively from a large user base.

The GenomicsDB documentation for users is hosted as a Github wiki:
https://github.com/GenomicsDB/GenomicsDB/wiki

GenomicsDB is open source and all participation is welcome. Please read
the `guidelines <contrib/README.md>`__ to help with contributions.

.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
