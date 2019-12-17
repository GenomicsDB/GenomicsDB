[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

| Master | Develop |
| --- | --- |
| [![Travis](https://travis-ci.org/GenomicsDB/GenomicsDB.svg?branch=master)](https://travis-ci.org/GenomicsDB/GenomicsDB) | [![Travis](https://travis-ci.org/GenomicsDB/GenomicsDB.svg?branch=develop)](https://travis-ci.org/GenomicsDB/GenomicsDB?branch=develop) |
| [![codecov](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/master/graph/badge.svg)](https://codecov.io/gh/GenomicsDB/GenomicsDB) | [![codecov](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/develop/graph/badge.svg)](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/develop) |

GenomicsDB, originally from [Intel Health and Lifesciences](https://github.com/Intel-HLS/GenomicsDB), is built on top of a fork of [htslib](https://github.com/samtools/htslib) and a tile-based array storage system for importing, querying and transforming variant data. Variant data is sparse by nature (sparse relative to the whole genome) and using sparse array data stores is a perfect fit for storing such data. GenomicsDB is highly performant data storage written in C++ for importing, querying and transforming variant data.
* Supported platforms : Linux and MacOS.
* Supported filesystems : POSIX, HDFS, S3, GCS and Azure.

Included are
* JVM/Spark wrappers that allow for streaming [VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) buffers to/from the C++ layer among other functions. GenomicsDB jars with native libraries and only zlib dependencies are regularly published on [Maven Central](https://repo1.maven.org/maven2/org/genomicsdb/genomicsdb).
* Native tools for incremental ingestion of variants in the form of VCF/BCF/CSV into GenomicsDB for performance.
* MPI and Spark support for parallel querying of GenomicsDB.

GenomicsDB is packaged into [gatk4](https://software.broadinstitute.org/gatk/documentation/article?id=11091) and benefits qualitatively from a large user base.

The GenomicsDB documentation for users is hosted as a Github wiki:
https://github.com/GenomicsDB/GenomicsDB/wiki

GenomicsDB is open source and all participation is welcome. Please read the [guidelines](contrib/README.md) to help with contributions.
