[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![readthedocs](https://readthedocs.org/projects/genomicsdb/badge/?version=latest)](https://genomicsdb.readthedocs.io/en/latest/)
[![Maven Central](https://img.shields.io/maven-central/v/org.genomicsdb/genomicsdb.svg)](https://mvnrepository.com/artifact/org.genomicsdb)

| Master | Develop |
| --- | --- |
| [![actions](https://github.com/GenomicsDB/GenomicsDB/workflows/build/badge.svg?branch=master)](https://github.com/GenomicsDB/GenomicsDB/actions?query=branch%3Amaster) | [![actions](https://github.com/GenomicsDB/GenomicsDB/workflows/build/badge.svg?branch=develop)](https://github.com/GenomicsDB/GenomicsDB/actions?query=branch%3Adevelop) |
| [![codecov](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/master/graph/badge.svg)](https://codecov.io/gh/GenomicsDB/GenomicsDB) | [![codecov](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/develop/graph/badge.svg)](https://codecov.io/gh/GenomicsDB/GenomicsDB/tree/develop) |

GenomicsDB is built on top of a fork of [htslib](https://github.com/samtools/htslib) and a tile-based array storage system for importing, querying and transforming variant data. Variant data is sparse by nature (sparse relative to the whole genome) and using sparse array data stores is a perfect fit for storing such data. GenomicsDB is a highly performant scalable data storage written in C++ for importing, querying and transforming genomic variant data. See [genomicsdb.readthedocs.io](https://genomicsdb.readthedocs.io/en/latest/) for documentation and usage.
* Supported platforms : Linux and MacOS.
* Supported filesystems : POSIX, HDFS, EMRFS(S3), GCS and Azure Blob.

Included are
* JVM/Spark wrappers that allow for streaming [VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) buffers to/from the C++ layer among other functions. GenomicsDB jars with native libraries and only zlib dependencies are regularly published on [Maven Central](https://repo1.maven.org/maven2/org/genomicsdb/genomicsdb).
* Native tools for incremental ingestion of variants in the form of VCF/BCF/CSV into GenomicsDB for performance.
* MPI and Spark support for parallel querying of GenomicsDB.

GenomicsDB is packaged into [gatk4](https://software.broadinstitute.org/gatk/documentation/article?id=11091) and benefits qualitatively from a large user base.

## External Contributions
GenomicsDB is open source and all participation is welcome.
GenomicsDB is released under the MIT License and all external
contributors are expected to grant an MIT License for their contributions.

### Checklist before creating Pull Request
Please ensure that the code is well documented in Javadoc style for Java/Scala.  For Java/C/C++ code formatting, roughly adhere to the Google Style Guides. See [GenomicsDB Style Guide](Style.md)

