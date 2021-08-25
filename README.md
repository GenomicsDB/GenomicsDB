[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Maven Central](https://img.shields.io/maven-central/v/org.genomicsdb/genomicsdb.svg)](https://mvnrepository.com/artifact/org.genomicsdb)

| Master | Develop |
| --- | --- |
| [![actions](https://github.com/GenomicsDB/GenomicsDB/workflows/build/badge.svg?branch=master)](https://github.com/GenomicsDB/GenomicsDB/actions?query=branch%3Amaster) | [![actions](https://github.com/GenomicsDB/GenomicsDB/workflows/build/badge.svg?branch=develop)](https://github.com/GenomicsDB/GenomicsDB/actions?query=branch%3Adevelop) |
| [![codecov](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/master/graph/badge.svg)](https://codecov.io/gh/GenomicsDB/GenomicsDB) | [![codecov](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/develop/graph/badge.svg)](https://codecov.io/gh/GenomicsDB/GenomicsDB/branch/develop) |

GenomicsDB, originally from [Intel Health and Lifesciences](https://github.com/Intel-HLS/GenomicsDB), is built on top of a fork of [htslib](https://github.com/samtools/htslib) and a tile-based array storage system for importing, querying and transforming variant data. Variant data is sparse by nature (sparse relative to the whole genome) and using sparse array data stores is a perfect fit for storing such data. GenomicsDB is a highly performant scalable data storage written in C++ for importing, querying and transforming genomic variant data.
* Supported platforms : Linux and MacOS.
* Supported filesystems : POSIX, HDFS, EMRFS(S3), GCS and Azure Blob.

Included are
* JVM/Spark wrappers that allow for streaming [VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) buffers to/from the C++ layer among other functions. GenomicsDB jars with native libraries and only zlib dependencies are regularly published on [Maven Central](https://repo1.maven.org/maven2/org/genomicsdb/genomicsdb).
* Native tools for incremental ingestion of variants in the form of VCF/BCF/CSV into GenomicsDB for performance.
* MPI and Spark support for parallel querying of GenomicsDB.

GenomicsDB is packaged into [gatk4](https://software.broadinstitute.org/gatk/documentation/article?id=11091) and benefits qualitatively from a large user base.

The GenomicsDB documentation for users is hosted as a Github wiki:
https://github.com/GenomicsDB/GenomicsDB/wiki

## External Contributions
GenomicsDB is open source and all participation is welcome.
GenomicsDB is released under the MIT License and all external
contributors are expected to grant an MIT License for their contributions.

### Checklist before creating Pull Request
Please ensure that the code is well documented in Javadoc style for Java/Scala or roughly adhere to [Google C++ Style](https://google.github.io/styleguide/cppguide.html) for C/C++ for consistency/readabilty.

```
Use spaces instead of tabs.
Use 2 spaces for indenting.
Add brackets even for one line blocks e.g. 
        if (x>0)
           do_foo();
 should ideally be 
       if (x>0) {
         do_foo();
       }
Pad header e.g.
        if(x>0) should be if (x>0)
        while(x>0) should be while (x>0)
One half indent for class modifiers.
```
