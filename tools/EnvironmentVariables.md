## Environment Variables

The use of these Environment Variables will alter the behavior of GenomicsDB tools. These variables are meant to be experimental and may be removed in future when the functionality is main streamed.

* GENOMICSDB_SHARED_POSIXFS_OPTIMIZATIONS
     Relevant only for PosixFS. No read/write locks are maintained for GenomicsDB arrays. If used, it is the responsibility of the client to ensure that creating/updating/deleting arrays and array fragments are done with utmost care when shared.

The following variables can be used to experiment with different codecs and filters with tiles associated with GenomicsDB array fragments.
* GENOMICSDB_OFFSETS_DISABLE_DELTA_ENCODE
* GENOMICSDB_COORDS_DISABLE_DELTA_ENCODE
* GENOMICSDB_GT_ENABLE_BIT_SHUFFLE
* GENOMICSDB_GT_ENABLE_LZ4_COMPRESSION

* VCF2GENOMICSDB_INIT_PROGRESS_UPDATE_SAMPLE_SIZE
     Default for monitoring progress of merging headers is 2048 samples. Use this env to override default


