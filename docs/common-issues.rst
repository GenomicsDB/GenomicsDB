###############################
Common Issues
###############################

1.  Verify that your JSON configuration files are syntactically correct before invoking GenomicsDB tools

    The GenomicsDB tools assume that you provide a syntactically correct JSON - always validate your JSON files using tools such as *json_verify* before invoking GenomicsDB tools.

    .. code-block:: none  

        cat <json_file> | json_verify

.. ------------------------------

2.  If you are using GenomicsDB version < 0.8.0, you cannot move/copy the GenomicsDB workspace/array directory arbitrarily after loading. If you do copy/move files to a different machine, ensure that the absolute path of the workspace/array directory is the same as that on the machine on which data is loaded/imported. This limitation is removed in version 0.8.0.

.. ------------------------------

3.  I have prepared my JSON files correctly, yet I get an exception with the following message:
    
    .. code-block:: none 

        Could not open vid mapping file "~/<directory>/vid_mapping_file.json OR
        Could not open callsets file "~/<directory>/callset_mapping_file.json


    The character '~' is interpreted by shells (bash, tcsh) as the home directory of the user - it is NOT interpreted by file I/O system calls invoked by the GenomicsDB executables/libraries. Hence, you must specify the path of the file without special characters that are interpreted by shells.
  
    Correct examples:

    .. code-block:: none  

        /home/<user>/<directory>/vid_mapping_file.json OR
        ../../<directory>/vid_mapping_file.json      

.. ------------------------------

4.  I get an exception with the following message:

    .. code-block:: none  

        Unhandled overlapping variants at columns <col1> and <col2> for row <row>


    GenomicsDB cannot deal with overlapping variants within a single sample - `more details documented here`_. Workarounds exist for dealing with overlapping deletions and gVCF reference blocks (intervals with \<NON_REF\> as the only alternate allele). When overlapping variants which are neither deletions nor reference blocks are found in the input VCF file, then the above exception is thrown. Check whether bcftools can help you - see point 2 `in the section on organizing your data`_.

    .. _more details documented here: https://github.com/GenomicsDB/GenomicsDB/wiki/Overlapping-variant-calls-in-a-sample
    .. _in the section on organizing your data: https://github.com/GenomicsDB/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB#organizing-your-data

.. ------------------------------

5.  I have both *row_partitions* and *column_partitions* in my loader JSON file and I get an exception message:

    .. code-block:: none  

        Cannot have both "row_partitions" and "column_partitions" simultaneously in the JSON file


    A GenomicsDB array can be partitioned by rows or columns but not both simultaneously - see `this page`_ for more information.

    .. _this page: https://github.com/GenomicsDB/GenomicsDB/wiki/GenomicsDB-setup-in-a-multi-node-cluster

.. ------------------------------

6.  I have setup all my JSON files correctly, but the import program finishes almost immediately without importing any data from my VCFs:

    There could be many reasons, but here are the common issues we have seen users running into:

    * *Contig/chromosome names don't match in the vid_mapping_file and the input VCFs*: The contig/chromosome names in the VCF and the vid mapping JSON file MUST match EXACTLY. For example, if the vid file has a contig named _"1"_ (as per the 1000 genomes naming convention) while the VCF has a contig named _"chr1"_ (as per the UCSC convention), GenomicsDB will ignore all data corresponding to _"chr1"_.

.. ------------------------------

7.  I have setup all my JSON files correctly, but the import program doesn't load data for some of the samples:
  
    There could be many reasons, but here are the common issues we have seen users running into:

    * *Incorrect value(s) of idx_in_file in the callset_mapping_file*: Note that *row_idx* is the globally unique value of the TileDB row index corresponding to a given sample/CallSet. *idx_in_file* is useful mostly for multi-sample VCFs and specifies the index of the sample in a given VCF. For single sample VCFs, this field should be 0 (or omitted altogether).

    * *Incorrect partition bounds in the loader JSON or incorrectly specified partition index in the command line*: Read the section on `running the program`_ in the `import data`_ wiki section. Also, please re-check your partition bounds in the loader JSON.

    .. _running the program: https://github.com/GenomicsDB/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB#running-the-program
    .. _import data: https://github.com/GenomicsDB/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB

.. ------------------------------

8.  I see an incorrect cell order found error as:

    .. code-block:: bash

        $ vcf2genomicsdb loader.json 
        terminate called after throwing an instance of 'VCF2GenomicsDBException'
        what():  VCF2GenomicsDBException : Incorrect cell order found - cells must be in column major order. Previous cell: [ 0, 114111 ] current cell: [ 0, 114111 ] 
        Aborted
    

    The error occurs if alleles at the same position span across multiple lines, for example:

    .. code-block:: none  

        chrX	114112	.	TCT	T	999	PASS	.	GT:DP:GQ:MIN_DP:PL	0/0:0:0:0:0,0,0
        chrX	114112	.	TCT	TTT	999	PASS	.	GT:DP:GQ:MIN_DP:PL	0/0:0:0:0:0,0,0


    The fix is to run bcftools norm as described `here`_ which will merge the alleles as

    .. _here: https://github.com/Intel-HLS/GenomicsDB/wiki/Useful-external-tools#useful-bcftools-commands

    .. code-block:: none  

        chrX	114112	.	TCT	T,TTT	999	PASS	.	GT:DP:GQ:MIN_DP:PL	0/2:0:0:0:0,0,0,0,0,0

.. ------------------------------

9.  I see an error message:

    .. code-block:: none  

        Cannot open VCF/BCF file <path.vcf.gz>


    even when the file and its index exist. What's going on?

    If you are importing data from many files (\>1000), then it's likely that you are hitting the limit on the number of open files set in your machine(s). Find out how to increase the limit.