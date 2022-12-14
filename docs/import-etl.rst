.. _Import / ETL:

###############################
Import / ETL
###############################
GenomicsDB supports importing genomics data in several common formats, and with a variety of methods.

GenomicsDB can ingest data in VCF, gVCF. and CSV formats. When importing VCFs or gVCFs, 
you need to run a few preprocessing steps before importing. 

The primary and suggested method of importing is to use the native vcf2genomicsdb importer tool (see CLI Tools).
There is also a Java API for importing (see here).


VCF/gVCF
*******************************
The import program can handle block compressed and indexed VCFs, gVCFs, BCFs and gBCFs. 
For brevity, we will only use the term VCF.

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


Organizing your data
===============================
* All your VCFs must be block compressed and indexed. `Bcftools`_ is one good option for compressing and indexing.

.. _Bcftools: https://github.com/samtools/bcftools

* The VCF format allows you to have multiple lines with the same position (identical chromosome+pos). The import program can NOT handle such VCFs. Make sure you use tools like bcftools to collapse multiple lines with the same position into a single line. Example commands are `provided here`_.

.. _provided here: https://github.com/GenomicsDB/GenomicsDB/wiki/Useful-external-tools#collapse-multiple-lines-with-the-same-genomic-position-into-a-single-line
:was
* In a multi-node environment, you must decide:
   * How to `partition your data in GenomicsDB`_.

   .. _partition your data in GenomicsDB: https://github.com/GenomicsDB/GenomicsDB/wiki/GenomicsDB-setup-in-a-multi-node-cluster

   * How your VCF files are accessed by the import program:
      * On a shared filesystem (NFS, Lustre etc) accessible from all nodes.
      *  If you are partitioning data by rows, your files can be scattered across local filesystems on multiple machines; each filesystem accessible only by the node on which it is mounted. Currently, such a setup isn't supported for column partitioning.




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




Best Practices
*******************************
* best practices
