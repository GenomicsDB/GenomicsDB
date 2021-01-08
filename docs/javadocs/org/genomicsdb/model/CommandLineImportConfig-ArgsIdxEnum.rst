.. java:import:: gnu.getopt Getopt

.. java:import:: gnu.getopt LongOpt

.. java:import:: htsjdk.tribble AbstractFeatureReader

.. java:import:: htsjdk.tribble FeatureReader

.. java:import:: htsjdk.tribble.readers LineIterator

.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: htsjdk.variant.vcf VCFCodec

.. java:import:: htsjdk.variant.vcf VCFHeader

.. java:import:: htsjdk.variant.vcf VCFUtils

.. java:import:: java.util.function Consumer

.. java:import:: java.util.stream IntStream

.. java:import:: java.io IOException

.. java:import:: java.net URI

.. java:import:: java.net URISyntaxException

CommandLineImportConfig.ArgsIdxEnum
===================================

.. java:package:: org.genomicsdb.model
   :noindex:

.. java:type:: public enum ArgsIdxEnum
   :outertype: CommandLineImportConfig

Enum Constants
--------------
ARGS_IDX_AFTER_LAST_ARG_IDX
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_AFTER_LAST_ARG_IDX
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_BATCHSIZE
^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_BATCHSIZE
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_CALLSET_OUTPUT
^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_CALLSET_OUTPUT
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_FAIL_IF_UPDATING
^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_FAIL_IF_UPDATING
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_INCREMENTAL_IMPORT
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_INCREMENTAL_IMPORT
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_PASS_AS_BCF
^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_PASS_AS_BCF
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_READ_INPUT_VCF_USING_HTSLIB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_READ_INPUT_VCF_USING_HTSLIB
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_SEGMENT_SIZE
^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_SEGMENT_SIZE
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_SIZE_PER_COLUMN_PARTITION
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_SIZE_PER_COLUMN_PARTITION
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_USE_SAMPLES_IN_ORDER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_USE_SAMPLES_IN_ORDER
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_VCF_HEADER_OUTPUT
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_VCF_HEADER_OUTPUT
   :outertype: CommandLineImportConfig.ArgsIdxEnum

ARGS_IDX_VIDMAP_OUTPUT
^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: public static final CommandLineImportConfig.ArgsIdxEnum ARGS_IDX_VIDMAP_OUTPUT
   :outertype: CommandLineImportConfig.ArgsIdxEnum

Methods
-------
idx
^^^

.. java:method:: public int idx()
   :outertype: CommandLineImportConfig.ArgsIdxEnum

