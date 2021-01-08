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

CommandLineImportConfig
=======================

.. java:package:: org.genomicsdb.model
   :noindex:

.. java:type:: public class CommandLineImportConfig extends ImportConfig

Fields
------
DEFAULT_SIZE_PER_COLUMN_PARTITION
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. java:field:: protected static final long DEFAULT_SIZE_PER_COLUMN_PARTITION
   :outertype: CommandLineImportConfig

Constructors
------------
CommandLineImportConfig
^^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public CommandLineImportConfig(String command, String[] commandArgs)
   :outertype: CommandLineImportConfig

