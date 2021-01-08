.. java:import:: htsjdk.variant.variantcontext VariantContext

.. java:import:: org.apache.hadoop.conf Configuration

.. java:import:: org.apache.spark SparkConf

.. java:import:: org.apache.spark.api.java JavaPairRDD

.. java:import:: org.apache.spark.api.java JavaRDD

.. java:import:: org.apache.spark.api.java JavaSparkContext

.. java:import:: java.util List

GenomicsDBJavaSparkFactory
==========================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public final class GenomicsDBJavaSparkFactory

   This factory class exposes how a JavaRDD of variant contexts (htsjdk) can be retrieved from GenomicsDB. In case of the newAPIHadoopRDD(), GenomicsDB returns a JavaPairRDD where the genomics positions are the key. However, this is seldom used in the variant contexts as downstream applications in HellBender code uses only the values and ignores the key

Methods
-------
main
^^^^

.. java:method:: public static void main(String[] args)
   :outertype: GenomicsDBJavaSparkFactory

usingNewAPIHadoopRDD
^^^^^^^^^^^^^^^^^^^^

.. java:method:: @SuppressWarnings public static void usingNewAPIHadoopRDD(String[] args)
   :outertype: GenomicsDBJavaSparkFactory

