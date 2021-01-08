.. java:import:: org.apache.spark.sql.sources.v2 ReadSupport

.. java:import:: org.apache.spark.sql.sources.v2.reader DataSourceReader

.. java:import:: org.apache.spark.sql.sources.v2 DataSourceOptions

.. java:import:: org.apache.spark.sql.types StructType

.. java:import:: org.apache.spark.sql.sources DataSourceRegister

GenomicsDBDataSourceV2
======================

.. java:package:: org.genomicsdb.spark
   :noindex:

.. java:type:: public class GenomicsDBDataSourceV2 implements ReadSupport, DataSourceRegister

Constructors
------------
GenomicsDBDataSourceV2
^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public GenomicsDBDataSourceV2()
   :outertype: GenomicsDBDataSourceV2

Methods
-------
createReader
^^^^^^^^^^^^

.. java:method:: public DataSourceReader createReader(DataSourceOptions options)
   :outertype: GenomicsDBDataSourceV2

createReader
^^^^^^^^^^^^

.. java:method:: public DataSourceReader createReader(StructType schema, DataSourceOptions options)
   :outertype: GenomicsDBDataSourceV2

shortName
^^^^^^^^^

.. java:method:: public String shortName()
   :outertype: GenomicsDBDataSourceV2

