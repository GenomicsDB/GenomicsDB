.. java:import:: htsjdk.samtools.util Locatable

ChromosomeInterval
==================

.. java:package:: org.genomicsdb.importer.model
   :noindex:

.. java:type:: public class ChromosomeInterval implements Locatable

   Utility class to represent a chromosome interval Contains 3 members - chr name, start, end (1-based)

Constructors
------------
ChromosomeInterval
^^^^^^^^^^^^^^^^^^

.. java:constructor:: public ChromosomeInterval(String name, long begin, long end)
   :outertype: ChromosomeInterval

Methods
-------
getContig
^^^^^^^^^

.. java:method:: @Override public String getContig()
   :outertype: ChromosomeInterval

getEnd
^^^^^^

.. java:method:: @Override public int getEnd()
   :outertype: ChromosomeInterval

getStart
^^^^^^^^

.. java:method:: @Override public int getStart()
   :outertype: ChromosomeInterval

