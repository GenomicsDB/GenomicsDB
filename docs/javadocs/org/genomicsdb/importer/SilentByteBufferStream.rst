.. java:import:: java.io IOException

.. java:import:: java.io OutputStream

SilentByteBufferStream
======================

.. java:package:: org.genomicsdb.importer
   :noindex:

.. java:type::  class SilentByteBufferStream extends OutputStream

   Buffer stream implementation - it's silent in the sense that when the buffer is full, it doesn't raise an exception but just marks a a flag as full. It's up to the caller to check the flag and retry later Why? Most likely, it's faster to check a flag rather than throw and catch an exception

Constructors
------------
SilentByteBufferStream
^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public SilentByteBufferStream()
   :outertype: SilentByteBufferStream

   Constructor - uses default value of buffer capacity (20KiB)

SilentByteBufferStream
^^^^^^^^^^^^^^^^^^^^^^

.. java:constructor:: public SilentByteBufferStream(long capacity)
   :outertype: SilentByteBufferStream

   Constructor - uses specified buffer capacity

   :param capacity: size of buffer in bytes

Methods
-------
close
^^^^^

.. java:method:: @Override public void close() throws IOException
   :outertype: SilentByteBufferStream

flush
^^^^^

.. java:method:: @Override public void flush() throws IOException
   :outertype: SilentByteBufferStream

getBuffer
^^^^^^^^^

.. java:method:: public byte[] getBuffer()
   :outertype: SilentByteBufferStream

   Get byte buffer for this stream

   :return: byte buffer for this stream

getMarker
^^^^^^^^^

.. java:method:: public long getMarker()
   :outertype: SilentByteBufferStream

   Get marker value

   :return: marker value

getNumValidBytes
^^^^^^^^^^^^^^^^

.. java:method:: public long getNumValidBytes()
   :outertype: SilentByteBufferStream

   Get number of valid bytes

   :return: number of valid bytes

overflow
^^^^^^^^

.. java:method:: public boolean overflow()
   :outertype: SilentByteBufferStream

   Returns if the buffer has overflowed

   :return: true if the buffer has overflowed

resize
^^^^^^

.. java:method:: public void resize(long newSize)
   :outertype: SilentByteBufferStream

   Resizes buffer to new size - data is retained

   :param newSize: new capacity of the buffer

setMarker
^^^^^^^^^

.. java:method:: public void setMarker(long value)
   :outertype: SilentByteBufferStream

   Caller code can use this function to mark a certain point in the buffer This is generally used to mark the position in the buffer after the last complete VariantContext object written

   :param value: set marker value

setNumValidBytes
^^^^^^^^^^^^^^^^

.. java:method:: public void setNumValidBytes(long value)
   :outertype: SilentByteBufferStream

   Set number of valid bytes

   :param value: number of valid bytes

setOverflow
^^^^^^^^^^^

.. java:method:: public void setOverflow(boolean value)
   :outertype: SilentByteBufferStream

   Set overflow value

   :param value: overflow value

size
^^^^

.. java:method:: public int size()
   :outertype: SilentByteBufferStream

   Returns buffer capacity in bytes

   :return: buffer capacity in bytes

write
^^^^^

.. java:method:: @Override public void write(byte[] b, int off, int len) throws IOException
   :outertype: SilentByteBufferStream

write
^^^^^

.. java:method:: @Override public void write(byte[] b) throws IOException
   :outertype: SilentByteBufferStream

write
^^^^^

.. java:method:: @Override public void write(int b) throws IOException
   :outertype: SilentByteBufferStream

