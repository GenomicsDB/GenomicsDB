/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.genomicsdb.spark;

import org.apache.spark.sql.types.StructType;

import java.util.Map;
import java.util.ArrayList;

/**
 * This interface declares some methods to be used by classes that define quantums of data that are
 * divided from a large amount of data to be distributed to many workers. Currently intended to be
 * implemeented by GenomicsDBInputSplit and GenomicsDBInputPartition
 */
public interface GenomicsDBInputInterface {
  default void setHost(String h) {}

  default void setGenomicsDBConf(GenomicsDBConfiguration g) {}

  default void setGenomicsDBSchema(StructType s) {}

  default void setGenomicsDBVidSchema(Map<String, GenomicsDBVidSchema> vMap) {}

  default void setPartitionInfo(GenomicsDBPartitionInfo p) {}

  default void setQueryInfoList(ArrayList<GenomicsDBQueryInfo> q) {}
}
