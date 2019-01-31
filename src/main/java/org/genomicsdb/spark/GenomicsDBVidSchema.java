/*
 * The MIT License (MIT)
 * Copyright (c) 2019 Omics Data Automation
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

import java.io.Serializable;

public class GenomicsDBVidSchema implements Serializable {

  boolean isInfo;
  Class<?> clazz;
  String length;

  public GenomicsDBVidSchema(boolean isInfo, Class<?> clazz, String length) {
    this.isInfo = isInfo;
    this.clazz = clazz;
    this.length = length;
  }

  public boolean isInfo() {
    return isInfo;
  }

  public Class<?> getFieldClass() {
    return clazz;
  }

  public String getLength() {
    return length;
  }
}
