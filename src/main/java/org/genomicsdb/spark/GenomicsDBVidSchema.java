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

/**
 * The base class for a field from the vid mapping file.
 **/
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

  /** 
   * Set the appropriate types to the default schema.
   * Could consider reflection in the future. 
   * Include unboxing Double to float
   * @return string datatype for java/scala compat.
   **/
  public String getDataType(){
    String typeName;
    if (this.getFieldClass().equals(String.class)) {
      typeName = "String";
    } else if (this.getFieldClass().equals(Integer.class)) {
      typeName = "Integer";
    } 
    // TODO: switching to double below because that is what VC
    // seems to be giving us. Otherwise spark will complain later
    // when we try to create the dataframe about unboxing from double to float    
    else if (this.getFieldClass().equals(Double.class)) {
      typeName = "Double";
    } else if (this.getFieldClass().equals(Boolean.class)) {
      typeName = "Boolean";
    } else {
      throw new RuntimeException("Unsupported field type in vid!");
    }
    String dtypeDDL;
    // single item INFO fields
    if (this.isInfo() && (this.getLength().equals("1") || this.getFieldClass().equals(String.class))) {  
      dtypeDDL = typeName;
    }
    // multi-item INFO fields OR
    // single item FORMAT fields
    else if (this.isInfo() || this.getLength().equals("1") || this.getFieldClass().equals(String.class)) {
      dtypeDDL = "ARRAY<" + typeName + ">";
    }
    // multi-item FORMAT fields
    else {
      dtypeDDL = "ARRAY<ARRAY<" + typeName + ">>";
    }
    return dtypeDDL;
  }

}
