/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 * Copyright (c) 2020 Omics Data Automation, Inc.
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

#include "jni_mpi_init.h"
#include "genomicsdb_GenomicsDBLibLoader.h"

JNIMpiInit g_jni_mpi_init;

std::string get_system_property(JNIEnv* env, const std::string& name ) {
  jclass java_system_class = env->FindClass("java/lang/System");
  jmethodID java_property_method = env->GetStaticMethodID(java_system_class, "getProperty", "(Ljava/lang/String;)Ljava/lang/String;");
  jstring java_property_name = env->NewStringUTF(name.c_str());
  jstring java_property_value = static_cast<jstring>(env->CallStaticObjectMethod(java_system_class, java_property_method, java_property_name ));
  if (java_property_value) {
    const char* property = env->GetStringUTFChars(java_property_value, 0);
    std::string property_string = property?std::string(property):"";
    env->ReleaseStringUTFChars(java_property_value, property);
    return property_string;
  }
  return "";
}

JNIEXPORT jint JNICALL Java_org_genomicsdb_GenomicsDBLibLoader_jniGenomicsDBOneTimeInitialize
  (JNIEnv * env, jclass obj)
{
  std::string gatk_stacktrace_on_user_exception = get_system_property(env, "GATK_STACKTRACE_ON_USER_EXCEPTION");
  if (!gatk_stacktrace_on_user_exception.empty()) {
    setenv("GENOMICSDB_PRINT_STACKTRACE", gatk_stacktrace_on_user_exception.c_str(), 1);
  }
  g_jni_mpi_init.initialize();
  return 0;
}
