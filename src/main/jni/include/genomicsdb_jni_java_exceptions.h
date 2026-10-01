/**
 * The MIT License (MIT)
 * Copyright (c) 2026 Tim Fennell
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

#ifndef GENOMICSDB_JNI_JAVA_EXCEPTIONS_H
#define GENOMICSDB_JNI_JAVA_EXCEPTIONS_H

#include <jni.h>

#include <exception>
#include <string>

/**
 * Throws an org.genomicsdb.exception.GenomicsDBException in the JVM with the message of the C++ exception being
 * handled. A Java exception that is already pending, e.g. from a failed JNI call, is left as the one Java sees.
 * Call only from a catch block.
 */
inline void throw_current_exception_to_java(JNIEnv* env) {
  if (env->ExceptionCheck()) {
    return;
  }
  std::string message;
  try {
    throw;
  } catch (const std::exception& exception) {
    message = exception.what();
  } catch (...) {
    message = "Unknown native exception";
  }
  auto exception_class = env->FindClass("org/genomicsdb/exception/GenomicsDBException");
  if (exception_class) {
    env->ThrowNew(exception_class, message.c_str());
  }
}

/**
 * Runs the body of a JNI entry point and turns any C++ exception it throws into a Java exception, since a C++
 * exception unwinding out of native code into the JVM aborts the whole process. Returns what the body returns or,
 * after an exception, a value-initialised result that Java never sees, as the Java exception is pending.
 */
template<typename Body>
auto with_java_exceptions(JNIEnv* env, Body&& body) -> decltype(body()) {
  using Result = decltype(body());
  try {
    return body();
  } catch (...) {
    throw_current_exception_to_java(env);
    return Result();
  }
}

#endif
