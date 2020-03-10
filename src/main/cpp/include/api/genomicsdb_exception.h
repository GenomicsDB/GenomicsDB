/**
 * @file genomicsdb_exception.h
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
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
 *
 * @section DESCRIPTION
 *
 * GenomicsDBException is a catch all for all underlying genomicsdb library exceptions.
 *
 **/

#ifndef GENOMICSDB_EXCEPTION_H
#define GENOMICSDB_EXCEPTION_H

#include <exception>
#include <string>

// Override project visibility set to hidden for api
#if (defined __GNUC__ && __GNUC__ >= 4) || defined __INTEL_COMPILER
#  define GENOMICSDB_EXCEPTION_EXPORT __attribute__((visibility("default")))
#else
#  define GENOMICSDB_EXCEPTION_EXPORT
#endif

class GenomicsDBException : public std::exception {
 public:
  GENOMICSDB_EXCEPTION_EXPORT GenomicsDBException(const std::string m="GenomicsDB Exception: ") : m_msg(m) {}
  GENOMICSDB_EXCEPTION_EXPORT ~GenomicsDBException() {}
  /** Returns the exception message. */
  GENOMICSDB_EXCEPTION_EXPORT const char* what() const noexcept {
    return m_msg.c_str();
  }
 private:
  std::string m_msg;
};

#endif /* GENOMICSDB_EXCEPTION_H */
