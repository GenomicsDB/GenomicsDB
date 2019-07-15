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
#ifndef GENOMICSDB_VID_MAPPER_PB_H
#define GENOMICSDB_VID_MAPPER_PB_H

#include "vid_mapper.h"
#include "genomicsdb_vid_mapping.pb.h"
#include "genomicsdb_callsets_mapping.pb.h"
#include <google/protobuf/stubs/common.h>

enum {
  GENOMICSDB_VID_MAPPER_SUCCESS = 0x0,
  GENOMICSDB_VID_MAPPER_FAILURE = 0x1
};

/*
 * Constructor is called for the global variable when the library is loaded
 * Destructor is called when the library is unloaded;
 */
class GenomicsDBProtoBufInitAndCleanup {
 public:
  GenomicsDBProtoBufInitAndCleanup() {
    // Verify that the version of the library that we linked against is
    // compatible with the version of the headers we compiled against.
    GOOGLE_PROTOBUF_VERIFY_VERSION;
  }
  ~GenomicsDBProtoBufInitAndCleanup() {
    GenomicsDBProtoBufInitAndCleanup::shutdown_protobuf_library();
  }
  static void shutdown_protobuf_library() {
    google::protobuf::ShutdownProtobufLibrary();
  }
};
extern GenomicsDBProtoBufInitAndCleanup g_genomicsdb_protobuf_init_and_cleanup;

class ProtoBufBasedVidMapperException : public std::exception {
 public:
  /**
   * Default constructor
   */
  ProtoBufBasedVidMapperException(const std::string m="");

  /*
   * Destructor: must be called before end of scope of a
   * VidMapper exception object
   */
  ~ProtoBufBasedVidMapperException();

  /*
   * ACCESSORS
   */

  /**
   * Returns the exception message.
   */
  const char* what() const noexcept {
    return msg_.c_str();
  }

 private:
  std::string msg_;
};

#endif  // GENOMICSDB_VID_MAPPER_PB_H
