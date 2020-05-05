/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 * Copyright (c) 2018-2019 Omics Data Automation, Inc.
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

#include "genomicsdb_GenomicsDBQueryStream.h"
#include "genomicsdb_bcf_generator.h"
#include "genomicsdb_export_config.pb.h"

#define GET_BCF_READER_FROM_HANDLE(X) (reinterpret_cast<GenomicsDBBCFGenerator*>(static_cast<std::uintptr_t>(X)))

void handleJNIException(JNIEnv *env, std::exception& exception) {
  std::string msg = std::string("GenomicsDB JNI Error: ") + exception.what();
  jclass io_exception_class = env->FindClass("java/io/IOException");
  if (io_exception_class) {
    jboolean flag = env->ExceptionCheck();
    if (flag) {
      env->ExceptionClear();
    }
    env->ThrowNew(io_exception_class, msg.c_str());
  } else {
    throw std::runtime_error(msg);
  }
}

JNIEXPORT jlong JNICALL Java_org_genomicsdb_reader_GenomicsDBQueryStream_jniGenomicsDBInit
  (JNIEnv* env, jobject curr_obj, jstring loader_configuration_file, 
   jbyteArray query_buffer, jstring chr, jlong start, jlong end,
   jint rank, jlong buffer_capacity, jlong segment_size,
   jboolean is_bcf, jboolean produce_header_only,
   jboolean use_missing_values_only_not_vector_end, jboolean keep_idx_fields_in_bcf_header)
{
  //Java string to char*
  auto loader_configuration_file_cstr = env->GetStringUTFChars(loader_configuration_file, NULL);
  auto chr_cstr = env->GetStringUTFChars(chr, NULL);
  // protobuf stuff from: https://askldjd.com/2013/02/19/protobuf-over-jni/
  genomicsdb_pb::ExportConfiguration query_config_pb;
  jbyte *bufferElems = env->GetByteArrayElements(query_buffer, 0);
  int len = env->GetArrayLength(query_buffer);
  auto output_format = is_bcf ? "bu" : "";
  GenomicsDBBCFGenerator *bcf_reader_obj;
  try {
    query_config_pb.ParseFromArray(reinterpret_cast<void*>(bufferElems), len);
    //Create object
    bcf_reader_obj = new GenomicsDBBCFGenerator(loader_configuration_file_cstr, &query_config_pb,
        chr_cstr, start, end,
        rank, buffer_capacity, segment_size, output_format,
        produce_header_only,
        is_bcf && use_missing_values_only_not_vector_end, is_bcf && keep_idx_fields_in_bcf_header);
  } catch (std::exception &exception) {
    handleJNIException(env, exception);
    bcf_reader_obj = NULL;
  }
  //Cleanup
  env->ReleaseStringUTFChars(loader_configuration_file, loader_configuration_file_cstr);
  env->ReleaseStringUTFChars(chr, chr_cstr);
  env->ReleaseByteArrayElements(query_buffer, bufferElems, JNI_ABORT);
  //Cast pointer to 64-bit int and return to Java
  return static_cast<jlong>(reinterpret_cast<std::uintptr_t>(bcf_reader_obj));
}

JNIEXPORT jlong JNICALL Java_org_genomicsdb_reader_GenomicsDBQueryStream_jniGenomicsDBClose
  (JNIEnv* env, jobject curr_obj, jlong handle)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  if(bcf_reader_obj) //not NULL
    delete bcf_reader_obj;
  return 0;
}

JNIEXPORT jlong JNICALL Java_org_genomicsdb_reader_GenomicsDBQueryStream_jniGenomicsDBGetNumBytesAvailable
  (JNIEnv* env, jobject curr_obj, jlong handle)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  return (bcf_reader_obj) ? bcf_reader_obj->get_buffer_capacity() : 0;
}

JNIEXPORT jbyte JNICALL Java_org_genomicsdb_reader_GenomicsDBQueryStream_jniGenomicsDBReadNextByte
  (JNIEnv* env, jobject curr_obj, jlong handle)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  return (bcf_reader_obj) ? bcf_reader_obj->read_next_byte() : -1;
}

JNIEXPORT jint JNICALL Java_org_genomicsdb_reader_GenomicsDBQueryStream_jniGenomicsDBRead
  (JNIEnv* env, jobject curr_obj, jlong handle, jbyteArray java_byte_array, jint offset, jint n)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  if(bcf_reader_obj == 0)
    return 0;
  auto total_num_bytes_read = 0ull;
  while(total_num_bytes_read < static_cast<uint64_t>(n) && !(bcf_reader_obj->end()))
  {
    auto& buffer_obj = bcf_reader_obj->get_read_batch();
    auto num_bytes_to_copy = std::min<size_t>(buffer_obj.get_num_remaining_bytes(), static_cast<size_t>(n)-total_num_bytes_read);
    //Handle this as a special case as read_and_advance will not advance anything if num_bytes_to_copy == 0u
    if(num_bytes_to_copy == 0u)
      num_bytes_to_copy = SIZE_MAX;     //forces jni_bcf_reader to produce the next batch of records
    else
    {
      env->SetByteArrayRegion(java_byte_array, offset+total_num_bytes_read, num_bytes_to_copy,
          reinterpret_cast<const jbyte*>(buffer_obj.get_pointer_at_read_position()));
      total_num_bytes_read += num_bytes_to_copy;
    }
    bcf_reader_obj->read_and_advance(0, 0u, num_bytes_to_copy);
  }
  return total_num_bytes_read;
}

JNIEXPORT jlong JNICALL Java_org_genomicsdb_reader_GenomicsDBQueryStream_jniGenomicsDBSkip
  (JNIEnv* env, jobject curr_obj, jlong handle, jlong n)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  return (bcf_reader_obj) ? bcf_reader_obj->read_and_advance(0, 0, n) : 0;
}
