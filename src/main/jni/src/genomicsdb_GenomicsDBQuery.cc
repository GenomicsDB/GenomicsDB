/**
 * The MIT License (MIT)
 * Copyright (c) 2019 Omics Data Automation, Inc.
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

#include "genomicsdb_GenomicsDBQuery.h"
#include "genomicsdb.h"

#include <mutex>

bool is_initialized_ = false;
std::mutex initializing_mutex_;

//java.util.ArrayList
jclass java_ArrayList_ ;
jmethodID java_ArrayList_init_;
jmethodID java_ArrayList_size_;
jmethodID java_ArrayList_get_;
jmethodID java_ArrayList_add_;

//java.util.Hashmap
jclass java_HashMap_;
jmethodID java_HashMap_init_;
jmethodID java_HashMap_put_;

//org.genomicsdb.reader.GenomicsDBQuery$VariantCalls
jclass java_VariantCalls_;
jmethodID java_VariantCalls_init_default_;
jmethodID java_VariantCalls_init_;

//org.genomicsdb.reader.GenomicsDBQuery$VariantCall
jclass java_VariantCall_;
jmethodID java_VariantCall_init_;

//org.genomicsdb.reader.GenomicsDBQuery$Pair
jclass java_Pair_;
jmethodID java_Pair_init_;
jmethodID java_Pair_getStart_;
jmethodID java_Pair_getEnd_;

#define INIT(VAR,CODE)                                                                                 \
  do {                                                                                                 \
    VAR = CODE;                                                                                        \
    if (VAR == NULL) {                                                                                 \
      throw GenomicsDBException("genomicsdb_GenomicsDBQuery.cc#initialize:"+std::to_string(__LINE__)); \
    }                                                                                                  \
  } while (false)

void initialize(JNIEnv *env) {
  std::lock_guard<std::mutex> guard(initializing_mutex_);
  if (!is_initialized_) {
    //java.util.ArrayList
    INIT(java_ArrayList_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList"))));
    INIT(java_ArrayList_init_, env->GetMethodID(java_ArrayList_, "<init>", "()V"));
    INIT(java_ArrayList_size_, env->GetMethodID (java_ArrayList_, "size", "()I"));
    INIT(java_ArrayList_get_, env->GetMethodID(java_ArrayList_, "get", "(I)Ljava/lang/Object;"));
    INIT(java_ArrayList_add_, env->GetMethodID(java_ArrayList_, "add", "(Ljava/lang/Object;)Z"));

    //java.util.Hashmap
    INIT(java_HashMap_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/HashMap"))));
    INIT(java_HashMap_init_, env->GetMethodID(java_HashMap_, "<init>", "()V"));
    INIT(java_HashMap_put_, env->GetMethodID(java_HashMap_, "put", "(Ljava/lang/Object;Ljava/lang/Object;)Ljava/lang/Object;"));
    
    //org.genomicsdb.reader.GenomicsDBQuery$VariantCalls
    INIT(java_VariantCalls_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("org/genomicsdb/reader/GenomicsDBQuery$VariantCalls"))));
    INIT(java_VariantCalls_init_default_,  env->GetMethodID(java_VariantCalls_, "<init>", "()V"));
    INIT(java_VariantCalls_init_,  env->GetMethodID(java_VariantCalls_, "<init>", "(JJ)V"));

    //org.genomicsdb.reader.GenomicDBQuery$VariantCall
    INIT(java_VariantCall_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("org/genomicsdb/reader/GenomicsDBQuery$VariantCall"))));
    INIT(java_VariantCall_init_, env->GetMethodID(java_VariantCall_, "<init>", "(ILjava/lang/String;JJLjava/util/Map;)V"));

    //org.genomicsdb.reader.GenomicsDBQuery$Pair<L, R>
    INIT(java_Pair_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("org/genomicsdb/reader/GenomicsDBQuery$Pair"))));
    INIT(java_Pair_init_, env->GetMethodID(java_Pair_, "<init>", "(JJ)V"));
    INIT(java_Pair_getStart_, env->GetMethodID(java_Pair_, "getStart", "()J"));
    INIT(java_Pair_getEnd_, env->GetMethodID(java_Pair_, "getEnd", "()J"));
    
    is_initialized_ = true;
  }
}

void JNI_OnUnload(JavaVM *vm, void *reserved) {
  std::lock_guard<std::mutex> guard(initializing_mutex_);
  if (is_initialized_) {
    JNIEnv *env;
    if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION_1_8) == JNI_OK) {
      env->DeleteGlobalRef(java_Pair_);
      env->DeleteGlobalRef(java_VariantCall_);
      env->DeleteGlobalRef(java_VariantCalls_);
      env->DeleteGlobalRef(java_HashMap_);
      env->DeleteGlobalRef(java_ArrayList_);

      is_initialized_ = false;
    }
  }
}

#if(0)
# Useful for debugging JNI
void get_class_name(JNIEnv *env, jobject obj) {
  jclass cls = env->GetObjectClass(obj);
  jmethodID mid = env->GetMethodID(cls, "getClass", "()Ljava/lang/Class;");
  jobject clsObj = env->CallObjectMethod(obj, mid);
  cls = env->GetObjectClass(clsObj);
  mid = env->GetMethodID(cls, "getName", "()Ljava/lang/String;");
  jstring strObj = (jstring)env->CallObjectMethod(clsObj, mid);
  const char* str = env->GetStringUTFChars(strObj, NULL);
  // Print the class name
  printf("\nCalling class is: %s\n", str);
  env->ReleaseStringUTFChars(strObj, str);
}
#endif

std::vector<std::string> to_string_vector(JNIEnv *env, jobject arrayList) {
  initialize(env);
  jint size = env->CallIntMethod(arrayList, java_ArrayList_size_);
  std::vector<std::string> result;
  result.reserve(size);
  for (jint i=0; i<size; i++) {
    jstring element = static_cast<jstring>(env->CallObjectMethod(arrayList, java_ArrayList_get_, i));
    auto element_cstr = env->GetStringUTFChars(element, nullptr);
    result.emplace_back(element_cstr);
    env->ReleaseStringUTFChars(element, element_cstr);
    env->DeleteLocalRef(element);
  }
  return std::move(result);
}

genomicsdb_ranges_t to_genomicsdb_ranges_vector(JNIEnv *env, jobject arrayList) {
  jint size = env->CallIntMethod(arrayList, java_ArrayList_size_);
  genomicsdb_ranges_t result;
  result.reserve(size);
  for (jint i=0; i<size; i++) {
    jobject pair = env->CallObjectMethod(arrayList, java_ArrayList_get_, i);
    jlong low = env->CallLongMethod(pair, java_Pair_getStart_);
    jlong high  = env->CallLongMethod(pair, java_Pair_getEnd_);
    result.emplace_back(std::make_pair<uint64_t, uint64_t>((uint64_t)low,
                                                           (uint64_t)high));
    env->DeleteLocalRef(pair);
  }
  return result;
}

jobject to_java_map(JNIEnv *env, jobject obj, std::vector<genomic_field_t> genomic_fields) {
  jobject java_Map = env->NewObject(java_HashMap_, java_HashMap_init_);
  for (std::vector<genomic_field_t>::iterator it = genomic_fields.begin() ; it != genomic_fields.end(); ++it) {
    genomic_field_t field = *it;
    jstring key = env->NewStringUTF(field.first.c_str());
    jstring value = env->NewStringUTF(field.second.c_str());
    env->CallObjectMethod(java_Map, java_HashMap_put_, key, value);
    env->DeleteLocalRef(key);
    env->DeleteLocalRef(value);
  }
  return java_Map;
}

#if(0)
// TODO: Some version of this will be needed when we implement java bindings for GenomicsDB::query_variants()
jobject to_java_VariantCalls(JNIEnv *env, jclass cls, const char *array_name, GenomicsDB *genomicsdb, GenomicsDBVariantCalls variant_calls) {
  get_class_name(env, cls);
  auto java_VariantCalls =  env->NewObject(java_VariantCalls_, java_VariantCalls_init_default_);
  for (auto i=0ul; i<variant_calls.size(); i++) {
    const genomicsdb_variant_call_t* variant_call = variant_calls.at(i);
    genomic_interval_t genomic_interval = genomicsdb->get_genomic_interval(variant_call);
    jstring java_contigName = env->NewStringUTF(genomic_interval.contig_name.c_str());
    jobject java_genomicFields = to_java_map(env, cls, genomicsdb->get_genomic_fields(array_name, variant_call));
    env->CallObjectMethod(java_VariantCalls_, java_VariantCalls_init_, cls,
                   java_contigName,
                   genomic_interval.interval.first,
                   genomic_interval.interval.second,
                   java_genomicFields);
    env->DeleteLocalRef(java_contigName);
  }
  return java_VariantCalls;
}
#endif

void handleJNIException(JNIEnv *env, GenomicsDBException& exception) {
  jclass genomicsdb_java_exception_class = env->FindClass("org/genomicsdb/exception/GenomicsDBException");
  jboolean flag = env->ExceptionCheck();
  if (flag) {
    env->ExceptionClear();
  }
  std::string msg = std::string("JNI Error: ") + exception.what();
  env->ThrowNew(genomicsdb_java_exception_class, msg.c_str());
}

JNIEXPORT jstring JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniVersion(JNIEnv *env, jclass cls) {
  return env->NewStringUTF(genomicsdb_version().c_str());
}

JNIEXPORT jlong JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniConnect(JNIEnv *env,
                                                      jclass cls,
                                                      jstring workspace,
                                                      jstring vid_mapping_file,
                                                      jstring callset_mapping_file,
                                                      jstring reference_genome,
                                                      jobject attributes,
                                                      jlong segment_size) {
  // Convert
  auto workspace_cstr = env->GetStringUTFChars(workspace, NULL);
  auto vid_mapping_file_cstr = env->GetStringUTFChars(vid_mapping_file, NULL);
  auto callset_mapping_file_cstr = env->GetStringUTFChars(callset_mapping_file, NULL);
  auto reference_genome_cstr = env->GetStringUTFChars(reference_genome, NULL);

  GenomicsDB *genomicsdb = NULL;
  try {
    genomicsdb =  new GenomicsDB(workspace_cstr,
                                 callset_mapping_file_cstr,
                                 vid_mapping_file_cstr,
                                 reference_genome_cstr,
                                 to_string_vector(env, attributes),
                                 segment_size);
  } catch (GenomicsDBException& e) {
    handleJNIException(env, e);
  }
  
  // Cleanup
  env->ReleaseStringUTFChars(workspace, workspace_cstr);
  env->ReleaseStringUTFChars(vid_mapping_file, vid_mapping_file_cstr);
  env->ReleaseStringUTFChars(callset_mapping_file, callset_mapping_file_cstr);
  env->ReleaseStringUTFChars(reference_genome, reference_genome_cstr);
  
  return static_cast<jlong>(reinterpret_cast<uintptr_t>(genomicsdb));
}

JNIEXPORT void JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniDisconnect(JNIEnv *env,
                                                         jclass cls,
                                                         jlong handle) {
  delete reinterpret_cast<GenomicsDB *>(static_cast<std::uintptr_t>(handle));
}

class VariantCallProcessor : public GenomicsDBVariantCallProcessor {
 public:
  VariantCallProcessor(JNIEnv *env, jclass cls) {
    env_ = env;
    cls_ = cls;
    intervals_list_ = env->NewObject(java_ArrayList_, java_ArrayList_init_);
  }
  ~VariantCallProcessor() {
    finalize_interval();
    env_->DeleteLocalRef(intervals_list_);
  }
  jobject get_intervals_list() {
    if (intervals_list_) {
      return env_->NewGlobalRef(intervals_list_);
    } else {
      return NULL;
    }
  }
  void process(interval_t interval) {
    finalize_interval();
    current_calls_list_ =  env_->NewObject(java_VariantCalls_, java_VariantCalls_init_, (jlong)interval.first, (jlong)interval.second);
  }
  void process(uint32_t row, genomic_interval_t interval, std::vector<genomic_field_t> fields) {
    jstring java_contigName = env_->NewStringUTF(interval.contig_name.c_str());
    jobject java_fields = to_java_map(env_, cls_, fields);
    env_->CallObjectMethod(current_calls_list_, java_VariantCall_init_,
                          row,
                          java_contigName,
                          (jlong)interval.interval.first,
                          (jlong)interval.interval.second,
                          java_fields);
    env_->DeleteLocalRef(java_contigName);
    env_->DeleteLocalRef(java_fields);
  }
 private:
  void finalize_interval() {
    if (current_calls_list_) {
      env_->CallBooleanMethod(intervals_list_, java_ArrayList_add_, current_calls_list_);
      env_->DeleteLocalRef(current_calls_list_);
    }
    current_calls_list_ = NULL;
  }
  jobject current_calls_list_ = NULL;
  jobject intervals_list_ = NULL;

  JNIEnv *env_;
  jclass cls_;
};

JNIEXPORT jobject JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniQueryVariantCalls(JNIEnv *env,
                                                                jclass cls,
                                                                jlong handle,
                                                                jstring array_name,
                                                                jobject column_ranges,
                                                                jobject row_ranges) {
  // Convert
  GenomicsDB *genomicsdb = reinterpret_cast<GenomicsDB *>(static_cast<uintptr_t>(handle));
  auto array_name_cstr = env->GetStringUTFChars(array_name, NULL);

  VariantCallProcessor processor(env, cls);
  GenomicsDBVariantCalls variant_calls = genomicsdb->query_variant_calls(processor, array_name_cstr,
                                                                         to_genomicsdb_ranges_vector(env, column_ranges),
                                                                         to_genomicsdb_ranges_vector(env, row_ranges));

  // auto result = to_java_VariantCalls(env, cls, array_name_cstr, genomicsdb, variant_calls);
  
  env->ReleaseStringUTFChars(array_name, array_name_cstr);
  return processor.get_intervals_list();
}
