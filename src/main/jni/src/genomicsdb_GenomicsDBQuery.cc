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

#include <assert.h>
#include <mutex>

//java.util.ArrayList
static jclass java_ArrayList_ ;
static jmethodID java_ArrayList_init_;
static jmethodID java_ArrayList_size_;
static jmethodID java_ArrayList_get_;
static jmethodID java_ArrayList_add_;

//java.util.Hashmap
static jclass java_HashMap_;
static jmethodID java_HashMap_init_;
static jmethodID java_HashMap_put_;

//org.genomicsdb.reader.GenomicsDBQuery$Interval
static jclass java_Interval_;
static jmethodID java_Interval_init_default_;
static jmethodID java_Interval_init_;
static jmethodID java_Interval_addCall_;

//org.genomicsdb.reader.GenomicsDBQuery$VariantCall
static jclass java_VariantCall_;
static jmethodID java_VariantCall_init_;

//org.genomicsdb.reader.GenomicsDBQuery$Pair
static jclass java_Pair_;
static jmethodID java_Pair_init_;
static jmethodID java_Pair_getStart_;
static jmethodID java_Pair_getEnd_;

#define INIT(VAR,CODE)                                                                                 \
  do {                                                                                                 \
    VAR = CODE;                                                                                        \
    if (VAR == NULL) {                                                                                 \
      throw GenomicsDBException("genomicsdb_GenomicsDBQuery.cc#JNI_OnLoad:"+std::to_string(__LINE__)); \
    }                                                                                                  \
  } while (false)

#define GENOMICSDB_JNI_VERSION JNI_VERSION_1_8

JNIEXPORT void JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniInitialize(JNIEnv *env, jclass cls) {
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
    
    //org.genomicsdb.reader.GenomicsDBQuery$Interval
    INIT(java_Interval_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("org/genomicsdb/reader/GenomicsDBQuery$Interval"))));
    INIT(java_Interval_init_default_,  env->GetMethodID(java_Interval_, "<init>", "()V"));
    INIT(java_Interval_init_,  env->GetMethodID(java_Interval_, "<init>", "(JJ)V"));
    INIT(java_Interval_addCall_, env->GetMethodID(java_Interval_, "addCall", "(Lorg/genomicsdb/reader/GenomicsDBQuery$VariantCall;)V"));

    //org.genomicsdb.reader.GenomicDBQuery$VariantCall
    INIT(java_VariantCall_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("org/genomicsdb/reader/GenomicsDBQuery$VariantCall"))));
    INIT(java_VariantCall_init_, env->GetMethodID(java_VariantCall_, "<init>", "(JJLjava/lang/String;Ljava/lang/String;JJLjava/util/Map;)V"));

    //org.genomicsdb.reader.GenomicsDBQuery$Pair<L, R>
    INIT(java_Pair_, static_cast<jclass>(env->NewGlobalRef(env->FindClass("org/genomicsdb/reader/GenomicsDBQuery$Pair"))));
    INIT(java_Pair_init_, env->GetMethodID(java_Pair_, "<init>", "(JJ)V"));
    INIT(java_Pair_getStart_, env->GetMethodID(java_Pair_, "getStart", "()J"));
    INIT(java_Pair_getEnd_, env->GetMethodID(java_Pair_, "getEnd", "()J"));
}

void JNI_OnUnload(JavaVM *vm, void *reserved) {
  JNIEnv *env;
  if (vm->GetEnv(reinterpret_cast<void**>(&env), GENOMICSDB_JNI_VERSION) == JNI_OK) {
    env->DeleteGlobalRef(java_Pair_);
    env->DeleteGlobalRef(java_VariantCall_);
    env->DeleteGlobalRef(java_Interval_);
    env->DeleteGlobalRef(java_HashMap_);
    env->DeleteGlobalRef(java_ArrayList_);
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
  return result;
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

jobject to_java_map(JNIEnv *env, jobject obj, std::vector<genomic_field_t> genomic_fields, const std::shared_ptr<std::map<std::string, genomic_field_type_t>> genomic_field_types) {
  jobject java_Map = env->NewObject(java_HashMap_, java_HashMap_init_);
  for (std::vector<genomic_field_t>::iterator it = genomic_fields.begin() ; it != genomic_fields.end(); ++it) {
    genomic_field_t field = *it;
    jstring key = env->NewStringUTF(field.name.c_str());
    if (genomic_field_types->find(field.name) == genomic_field_types->end()) {
      throw GenomicsDBException("Genomic Field="+field.name+" does not seem to have an associated type");
    }
    // TODO: Return appropriate objects based on type
    const genomic_field_type_t field_type =  genomic_field_types->at(field.name);
    jstring value = env->NewStringUTF(field.to_string(field_type).c_str());
    env->CallObjectMethod(java_Map, java_HashMap_put_, key, value);
    env->DeleteLocalRef(key);
    env->DeleteLocalRef(value);
  }
  return java_Map;
}

#if(0)
// TODO: Some version of this will be needed when we implement java bindings for GenomicsDB::query_variants()
jobject to_java_Interval(JNIEnv *env, jclass cls, const char *array_name, GenomicsDB *genomicsdb, GenomicsDBVariantCalls variant_calls) {
  get_class_name(env, cls);
  auto java_Interval =  env->NewObject(java_Interval_, java_Interval_init_default_);
  for (auto i=0ul; i<variant_calls.size(); i++) {
    const genomicsdb_variant_call_t* variant_call = variant_calls.at(i);
    genomic_interval_t genomic_interval = genomicsdb->get_genomic_interval(variant_call);
    jstring java_contigName = env->NewStringUTF(genomic_interval.contig_name.c_str());
    jobject java_genomicFields = to_java_map(env, cls, genomicsdb->get_genomic_fields(array_name, variant_call));
    env->CallObjectMethod(java_Interval_, java_Interval_init_, cls,
                   java_contigName,
                   genomic_interval.interval.first,
                   genomic_interval.interval.second,
                   java_genomicFields);
    env->DeleteLocalRef(java_contigName);
  }
  return java_Interval;
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


JNIEXPORT jlong JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniConnectJSON(JNIEnv *env,
                                                      jclass cls,
                                                      jstring query_json_file,
                                                      jstring loader_json_file) {
  // Convert
  auto query_json_file_cstr = env->GetStringUTFChars(query_json_file, NULL);
  auto loader_json_file_cstr = env->GetStringUTFChars(loader_json_file, NULL);

  GenomicsDB *genomicsdb = NULL;
  try {
    genomicsdb =  new GenomicsDB(query_json_file_cstr,
                                 loader_json_file_cstr);
  } catch (GenomicsDBException& e) {
    handleJNIException(env, e);
  }
  
  // Cleanup
  env->ReleaseStringUTFChars(query_json_file, query_json_file_cstr);
  env->ReleaseStringUTFChars(loader_json_file, loader_json_file_cstr);

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
  }
  jobject get_intervals_list() {
    return intervals_list_;
  }
  void process(const interval_t& interval) {
    finalize_interval();
    current_calls_list_ =  env_->NewObject(java_Interval_, java_Interval_init_, (jlong)interval.first, (jlong)interval.second);
  }
  void process(const std::string& sample_name,
               const int64_t* coordinates,
               const genomic_interval_t& interval,
               const std::vector<genomic_field_t>& fields) {
    int64_t row = coordinates[0];
    int64_t col = coordinates[1];
    jstring java_sample_name = env_->NewStringUTF(sample_name.c_str());
    jstring java_contig_name = env_->NewStringUTF(interval.contig_name.c_str());
    jobject java_fields = to_java_map(env_, cls_, fields, get_genomic_field_types());
    jobject java_variant_call = env_->NewObject(java_VariantCall_, java_VariantCall_init_,
                                                (jlong)row, (jlong)col,
                                                java_sample_name,
                                                java_contig_name,
                                                (jlong)interval.interval.first,
                                                (jlong)interval.interval.second,
                                                java_fields);
    env_->DeleteLocalRef(java_sample_name);
    env_->DeleteLocalRef(java_contig_name);
    env_->DeleteLocalRef(java_fields);
    if (java_variant_call) {
      assert(current_calls_list_);
      env_->CallObjectMethod(current_calls_list_, java_Interval_addCall_, java_variant_call);
      env_->DeleteLocalRef(java_variant_call);
    } else {
      throw GenomicsDBException("Could not create Java VariantCall via JNI: "+std::to_string(__LINE__));
    }
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

  if (variant_calls.size() > 0) {
    // auto result = to_java_Interval(env, cls, array_name_cstr, genomicsdb, variant_calls);
    throw GenomicsDBException("NYI: processing results of genomicsdb_GenomicsDBQuery.cc#jniQueryInterval :"+std::to_string(__LINE__));
  }
  
  env->ReleaseStringUTFChars(array_name, array_name_cstr);
  return processor.get_intervals_list();
}

JNIEXPORT void JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniGenerateVCF(JNIEnv *env,
                                                          jclass cls,
                                                          jlong handle,
                                                          jstring array_name,
                                                          jobject column_ranges,
                                                          jobject row_ranges,
                                                          jstring output,
                                                          jstring output_format,
                                                          jboolean overwrite) {
  // Convert
  GenomicsDB *genomicsdb = reinterpret_cast<GenomicsDB *>(static_cast<uintptr_t>(handle));
  auto array_name_cstr = env->GetStringUTFChars(array_name, NULL);
  auto output_cstr =  env->GetStringUTFChars(output, NULL);
  auto output_format_cstr = env->GetStringUTFChars(output_format, NULL);

  genomicsdb->generate_vcf(array_name_cstr,
                           to_genomicsdb_ranges_vector(env, column_ranges),
                           to_genomicsdb_ranges_vector(env, row_ranges),
                           output_cstr, output_format_cstr, overwrite);

  // Cleanup
  env->ReleaseStringUTFChars(array_name, array_name_cstr);
  env->ReleaseStringUTFChars(output_format, output_format_cstr);
  env->ReleaseStringUTFChars(output, output_cstr);
}

JNIEXPORT void JNICALL
Java_org_genomicsdb_reader_GenomicsDBQuery_jniGenerateVCF1(JNIEnv *env,
                                                          jclass cls,
                                                          jlong handle,
                                                          jstring output,
                                                          jstring output_format,
                                                          jboolean overwrite) {
  // Convert
  GenomicsDB *genomicsdb = reinterpret_cast<GenomicsDB *>(static_cast<uintptr_t>(handle));
  auto output_cstr =  env->GetStringUTFChars(output, NULL);
  auto output_format_cstr = env->GetStringUTFChars(output_format, NULL);

  genomicsdb->generate_vcf(output_cstr, output_format_cstr, overwrite);

  // Cleanup
  env->ReleaseStringUTFChars(output_format, output_format_cstr);
  env->ReleaseStringUTFChars(output, output_cstr);
 }
