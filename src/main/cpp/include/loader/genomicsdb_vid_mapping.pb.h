// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: genomicsdb_vid_mapping.proto

#ifndef PROTOBUF_genomicsdb_5fvid_5fmapping_2eproto__INCLUDED
#define PROTOBUF_genomicsdb_5fvid_5fmapping_2eproto__INCLUDED

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 3000000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 3000002 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/metadata.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)

// Internal implementation detail -- do not call these.
void protobuf_AddDesc_genomicsdb_5fvid_5fmapping_2eproto();
void protobuf_AssignDesc_genomicsdb_5fvid_5fmapping_2eproto();
void protobuf_ShutdownFile_genomicsdb_5fvid_5fmapping_2eproto();

class Chromosome;
class InfoField;
class VidMapping;

// ===================================================================

class InfoField : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:InfoField) */ {
 public:
  InfoField();
  virtual ~InfoField();

  InfoField(const InfoField& from);

  inline InfoField& operator=(const InfoField& from) {
    CopyFrom(from);
    return *this;
  }

  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields();
  }

  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields();
  }

  static const ::google::protobuf::Descriptor* descriptor();
  static const InfoField& default_instance();

  void Swap(InfoField* other);

  // implements Message ----------------------------------------------

  inline InfoField* New() const { return New(NULL); }

  InfoField* New(::google::protobuf::Arena* arena) const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const InfoField& from);
  void MergeFrom(const InfoField& from);
  void Clear();
  bool IsInitialized() const;

  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const {
    return InternalSerializeWithCachedSizesToArray(false, output);
  }
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  void InternalSwap(InfoField* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return _internal_metadata_.arena();
  }
  inline void* MaybeArenaPtr() const {
    return _internal_metadata_.raw_arena_ptr();
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // required string name = 1;
  bool has_name() const;
  void clear_name();
  static const int kNameFieldNumber = 1;
  const ::std::string& name() const;
  void set_name(const ::std::string& value);
  void set_name(const char* value);
  void set_name(const char* value, size_t size);
  ::std::string* mutable_name();
  ::std::string* release_name();
  void set_allocated_name(::std::string* name);

  // required string type = 2;
  bool has_type() const;
  void clear_type();
  static const int kTypeFieldNumber = 2;
  const ::std::string& type() const;
  void set_type(const ::std::string& value);
  void set_type(const char* value);
  void set_type(const char* value, size_t size);
  ::std::string* mutable_type();
  ::std::string* release_type();
  void set_allocated_type(::std::string* type);

  // repeated string vcf_field_class = 3;
  int vcf_field_class_size() const;
  void clear_vcf_field_class();
  static const int kVcfFieldClassFieldNumber = 3;
  const ::std::string& vcf_field_class(int index) const;
  ::std::string* mutable_vcf_field_class(int index);
  void set_vcf_field_class(int index, const ::std::string& value);
  void set_vcf_field_class(int index, const char* value);
  void set_vcf_field_class(int index, const char* value, size_t size);
  ::std::string* add_vcf_field_class();
  void add_vcf_field_class(const ::std::string& value);
  void add_vcf_field_class(const char* value);
  void add_vcf_field_class(const char* value, size_t size);
  const ::google::protobuf::RepeatedPtrField< ::std::string>& vcf_field_class() const;
  ::google::protobuf::RepeatedPtrField< ::std::string>* mutable_vcf_field_class();

  // optional string length = 4;
  bool has_length() const;
  void clear_length();
  static const int kLengthFieldNumber = 4;
  const ::std::string& length() const;
  void set_length(const ::std::string& value);
  void set_length(const char* value);
  void set_length(const char* value, size_t size);
  ::std::string* mutable_length();
  ::std::string* release_length();
  void set_allocated_length(::std::string* length);

  // @@protoc_insertion_point(class_scope:InfoField)
 private:
  inline void set_has_name();
  inline void clear_has_name();
  inline void set_has_type();
  inline void clear_has_type();
  inline void set_has_length();
  inline void clear_has_length();

  // helper for ByteSize()
  int RequiredFieldsByteSizeFallback() const;

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::uint32 _has_bits_[1];
  mutable int _cached_size_;
  ::google::protobuf::internal::ArenaStringPtr name_;
  ::google::protobuf::internal::ArenaStringPtr type_;
  ::google::protobuf::RepeatedPtrField< ::std::string> vcf_field_class_;
  ::google::protobuf::internal::ArenaStringPtr length_;
  friend void  protobuf_AddDesc_genomicsdb_5fvid_5fmapping_2eproto();
  friend void protobuf_AssignDesc_genomicsdb_5fvid_5fmapping_2eproto();
  friend void protobuf_ShutdownFile_genomicsdb_5fvid_5fmapping_2eproto();

  void InitAsDefaultInstance();
  static InfoField* default_instance_;
};
// -------------------------------------------------------------------

class Chromosome : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:Chromosome) */ {
 public:
  Chromosome();
  virtual ~Chromosome();

  Chromosome(const Chromosome& from);

  inline Chromosome& operator=(const Chromosome& from) {
    CopyFrom(from);
    return *this;
  }

  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields();
  }

  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields();
  }

  static const ::google::protobuf::Descriptor* descriptor();
  static const Chromosome& default_instance();

  void Swap(Chromosome* other);

  // implements Message ----------------------------------------------

  inline Chromosome* New() const { return New(NULL); }

  Chromosome* New(::google::protobuf::Arena* arena) const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const Chromosome& from);
  void MergeFrom(const Chromosome& from);
  void Clear();
  bool IsInitialized() const;

  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const {
    return InternalSerializeWithCachedSizesToArray(false, output);
  }
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  void InternalSwap(Chromosome* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return _internal_metadata_.arena();
  }
  inline void* MaybeArenaPtr() const {
    return _internal_metadata_.raw_arena_ptr();
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // required string name = 1;
  bool has_name() const;
  void clear_name();
  static const int kNameFieldNumber = 1;
  const ::std::string& name() const;
  void set_name(const ::std::string& value);
  void set_name(const char* value);
  void set_name(const char* value, size_t size);
  ::std::string* mutable_name();
  ::std::string* release_name();
  void set_allocated_name(::std::string* name);

  // required int64 length = 2;
  bool has_length() const;
  void clear_length();
  static const int kLengthFieldNumber = 2;
  ::google::protobuf::int64 length() const;
  void set_length(::google::protobuf::int64 value);

  // required int64 tiledb_column_offset = 3;
  bool has_tiledb_column_offset() const;
  void clear_tiledb_column_offset();
  static const int kTiledbColumnOffsetFieldNumber = 3;
  ::google::protobuf::int64 tiledb_column_offset() const;
  void set_tiledb_column_offset(::google::protobuf::int64 value);

  // @@protoc_insertion_point(class_scope:Chromosome)
 private:
  inline void set_has_name();
  inline void clear_has_name();
  inline void set_has_length();
  inline void clear_has_length();
  inline void set_has_tiledb_column_offset();
  inline void clear_has_tiledb_column_offset();

  // helper for ByteSize()
  int RequiredFieldsByteSizeFallback() const;

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::uint32 _has_bits_[1];
  mutable int _cached_size_;
  ::google::protobuf::internal::ArenaStringPtr name_;
  ::google::protobuf::int64 length_;
  ::google::protobuf::int64 tiledb_column_offset_;
  friend void  protobuf_AddDesc_genomicsdb_5fvid_5fmapping_2eproto();
  friend void protobuf_AssignDesc_genomicsdb_5fvid_5fmapping_2eproto();
  friend void protobuf_ShutdownFile_genomicsdb_5fvid_5fmapping_2eproto();

  void InitAsDefaultInstance();
  static Chromosome* default_instance_;
};
// -------------------------------------------------------------------

class VidMapping : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:VidMapping) */ {
 public:
  VidMapping();
  virtual ~VidMapping();

  VidMapping(const VidMapping& from);

  inline VidMapping& operator=(const VidMapping& from) {
    CopyFrom(from);
    return *this;
  }

  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields();
  }

  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields();
  }

  static const ::google::protobuf::Descriptor* descriptor();
  static const VidMapping& default_instance();

  void Swap(VidMapping* other);

  // implements Message ----------------------------------------------

  inline VidMapping* New() const { return New(NULL); }

  VidMapping* New(::google::protobuf::Arena* arena) const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const VidMapping& from);
  void MergeFrom(const VidMapping& from);
  void Clear();
  bool IsInitialized() const;

  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const {
    return InternalSerializeWithCachedSizesToArray(false, output);
  }
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  void InternalSwap(VidMapping* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return _internal_metadata_.arena();
  }
  inline void* MaybeArenaPtr() const {
    return _internal_metadata_.raw_arena_ptr();
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // repeated .InfoField infofields = 1;
  int infofields_size() const;
  void clear_infofields();
  static const int kInfofieldsFieldNumber = 1;
  const ::InfoField& infofields(int index) const;
  ::InfoField* mutable_infofields(int index);
  ::InfoField* add_infofields();
  ::google::protobuf::RepeatedPtrField< ::InfoField >*
      mutable_infofields();
  const ::google::protobuf::RepeatedPtrField< ::InfoField >&
      infofields() const;

  // repeated .Chromosome chromosomes = 2;
  int chromosomes_size() const;
  void clear_chromosomes();
  static const int kChromosomesFieldNumber = 2;
  const ::Chromosome& chromosomes(int index) const;
  ::Chromosome* mutable_chromosomes(int index);
  ::Chromosome* add_chromosomes();
  ::google::protobuf::RepeatedPtrField< ::Chromosome >*
      mutable_chromosomes();
  const ::google::protobuf::RepeatedPtrField< ::Chromosome >&
      chromosomes() const;

  // @@protoc_insertion_point(class_scope:VidMapping)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::uint32 _has_bits_[1];
  mutable int _cached_size_;
  ::google::protobuf::RepeatedPtrField< ::InfoField > infofields_;
  ::google::protobuf::RepeatedPtrField< ::Chromosome > chromosomes_;
  friend void  protobuf_AddDesc_genomicsdb_5fvid_5fmapping_2eproto();
  friend void protobuf_AssignDesc_genomicsdb_5fvid_5fmapping_2eproto();
  friend void protobuf_ShutdownFile_genomicsdb_5fvid_5fmapping_2eproto();

  void InitAsDefaultInstance();
  static VidMapping* default_instance_;
};
// ===================================================================


// ===================================================================

#if !PROTOBUF_INLINE_NOT_IN_HEADERS
// InfoField

// required string name = 1;
inline bool InfoField::has_name() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void InfoField::set_has_name() {
  _has_bits_[0] |= 0x00000001u;
}
inline void InfoField::clear_has_name() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void InfoField::clear_name() {
  name_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
  clear_has_name();
}
inline const ::std::string& InfoField::name() const {
  // @@protoc_insertion_point(field_get:InfoField.name)
  return name_.GetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void InfoField::set_name(const ::std::string& value) {
  set_has_name();
  name_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:InfoField.name)
}
inline void InfoField::set_name(const char* value) {
  set_has_name();
  name_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:InfoField.name)
}
inline void InfoField::set_name(const char* value, size_t size) {
  set_has_name();
  name_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:InfoField.name)
}
inline ::std::string* InfoField::mutable_name() {
  set_has_name();
  // @@protoc_insertion_point(field_mutable:InfoField.name)
  return name_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* InfoField::release_name() {
  // @@protoc_insertion_point(field_release:InfoField.name)
  clear_has_name();
  return name_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void InfoField::set_allocated_name(::std::string* name) {
  if (name != NULL) {
    set_has_name();
  } else {
    clear_has_name();
  }
  name_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), name);
  // @@protoc_insertion_point(field_set_allocated:InfoField.name)
}

// required string type = 2;
inline bool InfoField::has_type() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void InfoField::set_has_type() {
  _has_bits_[0] |= 0x00000002u;
}
inline void InfoField::clear_has_type() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void InfoField::clear_type() {
  type_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
  clear_has_type();
}
inline const ::std::string& InfoField::type() const {
  // @@protoc_insertion_point(field_get:InfoField.type)
  return type_.GetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void InfoField::set_type(const ::std::string& value) {
  set_has_type();
  type_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:InfoField.type)
}
inline void InfoField::set_type(const char* value) {
  set_has_type();
  type_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:InfoField.type)
}
inline void InfoField::set_type(const char* value, size_t size) {
  set_has_type();
  type_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:InfoField.type)
}
inline ::std::string* InfoField::mutable_type() {
  set_has_type();
  // @@protoc_insertion_point(field_mutable:InfoField.type)
  return type_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* InfoField::release_type() {
  // @@protoc_insertion_point(field_release:InfoField.type)
  clear_has_type();
  return type_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void InfoField::set_allocated_type(::std::string* type) {
  if (type != NULL) {
    set_has_type();
  } else {
    clear_has_type();
  }
  type_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), type);
  // @@protoc_insertion_point(field_set_allocated:InfoField.type)
}

// repeated string vcf_field_class = 3;
inline int InfoField::vcf_field_class_size() const {
  return vcf_field_class_.size();
}
inline void InfoField::clear_vcf_field_class() {
  vcf_field_class_.Clear();
}
inline const ::std::string& InfoField::vcf_field_class(int index) const {
  // @@protoc_insertion_point(field_get:InfoField.vcf_field_class)
  return vcf_field_class_.Get(index);
}
inline ::std::string* InfoField::mutable_vcf_field_class(int index) {
  // @@protoc_insertion_point(field_mutable:InfoField.vcf_field_class)
  return vcf_field_class_.Mutable(index);
}
inline void InfoField::set_vcf_field_class(int index, const ::std::string& value) {
  // @@protoc_insertion_point(field_set:InfoField.vcf_field_class)
  vcf_field_class_.Mutable(index)->assign(value);
}
inline void InfoField::set_vcf_field_class(int index, const char* value) {
  vcf_field_class_.Mutable(index)->assign(value);
  // @@protoc_insertion_point(field_set_char:InfoField.vcf_field_class)
}
inline void InfoField::set_vcf_field_class(int index, const char* value, size_t size) {
  vcf_field_class_.Mutable(index)->assign(
    reinterpret_cast<const char*>(value), size);
  // @@protoc_insertion_point(field_set_pointer:InfoField.vcf_field_class)
}
inline ::std::string* InfoField::add_vcf_field_class() {
  // @@protoc_insertion_point(field_add_mutable:InfoField.vcf_field_class)
  return vcf_field_class_.Add();
}
inline void InfoField::add_vcf_field_class(const ::std::string& value) {
  vcf_field_class_.Add()->assign(value);
  // @@protoc_insertion_point(field_add:InfoField.vcf_field_class)
}
inline void InfoField::add_vcf_field_class(const char* value) {
  vcf_field_class_.Add()->assign(value);
  // @@protoc_insertion_point(field_add_char:InfoField.vcf_field_class)
}
inline void InfoField::add_vcf_field_class(const char* value, size_t size) {
  vcf_field_class_.Add()->assign(reinterpret_cast<const char*>(value), size);
  // @@protoc_insertion_point(field_add_pointer:InfoField.vcf_field_class)
}
inline const ::google::protobuf::RepeatedPtrField< ::std::string>&
InfoField::vcf_field_class() const {
  // @@protoc_insertion_point(field_list:InfoField.vcf_field_class)
  return vcf_field_class_;
}
inline ::google::protobuf::RepeatedPtrField< ::std::string>*
InfoField::mutable_vcf_field_class() {
  // @@protoc_insertion_point(field_mutable_list:InfoField.vcf_field_class)
  return &vcf_field_class_;
}

// optional string length = 4;
inline bool InfoField::has_length() const {
  return (_has_bits_[0] & 0x00000008u) != 0;
}
inline void InfoField::set_has_length() {
  _has_bits_[0] |= 0x00000008u;
}
inline void InfoField::clear_has_length() {
  _has_bits_[0] &= ~0x00000008u;
}
inline void InfoField::clear_length() {
  length_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
  clear_has_length();
}
inline const ::std::string& InfoField::length() const {
  // @@protoc_insertion_point(field_get:InfoField.length)
  return length_.GetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void InfoField::set_length(const ::std::string& value) {
  set_has_length();
  length_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:InfoField.length)
}
inline void InfoField::set_length(const char* value) {
  set_has_length();
  length_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:InfoField.length)
}
inline void InfoField::set_length(const char* value, size_t size) {
  set_has_length();
  length_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:InfoField.length)
}
inline ::std::string* InfoField::mutable_length() {
  set_has_length();
  // @@protoc_insertion_point(field_mutable:InfoField.length)
  return length_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* InfoField::release_length() {
  // @@protoc_insertion_point(field_release:InfoField.length)
  clear_has_length();
  return length_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void InfoField::set_allocated_length(::std::string* length) {
  if (length != NULL) {
    set_has_length();
  } else {
    clear_has_length();
  }
  length_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), length);
  // @@protoc_insertion_point(field_set_allocated:InfoField.length)
}

// -------------------------------------------------------------------

// Chromosome

// required string name = 1;
inline bool Chromosome::has_name() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void Chromosome::set_has_name() {
  _has_bits_[0] |= 0x00000001u;
}
inline void Chromosome::clear_has_name() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void Chromosome::clear_name() {
  name_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
  clear_has_name();
}
inline const ::std::string& Chromosome::name() const {
  // @@protoc_insertion_point(field_get:Chromosome.name)
  return name_.GetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void Chromosome::set_name(const ::std::string& value) {
  set_has_name();
  name_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:Chromosome.name)
}
inline void Chromosome::set_name(const char* value) {
  set_has_name();
  name_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:Chromosome.name)
}
inline void Chromosome::set_name(const char* value, size_t size) {
  set_has_name();
  name_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:Chromosome.name)
}
inline ::std::string* Chromosome::mutable_name() {
  set_has_name();
  // @@protoc_insertion_point(field_mutable:Chromosome.name)
  return name_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* Chromosome::release_name() {
  // @@protoc_insertion_point(field_release:Chromosome.name)
  clear_has_name();
  return name_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void Chromosome::set_allocated_name(::std::string* name) {
  if (name != NULL) {
    set_has_name();
  } else {
    clear_has_name();
  }
  name_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), name);
  // @@protoc_insertion_point(field_set_allocated:Chromosome.name)
}

// required int64 length = 2;
inline bool Chromosome::has_length() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void Chromosome::set_has_length() {
  _has_bits_[0] |= 0x00000002u;
}
inline void Chromosome::clear_has_length() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void Chromosome::clear_length() {
  length_ = GOOGLE_LONGLONG(0);
  clear_has_length();
}
inline ::google::protobuf::int64 Chromosome::length() const {
  // @@protoc_insertion_point(field_get:Chromosome.length)
  return length_;
}
inline void Chromosome::set_length(::google::protobuf::int64 value) {
  set_has_length();
  length_ = value;
  // @@protoc_insertion_point(field_set:Chromosome.length)
}

// required int64 tiledb_column_offset = 3;
inline bool Chromosome::has_tiledb_column_offset() const {
  return (_has_bits_[0] & 0x00000004u) != 0;
}
inline void Chromosome::set_has_tiledb_column_offset() {
  _has_bits_[0] |= 0x00000004u;
}
inline void Chromosome::clear_has_tiledb_column_offset() {
  _has_bits_[0] &= ~0x00000004u;
}
inline void Chromosome::clear_tiledb_column_offset() {
  tiledb_column_offset_ = GOOGLE_LONGLONG(0);
  clear_has_tiledb_column_offset();
}
inline ::google::protobuf::int64 Chromosome::tiledb_column_offset() const {
  // @@protoc_insertion_point(field_get:Chromosome.tiledb_column_offset)
  return tiledb_column_offset_;
}
inline void Chromosome::set_tiledb_column_offset(::google::protobuf::int64 value) {
  set_has_tiledb_column_offset();
  tiledb_column_offset_ = value;
  // @@protoc_insertion_point(field_set:Chromosome.tiledb_column_offset)
}

// -------------------------------------------------------------------

// VidMapping

// repeated .InfoField infofields = 1;
inline int VidMapping::infofields_size() const {
  return infofields_.size();
}
inline void VidMapping::clear_infofields() {
  infofields_.Clear();
}
inline const ::InfoField& VidMapping::infofields(int index) const {
  // @@protoc_insertion_point(field_get:VidMapping.infofields)
  return infofields_.Get(index);
}
inline ::InfoField* VidMapping::mutable_infofields(int index) {
  // @@protoc_insertion_point(field_mutable:VidMapping.infofields)
  return infofields_.Mutable(index);
}
inline ::InfoField* VidMapping::add_infofields() {
  // @@protoc_insertion_point(field_add:VidMapping.infofields)
  return infofields_.Add();
}
inline ::google::protobuf::RepeatedPtrField< ::InfoField >*
VidMapping::mutable_infofields() {
  // @@protoc_insertion_point(field_mutable_list:VidMapping.infofields)
  return &infofields_;
}
inline const ::google::protobuf::RepeatedPtrField< ::InfoField >&
VidMapping::infofields() const {
  // @@protoc_insertion_point(field_list:VidMapping.infofields)
  return infofields_;
}

// repeated .Chromosome chromosomes = 2;
inline int VidMapping::chromosomes_size() const {
  return chromosomes_.size();
}
inline void VidMapping::clear_chromosomes() {
  chromosomes_.Clear();
}
inline const ::Chromosome& VidMapping::chromosomes(int index) const {
  // @@protoc_insertion_point(field_get:VidMapping.chromosomes)
  return chromosomes_.Get(index);
}
inline ::Chromosome* VidMapping::mutable_chromosomes(int index) {
  // @@protoc_insertion_point(field_mutable:VidMapping.chromosomes)
  return chromosomes_.Mutable(index);
}
inline ::Chromosome* VidMapping::add_chromosomes() {
  // @@protoc_insertion_point(field_add:VidMapping.chromosomes)
  return chromosomes_.Add();
}
inline ::google::protobuf::RepeatedPtrField< ::Chromosome >*
VidMapping::mutable_chromosomes() {
  // @@protoc_insertion_point(field_mutable_list:VidMapping.chromosomes)
  return &chromosomes_;
}
inline const ::google::protobuf::RepeatedPtrField< ::Chromosome >&
VidMapping::chromosomes() const {
  // @@protoc_insertion_point(field_list:VidMapping.chromosomes)
  return chromosomes_;
}

#endif  // !PROTOBUF_INLINE_NOT_IN_HEADERS
// -------------------------------------------------------------------

// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_genomicsdb_5fvid_5fmapping_2eproto__INCLUDED
