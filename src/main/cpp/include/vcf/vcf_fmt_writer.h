#ifndef VCF_FMT_WRITER_H
#define VCF_FMT_WRITER_H

#include<string>
#include<iostream>
#include<fstream>
#include "fmt_writer.h"
#include "vcf.h"

template<class WriteTargetTy>   //can be std::string or std::ostream
class VCFWriterNoOverflow {
  public:
    //Why do we need to add the dummy template parameter WriteTargetTy_dummy?
    //See https://stackoverflow.com/a/17842695
    //Only one of these constructors will be valid depending on WriteTargetTy
    //VCFWriterNoOverflow() when WriteTargetTy is std::string
    template<typename WriteTargetTy_dummy=WriteTargetTy>
    VCFWriterNoOverflow(typename std::enable_if<std::is_same<WriteTargetTy_dummy, std::string>::value>::type* =nullptr)
      : m_ptr(&m_buffer) { }
    //VCFWriterNoOverflow(std::ostream&)
    template<typename WriteTargetTy_dummy=WriteTargetTy>
    VCFWriterNoOverflow(typename std::enable_if<std::is_base_of<std::ostream, WriteTargetTy_dummy>::value, std::ostream>::type& ref)
      : m_ptr(&ref) { }

    template<class T, bool assume_valid=false>
    inline bool write(const T v) {
      //if called with assume_valid=true, I hope the compiler gets rid of this if stmt
      if(!assume_valid && !is_bcf_valid_value<T>(v))
        return FmtWriter::write<char>(*m_ptr, '.');
      else
        return FmtWriter::write<T>(*m_ptr, v);
    };
    
    template<class T, bool assume_valid=false>
    inline bool write(const T* v, const size_t n) {
      if(assume_valid)
        return FmtWriter::write<T>(*m_ptr, v, n);
      else {
        for(auto i=0u;i<n;++i) {
          const auto val = v[i];
          if(is_bcf_valid_value<T>(val))
            FmtWriter::write<T>(*m_ptr, val);
          else {
            if(is_bcf_vector_end_value<T>(val))
              break;
            FmtWriter::write<char>(*m_ptr, '.');
          }
        }
        return true;
      }
    }

    template<class T, bool assume_valid=false>
    inline bool write(const T* v, const size_t n, const char sep) {
      if(assume_valid)
        return FmtWriter::write<T>(*m_ptr, v, n, sep);
      else {
        if(n > 0u) {
          const auto val = v[0u];
          if(is_bcf_valid_value<T>(val))
            FmtWriter::write<T>(*m_ptr, val);
          else {
            if(is_bcf_vector_end_value<T>(val))
              return true;
            FmtWriter::write<char>(*m_ptr, '.');
          }
        }
        for(auto i=1u;i<n;++i) {
          FmtWriter::write<char>(*m_ptr, sep);
          const auto val = v[i];
          if(is_bcf_valid_value<T>(val))
            FmtWriter::write<T>(*m_ptr, val);
          else {
            if(is_bcf_vector_end_value<T>(val))
              break;
            FmtWriter::write<char>(*m_ptr, '.');
          }
        }
        return true;
      }
    }

    inline bool write_GT_allele_index(const uint64_t row_query_idx, const int v) {
      if(v == get_bcf_gt_no_call_allele_index<int>())
        FmtWriter::write<char>(*m_ptr, '.');
      else
        FmtWriter::write<int>(*m_ptr, v);
      return true;
    }
 
    inline bool write_GT_phase(const uint64_t row_query_idx, const int v) {
      FmtWriter::write<char>(*m_ptr, (v == 0) ? '/' : '|');
      return true;
    }

    inline bool write_GT_empty(const uint64_t row_query_idx) {
      FmtWriter::write<char>(*m_ptr, '.');
      return true;
    }

    void reset(const size_t n) { m_buffer.resize(n); }

    std::string& get_string_buffer() {
      return m_buffer;
    }
  private:
    WriteTargetTy* m_ptr; 
    std::string m_buffer;
};

//FSB - fixed size buffer
class VCFWriterFSB {
  public:
    VCFWriterFSB(char* buffer, const size_t max_size, const size_t offset=0u) {
      m_ptr = buffer;
      m_size = max_size;
      m_offset = offset;
    }
    inline size_t get_offset() const { return m_offset; }
    void reset(const size_t n) { m_offset = n; }

    template<class T, bool assume_valid=false>
    inline bool write(const T v) {
      //if called assume_valid=true, I hope the compiler gets rid of this if stmt
      if(!assume_valid && !is_bcf_valid_value<T>(v))
        return FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
      else
        return FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, v);
    }

    template<class T, bool assume_valid=false>
    inline bool write(const T* v, const size_t n) {
      //if called assume_valid=true, I hope the compiler gets rid of this if stmt
      if(assume_valid)
        return FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, v, n);
      else {
        auto no_overflow = true;
        for(auto i=0u;i<n;++i) {
          const auto val = v[i];
          if(is_bcf_valid_value<T>(val))
            no_overflow = no_overflow && FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, val);
          else {
            if(is_bcf_vector_end_value<T>(val))
              break;
            no_overflow = no_overflow && FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
          }
        }
        return no_overflow;
      }
    }

    template<class T, bool assume_valid=false>
    inline bool write(const T* v, const size_t n, const char sep) {
      //if called assume_valid=true, I hope the compiler gets rid of this if stmt
      if(assume_valid)
        return FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, v, n, sep);
      else {
        auto no_overflow = true;
        if(n > 0u) {
          const auto val = v[0u];
          if(is_bcf_valid_value<T>(val))
            no_overflow = no_overflow && FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, val);
          else {
            if(is_bcf_vector_end_value<T>(val))
              return no_overflow;
            no_overflow = no_overflow && FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
          }
        }
        for(auto i=1u;i<n;++i) {
          no_overflow = no_overflow && FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, sep);
          const auto val = v[i];
          if(is_bcf_valid_value<T>(val))
            no_overflow = no_overflow && FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, val);
          else {
            if(is_bcf_vector_end_value<T>(val))
              return no_overflow;
            no_overflow = no_overflow && FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
          }
        }
        return no_overflow;
      }
    }

    inline bool write_GT_allele_index(const uint64_t row_query_idx, const int v) {
      if(v == get_bcf_gt_no_call_allele_index<int>())
        return FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
      else
        return FmtWriter::write_if_space_available<int>(m_ptr, m_size, m_offset, v);
    }
 
    inline bool write_GT_phase(const uint64_t row_query_idx, const int v) {
      return FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset,
          (v == 0) ? '/' : '|');
    }

    inline bool write_GT_empty(const uint64_t row_query_idx) {
      return FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
    }

  private:
    char* m_ptr;
    size_t m_size;
    size_t m_offset;
};

enum VCFWRITER_ENUM {
  STL_STRING_NO_LIMIT = 0u,  
  STL_OSTREAM_NO_LIMIT,
  FIXED_SIZE_BUFFER
};

#endif
