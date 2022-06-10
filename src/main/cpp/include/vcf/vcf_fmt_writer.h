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
      : m_ptr(&m_buffer),
      m_size_threshold(1024u*1024u) {
    }
    //VCFWriterNoOverflow(std::ostream&)
    template<typename WriteTargetTy_dummy=WriteTargetTy>
    VCFWriterNoOverflow(typename std::enable_if<std::is_base_of<std::ostream, WriteTargetTy_dummy>::value, std::ostream>::type& ref)
      : m_ptr(&ref),
      m_size_threshold(1024u*1024u) {
    }

    template<class T, bool assume_valid=false>
    inline bool write(const T v) {
      //if called with assume_valid=true, I hope the compiler gets rid of this if stmt
      if(!assume_valid && !is_bcf_valid_value<T>(v))
        return FmtWriter::write<char>(m_buffer, '.');
      else
        return FmtWriter::write<T>(m_buffer, v);
    };

    //Returns true if vector end
    template<class T>
    bool check_for_bcf_invalid_values_and_write(const T val) {
      if(is_bcf_vector_end_value<T>(val))
        return true;
      if(is_bcf_missing_value<T>(val))
        FmtWriter::write<char>(m_buffer, '.');
      else
        FmtWriter::write<T>(m_buffer, val);
      return false;
    }

    //Returns true if vector end
    template<class T>
    bool check_for_bcf_invalid_values_and_write(const T val, const char sep) {
      if(is_bcf_vector_end_value<T>(val))
        return true;
      FmtWriter::write<char>(m_buffer, sep);
      if(is_bcf_missing_value<T>(val))
        FmtWriter::write<char>(m_buffer, '.');
      else
        FmtWriter::write<T>(m_buffer, val);
      return false;
    }

    template<class T, bool assume_valid=false>
    inline bool write(const T* v, const size_t n) {
      if(assume_valid)
        return FmtWriter::write<T>(m_buffer, v, n);
      else {
        for(auto i=0u;i<n;++i) {
          if(check_for_bcf_invalid_values_and_write<T>(v[i]))
            return true;
        }
        return true;
      }
    }

    template<class T, bool assume_valid=false>
    inline bool write(const T* v, const size_t n, const char sep) {
      if(assume_valid)
        return FmtWriter::write<T>(m_buffer, v, n, sep);
      else {
        if(n > 0u) {
          if(check_for_bcf_invalid_values_and_write<T>(v[0u]))
            return true;
          for(auto i=1u;i<n;++i) {
            if(check_for_bcf_invalid_values_and_write<T>(v[i], sep))
              return true;
          }
        }
        return true;
      }
    }

    inline bool write_GT_allele_index(const uint64_t row_query_idx, const int v) {
      if(v == get_bcf_gt_no_call_allele_index<int>())
        FmtWriter::write<char>(m_buffer, '.');
      else
        FmtWriter::write<int>(m_buffer, v);
      return true;
    }
 
    inline bool write_GT_phase(const uint64_t row_query_idx, const int v) {
      FmtWriter::write<char>(m_buffer, (v == 0) ? '/' : '|');
      return true;
    }

    inline bool write_GT_empty(const uint64_t row_query_idx) {
      FmtWriter::write<char>(m_buffer, '.');
      return true;
    }

    template <typename element_type>
    inline bool write_empty(unsigned field_query_idx, const uint64_t row_query_idx) {
      return true;
    }

    template <typename element_type, bool is_first_element>
    inline bool write_missing(unsigned field_query_idx, const uint64_t row_query_idx) {
      if (!is_first_element) FmtWriter::write<char>(m_buffer, ',');
      FmtWriter::write<char>(m_buffer, '.');
      return true;
    }

    template <typename element_type, bool is_first_element>
    inline bool write_element(unsigned field_query_idx, const uint64_t row_query_idx, const element_type value) {
      if (!is_first_element) FmtWriter::write<char>(m_buffer, ',');
      FmtWriter::write<element_type>(m_buffer, value);
      return true;
    }

    void reset(const size_t n) { m_buffer.resize(n); }

    std::string& get_string_buffer() {
      return m_buffer;
    }

    //When T is std::string, this function is a nop - caller must explicitly decide when to clear the string etc
    template<typename WriteTargetTy_dummy=WriteTargetTy, typename std::enable_if<std::is_same<WriteTargetTy_dummy, std::string>::value, bool>::type x= false>
    void flush_buffer_if_large_and_if_flush_supported() {
    }

    //When T is ostream, write contents of buffer
    template<typename WriteTargetTy_dummy=WriteTargetTy, typename std::enable_if<std::is_base_of<std::ostream, WriteTargetTy_dummy>::value, bool>::type x= false>
    void flush_buffer_if_large_and_if_flush_supported() {
      if(m_buffer.size() > m_size_threshold) {
        m_ptr->write(m_buffer.data(), m_buffer.size());
        m_buffer.clear();
      }
    }

    //When T is std::string, this function is a nop - caller must explicitly decide when to clear the string etc
    template<typename WriteTargetTy_dummy=WriteTargetTy, typename std::enable_if<std::is_same<WriteTargetTy_dummy, std::string>::value, bool>::type x= false>
    void flush_buffer_if_flush_supported() {
    }

    //When T is ostream, write contents of buffer
    template<typename WriteTargetTy_dummy=WriteTargetTy, typename std::enable_if<std::is_base_of<std::ostream, WriteTargetTy_dummy>::value, bool>::type x= false>
    void flush_buffer_if_flush_supported() {
      m_ptr->write(m_buffer.data(), m_buffer.size());
      m_buffer.clear();
      m_ptr->flush();
    }
  private:
    WriteTargetTy* m_ptr; 
    size_t m_size_threshold;
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

    //Returns std::pair(no_overflow, is_vector_end)
    template<class T>
    std::pair<bool,bool> check_for_bcf_invalid_values_and_write(const T val) {
      if(is_bcf_vector_end_value<T>(val))
        return std::pair<bool, bool>(true, true);
      if(is_bcf_missing_value<T>(val))
        return std::pair<bool, bool>(FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.'), false);
      else
        return std::pair<bool, bool>(FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, val), false);
    }

    //Returns std::pair(no_overflow, is_vector_end)
    template<class T>
    std::pair<bool,bool> check_for_bcf_invalid_values_and_write(const T val, const char sep) {
      if(is_bcf_vector_end_value<T>(val))
        return std::pair<bool, bool>(true, true);
      if(!FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, sep))
        return std::pair<bool, bool>(false, false);
      if(is_bcf_missing_value<T>(val))
        return std::pair<bool, bool>(FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.'), false);
      else
        return std::pair<bool, bool>(FmtWriter::write_if_space_available<T>(m_ptr, m_size, m_offset, val), false);
    }


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
          auto no_overflow_vector_end_pair = check_for_bcf_invalid_values_and_write<T>(v[i]);
          no_overflow = no_overflow_vector_end_pair.first;
          if(!no_overflow || no_overflow_vector_end_pair.second)
            break;
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
        if(n > 0u) {
          auto no_overflow_vector_end_pair = check_for_bcf_invalid_values_and_write<T>(v[0u]);
          if(!no_overflow_vector_end_pair.first)
            return false;
          if(no_overflow_vector_end_pair.second) //vector end, no overflow
            return true;
          for(auto i=1u;i<n;++i) {
            auto no_overflow_vector_end_pair = check_for_bcf_invalid_values_and_write<T>(v[i], sep);
            if(!no_overflow_vector_end_pair.first)
              return false;
            if(no_overflow_vector_end_pair.second) //vector end, no overflow
              return true;
          }
        }
        return true;
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

    template <typename element_type>
    inline bool write_empty(unsigned field_query_idx, const uint64_t row_query_idx) {
      return true;
    }

    template <typename element_type, bool is_first_element>
    inline bool write_missing(unsigned field_query_idx, const uint64_t row_query_idx) {
      auto no_overflow = true;
      if (!is_first_element) no_overflow = FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, ',');
      no_overflow = no_overflow && FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, '.');
      return no_overflow;
    }

    template <typename element_type, bool is_first_element>
    inline bool write_element(unsigned field_query_idx, const uint64_t row_query_idx, const element_type value) {
      auto no_overflow = true;
      if (!is_first_element) no_overflow = FmtWriter::write_if_space_available<char>(m_ptr, m_size, m_offset, ',');
      no_overflow = no_overflow && FmtWriter::write_if_space_available<element_type>(m_ptr, m_size, m_offset, value);
      return no_overflow;
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
