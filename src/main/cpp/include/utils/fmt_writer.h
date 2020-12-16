/**
 * The MIT License (MIT)
 * Copyright (c) 2020 Omics Data Automation Inc
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

#ifndef FMT_LIB_WRITER_H
#define FMT_LIB_WRITER_H

#include <iterator>
#include <fmt/format.h>

//Code for fixed iterator copied from https://github.com/fmtlib/fmt/issues/764#issuecomment-395753853
//Useful when the buffer into which fmt must write is pre-allocated by the caller and its size cannot be changed
template <typename OutputIt>
class fmt_fixed_size_buffer_iterator {
  private:
    typedef std::iterator_traits<OutputIt> traits;

    OutputIt out_;
    std::size_t limit_;
    std::size_t count_;

  public:
    typedef std::output_iterator_tag iterator_category;
    typedef typename traits::value_type value_type;
    typedef typename traits::difference_type difference_type;
    typedef typename traits::pointer pointer;
    typedef typename traits::reference reference;

    fmt_fixed_size_buffer_iterator(OutputIt out, std::size_t limit)
      : out_(out), limit_(limit), count_(0) {}

    OutputIt base() const { return out_; }
    std::size_t count() const { return count_; }

    fmt_fixed_size_buffer_iterator& operator++() {
      if (count_++ < limit_)
        ++out_;
      return *this;
    }

    fmt_fixed_size_buffer_iterator operator++(int) {
      auto it = *this;
      ++*this;
      return it;
    }

    reference operator*() const {
      if (count_ >= limit_)
        throw std::runtime_error("end of buffer");
      return *out_;
    }
};

//fsb = fixed size buffer
typedef fmt_fixed_size_buffer_iterator<char*> fsb_iterator;

class FmtWriter {
  public:
    //Write to string - always succeeds
    template<typename T>
      static inline bool write(std::string& out_str, const T v) {
        fmt::format_to(std::back_inserter(out_str), "{}", v);
        return true;
      }

    template<typename T>
      static inline bool write(std::string& out_str, const T* v, const size_t n) {
        fmt::format_to(std::back_inserter(out_str), "{}", fmt::join(v, v+n, ""));
        return true;
      }

    template<typename T>
      static inline bool write(std::string& out_str, const T* v, const size_t n, const char sep) {
        fmt::format_to(std::back_inserter(out_str), "{}", fmt::join(v, v+n, STRING_VIEW(&sep, 1u)));
        return true;
      }

    //Write to stream - always succeeds
    template<typename T>
      static inline bool write(std::ostream& out_str, const T v) {
        fmt::print(out_str, "{}", v);
        return true;
      }

    template<typename T>
      static inline bool write(std::ostream& out_str, const T* v, const size_t n) {
        fmt::print(out_str, "{}", fmt::join(v, v+n, ""));
        return true;
      }

    template<typename T>
      static inline bool write(std::ostream& out_str, const T* v, const size_t n, const char sep) {
        fmt::print(out_str, "{}", fmt::join(v, v+n, STRING_VIEW(&sep, 1u)));
        return true;
      }

    //Write to pre-allocated buffer
    //May fail if data goes past end of buffer
    //Return true if success and update offset, else only return false (doesn't update offset)
    template<typename T>
      static inline bool write_if_space_available(char* buffer, const size_t limit, size_t& offset, const T v) {
        const auto remaining = limit-offset;
        auto result = fmt::format_to_n(fsb_iterator(buffer+offset, remaining), remaining, "{}", v);
        const auto no_overflow = result.size < remaining;
        offset += no_overflow*result.size;
        return no_overflow;
      }

    template<typename T>
      static inline bool write_if_space_available(char* buffer, const size_t limit, size_t& offset,
          const T* v, const size_t n) {
        const auto remaining = limit-offset;
        auto result = fmt::format_to_n(fsb_iterator(buffer+offset, remaining), remaining, "{}", fmt::join(v, v+n,
              ""));
        const auto no_overflow = result.size < remaining;
        offset += no_overflow*result.size;
        return no_overflow;
      }

    template<typename T>
      static inline bool write_if_space_available(char* buffer, const size_t limit, size_t& offset,
          const T* v, const size_t n, const char sep) {
        const auto remaining = limit-offset;
        auto result = fmt::format_to_n(fsb_iterator(buffer+offset, remaining), remaining, "{}", fmt::join(v, v+n,
              STRING_VIEW(&sep, 1u)));
        const auto no_overflow = result.size < remaining;
        offset += no_overflow*result.size;
        return no_overflow;
      }
};

//Template specialization for char*
template<>
inline bool FmtWriter::write<char>(std::string& out_str, const char* v, const size_t n) {
  out_str.append(v, n);
  return true;
}

//Template specialization for char*
template<>
inline bool FmtWriter::write<char>(std::ostream& out_str, const char* v, const size_t n) {
  out_str.write(v, n);
  return true;
}

//Template specialization for char*
template<>
inline bool FmtWriter::write_if_space_available<char>(char* buffer, const size_t limit, size_t& offset,
    const char* v, const size_t n) {
  const auto remaining = limit-offset;
  if(remaining >= n) {
    memcpy(buffer+offset, v, n);
    offset += n;
    return true;
  }
  else
    return false;
}

#endif
