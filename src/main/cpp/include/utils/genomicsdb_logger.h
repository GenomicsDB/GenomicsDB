/**
 * @file   genomicsdb_logger.h
 *
 * @section LICENSE
 *
 * The MIT License
 * 
 * @copyright Copyright (c) 2020 Omics Data Automation, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @section DESCRIPTION
 *
 * This file defines the Logger class
 */

#ifndef GENOMICSDB_LOGGER
#define GENOMICSDB_LOGGER

// Override spdlog level names to upper case for consistency with log4j from Java
#if !defined(SPDLOG_LEVEL_NAMES)
#define SPDLOG_LEVEL_NAMES { "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL", "OFF" }
#endif

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/fmt/fmt.h>

#include <exception>
#include <execinfo.h>
#include <list>
#include <mutex>
#include <sstream>

//TODO: Prototype from TileDB/utils.h for now.
bool is_env_set(const std::string& name);

class Logger {
 public:
  Logger();
  Logger(std::shared_ptr<spdlog::logger> logger);
  Logger(std::shared_ptr<spdlog::logger> logger, std::shared_ptr<spdlog::logger> string_logger);
  ~Logger();

  /** Direct, formatless logging of messages */
  void info(const std::string& msg, bool once_only=false);
  void debug(const std::string& msg, bool once_only=false);
  void debug_only(const std::string& msg, bool once_only=false);
  void warn(const std::string& msg, bool once_only=false);
  void error(const std::string& msg, bool once_only=false);

  /* log4j pattern logging of messages */
  template<typename... Args>
  void info(const char* fmt, const Args &... args) {
    m_logger->info(fmt, args...);
  }

  template<typename... Args>
  void debug(const char* fmt, const Args &... args) {
    m_logger->debug(fmt, args...);
  }

  template<typename... Args>
  void debug_only(const char* fmt, const Args &... args) {
#ifdef DEBUG
    m_logger->debug(fmt, args...);
#endif
  }

  template<typename... Args>
  void warn(const char* fmt, const Args &... args) {
    m_logger->warn(fmt, args...);
  }

  template<typename... Args>
  void error(const char* fmt, const Args &... args) {
    m_logger->error(fmt, args...);
  }

#define BACKTRACE_LENGTH 10
  void print_backtrace() {
    if (is_env_set("GENOMICSDB_PRINT_STACKTRACE")) {
      void *buffer[BACKTRACE_LENGTH];
      int nptrs = backtrace(buffer, BACKTRACE_LENGTH);
      char **strings = backtrace_symbols(buffer, nptrs);
      m_logger->error("Native Stack Trace:");
      for (auto i = 1; i < nptrs; i++) {
	m_string_logger->error(std::string("\t")+strings[i]);
      }
      free(strings);
    }
  }

  template<typename T, typename... Args>
  void fatal(const T& exception, const char* fmt, const Args &... args) {
    static_assert(std::is_base_of<std::exception, T>::value, "Template class to fatal() must derive from std::exception");
    m_logger->error(fmt, args...);
    print_backtrace();
    throw exception;
  }

  template<typename T>
  void fatal(const T& exception) {
    static_assert(std::is_base_of<std::exception, T>::value, "Template class to fatal() must derive from std::exception");
    m_logger->error(exception.what());
    print_backtrace();
    throw exception;
  }

  template<typename... Args>
  void info_once(const char* fmt, const Args &... args) {
    if (not_been_logged(format(fmt, args...))) {
      m_logger->info(fmt, args...);
    }
  }

  template<typename... Args>
  void warn_once(const char* fmt, const Args &... args) {
    if (not_been_logged(format(fmt, args...))) {
      m_logger->warn(fmt, args...);
    }
  }

  template<typename... Args>
  const std::string format(const char* fmt, const Args &... args) {
    return fmt::format(fmt, args...);
  }

  static std::shared_ptr<spdlog::logger> get_logger(const std::string& name);
  
 private:
  std::shared_ptr<spdlog::logger> m_logger;
  std::shared_ptr<spdlog::logger> m_string_logger;
  std::mutex m_once_only_mutex;
  std::list<std::string> m_once_only_list;

  void setup_string_logger();
  bool not_been_logged(const std::string& msg);
};

/** Global thread safe logger */
extern Logger logger;

#endif /* GENOMICSDB_LOGGER */
