/**
 * @file   logger.h
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

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>

#include <exception>
#include <list>
#include <mutex>
#include <sstream>

class Logger {
 public:
  Logger();
  Logger(std::shared_ptr<spdlog::logger> logger);
  Logger(std::shared_ptr<spdlog::logger> logger, std::shared_ptr<spdlog::logger> string_logger);
  ~Logger();

  /** Direct, formatless logging of messages */
  void info(const std::string& msg, bool once_only=false);
  void debug(const std::string& msg, bool once_only=false);
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
  void warn(const char* fmt, const Args &... args) {
    m_logger->warn(fmt, args...);
  }
  template<typename... Args>
  void error(const char* fmt, const Args &... args) {
    m_logger->error(fmt, args...);
  }
  template<typename... Args>
  void fatal(const std::exception& exception, const char* fmt, const Args &... args) {
    m_logger->error(fmt, args...);
    throw exception;
  }

  template<typename... Args>
  void info_once(const char* fmt, const Args &... args) {
    if (not_been_logged(get_message(fmt, args...))) {
      m_logger->info(fmt, args...);
    }
  }
  template<typename... Args>
  void warn_once(const char* fmt, const Args &... args) {
    if (not_been_logged(get_message(fmt, args...))) {
      m_logger->warn(fmt, args...);
    }
  }
  
 private:
  std::shared_ptr<spdlog::logger> m_logger;
  std::shared_ptr<spdlog::logger> m_string_logger;
  std::mutex m_once_only_mutex;
  std::list<std::string> m_once_only_list;

  template<typename... Args>
  const std::string get_message(const char* fmt, const Args &... args) {
    std::ostringstream oss;
    auto oss_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
    
    spdlog::logger oss_logger("oss", oss_sink);
    oss_logger.set_level(spdlog::level::info);
    oss_logger.set_pattern("%v");
    oss_logger.info(fmt, args...);
  
    return oss.str().substr(0, oss.str().length() - strlen(spdlog::details::os::default_eol));
  }
  void setup_string_logger();
  bool not_been_logged(const std::string& msg);
};

/** Global thread safe logger */
extern Logger logger;

#endif /* GENOMICSDB_LOGGER */
