/**
 * @file   logger.cc
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
 * This file implements the Logger class
 */

#include <spdlog/sinks/stdout_color_sinks.h>

#include "logger.h"

#define LOGGER_NAME "NativeGenomicsDB"
#define LOGGER_NAME_STRING "GenomicsDB.String"

/** Global thread safe logger */
Logger logger;

Logger::Logger() {
  m_logger = spdlog::get(LOGGER_NAME);
  if (m_logger == nullptr) {
    m_logger = spdlog::stderr_color_mt(LOGGER_NAME);
    // Follow default log4j pattern for now
    m_logger->set_pattern("%H:%M:%S.%e %4!l  %n - pid=%P tid=%t %v");
#ifdef NDEBUG
    m_logger->set_level(spdlog::level::info);
#else
    m_logger->set_level(spdlog::level::debug);
#endif
  }

  setup_string_logger();
}

Logger::Logger(std::shared_ptr<spdlog::logger> logger) {
  m_logger = logger;
  setup_string_logger();
}

Logger::Logger(std::shared_ptr<spdlog::logger> logger, std::shared_ptr<spdlog::logger> string_logger) {
  m_logger = logger;
  m_string_logger = string_logger;
}

Logger::~Logger() {
  spdlog::drop_all();
}

void Logger::info(const std::string& msg, bool once_only) {
  if (!once_only || not_been_logged(msg)) { 
    m_string_logger->info(msg.c_str());
  }
}

void Logger::debug(const std::string& msg, bool once_only) {
  if (!once_only || not_been_logged(msg)) { 
    m_string_logger->debug(msg.c_str());
  }
}

void Logger::debug_only(const std::string& msg, bool once_only) {
#ifdef DEBUG
  debug(msg, once_only);
#endif
}

void Logger::warn(const std::string& msg, bool once_only) {
  if (!once_only || not_been_logged(msg)) { 
    m_string_logger->warn(msg.c_str());
  }
}

void Logger::error(const std::string& msg, bool once_only) {
  if (!once_only || not_been_logged(msg)) { 
    m_string_logger->error(msg.c_str());
  }
}

void Logger::setup_string_logger() {
  m_string_logger = spdlog::get(LOGGER_NAME_STRING);
  if (m_string_logger == nullptr) {
    m_string_logger = spdlog::stdout_color_mt(LOGGER_NAME_STRING);
    // No pattern for string logger
    m_string_logger->set_pattern("%v");
#ifdef NDEBUG
    m_string_logger->set_level(spdlog::level::info);
#else
    m_string_logger->set_level(spdlog::level::debug);
#endif
  }
}

bool Logger::not_been_logged(const std::string& msg) {
  const std::lock_guard<std::mutex> lock(m_once_only_mutex);
  if (std::find(m_once_only_list.begin(), m_once_only_list.end(), msg) == m_once_only_list.end()) {
    m_once_only_list.push_back(msg);
    return true;
  } else {
    return false;
  }
}
