/**
 * @file src/test/cpp/src/test_logger.cc
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
 * Test the Logger class
 */

#include <catch2/catch.hpp>

#include "logger.h"
#include "genomicsdb_exception.h"

#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <spdlog/sinks/ostream_sink.h>

#define TEST_STR "This is a test"
#define TEST_STR_ONCE_ONLY "This is a once_only test"

#define TEST_STR_FMT "i={} str={}"

TEST_CASE("test logger", "[test_logger_basic]") {
  logger.info(TEST_STR);
  logger.warn(TEST_STR);
  logger.error(TEST_STR);
  logger.debug(TEST_STR);

  const std::string test_str(TEST_STR);
  logger.info(test_str);
  logger.warn(test_str);
  logger.error(test_str);
  logger.debug(test_str);
  logger.debug_only(test_str);

  logger.info(std::string(TEST_STR_ONCE_ONLY), true);
  logger.info(std::string(TEST_STR_ONCE_ONLY), true);

  logger.info(logger.format("Trying format {} {}", "Hello1", "Hello2"));
}

TEST_CASE("test logger format", "[test_logger_format]") {
  logger.info(TEST_STR_FMT, 1, TEST_STR);
  logger.warn(TEST_STR_FMT, 2, TEST_STR);
  logger.error(TEST_STR_FMT, 3, TEST_STR);
  logger.debug(TEST_STR_FMT, 4, TEST_STR);
  logger.debug_only(TEST_STR_FMT, 5, TEST_STR);

  logger.info_once(TEST_STR_FMT, 6, TEST_STR_ONCE_ONLY);
  logger.info_once(TEST_STR_FMT, 6, TEST_STR_ONCE_ONLY);

  logger.warn_once(TEST_STR_FMT, 7, TEST_STR_ONCE_ONLY);
  logger.warn_once(TEST_STR_FMT, 7, TEST_STR_ONCE_ONLY);

  CHECK(logger.format(TEST_STR_FMT, 0, TEST_STR).find(TEST_STR) != std::string::npos);

  GenomicsDBException exception("Test Exception");
  CHECK_THROWS_AS(logger.fatal(exception, TEST_STR_FMT, 1, TEST_STR), GenomicsDBException);
  CHECK_THROWS_AS(logger.fatal(exception), GenomicsDBException);
}

#define CHECK_LOGGER(X, Y)                      \
  do {                                          \
      logger.X(Y);                              \
      CHECK(oss.str() == check_str);            \
      oss.str("");                              \
      } while (false)

#define CHECK_STR_LOGGER(X, Y)                  \
  do {                                          \
      with_string_logger.X(Y);                  \
      CHECK(oss.str() == check_str);            \
      oss.str("");                              \
      } while (false)                             

TEST_CASE("test explicit logger", "[test_logger_explicit]") {
  std::ostringstream oss;
  auto oss_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
  
  spdlog::logger oss_logger("oss", oss_sink);
  oss_logger.set_level(spdlog::level::debug);
  oss_logger.set_pattern("%v");
  auto oss_shared_logger = std::make_shared<spdlog::logger>(oss_logger);

  const std::string check_str =  std::string(TEST_STR)+spdlog::details::os::default_eol;

  // String logger logs to console  with this constructor
  Logger logger(oss_shared_logger);

  CHECK_LOGGER(info, TEST_STR);
  CHECK_LOGGER(warn, TEST_STR);
#ifdef DEBUG
  CHECK_LOGGER(debug, TEST_STR);
  CHECK_LOGGER(debug_only, TEST_STR);
#endif
  CHECK_LOGGER(error, TEST_STR);

  logger.info(TEST_STR_FMT, 1, TEST_STR);
  CHECK(oss.str().find(TEST_STR) != std::string::npos);
  oss.str("");

  logger.info_once(TEST_STR_FMT, 2, TEST_STR);
  CHECK(oss.str().find(TEST_STR) != std::string::npos);
  oss.str("");
  logger.info_once(TEST_STR_FMT, 2, TEST_STR);
  CHECK(oss.str() == "");

  // String logger logs to console 
  const std::string test_str(TEST_STR);
  logger.info(test_str);
  CHECK(oss.str() == "");

  Logger with_string_logger(oss_shared_logger, oss_shared_logger);
  CHECK_STR_LOGGER(info, TEST_STR);
  CHECK_STR_LOGGER(info, test_str);

  with_string_logger.info(test_str, true);
  CHECK(oss.str() == check_str);
  oss.str("");
  with_string_logger.info(test_str, true);
  CHECK(oss.str() == "");
}

