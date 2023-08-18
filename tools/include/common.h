/**
 * The MIT License (MIT)
 * Copyright (c) 2021 Omics Data Automation, Inc.
 * Copyright (c) 2023 dātma, inc™
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
 **/

#ifndef GENOMICSDB_TOOLS_COMMON_H
#define GENOMICSDB_TOOLS_COMMON_H

#include "genomicsdb_logger.h"
#include "tiledb_utils.h"

#include <getopt.h>

#define OK   0
#define ERR -1

enum GenomicsDBArgsEnum {
  VERSION=1000,
};

#define LOG_TILEDB_ERRMSG if (strnlen(tiledb_errmsg, TILEDB_ERRMSG_MAX_LEN) > 0) \
                            g_logger.error("{}", tiledb_errmsg)

#endif /* GENOMICSDB_TOOLS_COMMON_H */
