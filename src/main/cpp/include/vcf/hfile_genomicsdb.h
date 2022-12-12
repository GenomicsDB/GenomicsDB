/**
 * The MIT License (MIT)
 * Copyright (c) 2020-2021 Omics Data Automation, Inc.
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

#ifndef __GENOMICSDB_HTSLIB_FS_PLUGIN
#define  __GENOMICSDB_HTSLIB_FS_PLUGIN

#include "hfile_internal.h"

void genomicsdb_htslib_plugin_initialize();

#ifdef __cplusplus
extern "C" {
#endif

/**@{*/
/** C Library export. */
#if (defined __GNUC__ && __GNUC__ >= 4) || defined __INTEL_COMPILER
#  define GENOMICSDB_EXPORT __attribute__((visibility("default")))
#else
#  define GENOMICSDB_EXPORT
#endif
/**@}*/

  GENOMICSDB_EXPORT void *genomicsdb_filesystem_init(const char *filename, int mode);

  GENOMICSDB_EXPORT size_t genomicsdb_filesize(void *context, const char *filename);
  
  GENOMICSDB_EXPORT ssize_t genomicsdb_filesystem_read(void *fs, const char *filename, off_t offset, void *buffer, size_t length);

  GENOMICSDB_EXPORT ssize_t genomicsdb_filesystem_write(void *fs, const char *filename, const void *buffer, size_t nbytes);
  
  GENOMICSDB_EXPORT int genomicsdb_filesystem_close(void *fs, const char *filename);

  const struct hFILE_scheme_handler *find_scheme_handler(const char *s);

  int hfile_plugin_init(struct hFILE_plugin *self);

#undef GENOMICSDB_EXPORT
#ifdef __cplusplus
}
#endif

#endif /*__GENOMICSDB_HTSLIB_FS_PLUGIN*/
