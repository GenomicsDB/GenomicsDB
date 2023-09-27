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

#include "hfile_internal.h"
#include "hfile_genomicsdb.h"

typedef struct {
  hFILE base;
  void *context;
  const char *filename;
  const char *mode;
  off_t offset;
  size_t length;
} hFILE_genomicsdb;

static ssize_t genomicsdb_read(hFILE* fpv, void *buffer, size_t nbytes) {
  hFILE_genomicsdb *fp = (hFILE_genomicsdb *)fpv;
  size_t avail = fp->length - fp->offset;
  if (nbytes > avail) nbytes = avail;
  if (nbytes) {
    ssize_t length = genomicsdb_filesystem_read(fp->context, fp->filename, fp->offset, buffer, nbytes);
    fp->offset += length;
    return length;
  }
  return nbytes;
}

static ssize_t genomicsdb_write(hFILE* fpv, const void *buffer, size_t nbytes) {
  hFILE_genomicsdb *fp = (hFILE_genomicsdb *)fpv;
  return genomicsdb_filesystem_write(fp->context, fp->filename, buffer, nbytes);
}

static off_t genomicsdb_seek(hFILE *fpv, off_t offset, int whence) {
  hFILE_genomicsdb *fp = (hFILE_genomicsdb *)fpv;
           
  size_t absoffset = (offset >= 0)? offset : -offset;
  size_t origin;

  switch (whence) {
    case SEEK_SET: origin = 0; break;
    case SEEK_CUR: origin = fp->offset; break;
    case SEEK_END: origin = fp->length; break;
    default: errno = EINVAL; return -1;
  }

  if ((offset  < 0 && absoffset > origin) ||
      (offset >= 0 && absoffset > fp->length - origin)) {
    errno = EINVAL;
    return -1;
  }

  fp->offset = origin + offset;
  return fp->offset;
}

static int genomicsdb_close(hFILE *fpv) {
  hFILE_genomicsdb *fp = (hFILE_genomicsdb *)fpv;
  return genomicsdb_filesystem_close(fp->context, fp->filename);
}

static const struct hFILE_backend genomicsdb_backend = {
    genomicsdb_read,
    genomicsdb_write,
    genomicsdb_seek,
    NULL,
    genomicsdb_close};

static hFILE *genomicsdb_open(const char *uri, const char *mode) {
  hFILE_genomicsdb *fp = (hFILE_genomicsdb *)hfile_init(sizeof(hFILE_genomicsdb), mode, 0);
  if (fp == NULL) {
    // TODO: error handling
    return NULL;
  }
  fp->context = genomicsdb_filesystem_init(uri, hfile_oflags(mode));
  if (fp->context == NULL) {
    free(fp);
    return NULL;
  }
  fp->filename = uri;
  fp->mode = mode;
  fp->offset = 0;
  fp->length = genomicsdb_filesize(fp->context, uri);
  fp->base.backend = &genomicsdb_backend;
  return &fp->base;
}

int hfile_plugin_init(struct hFILE_plugin *self) {
   static const struct hFILE_scheme_handler genomicsdb_handler =
       { genomicsdb_open, hfile_always_remote, "GenomicsDB Storage", 99 };
   
  if (self) {
    self->name = "GenomicsDB Storage";
  }
  hfile_add_scheme_handler("az", &genomicsdb_handler);
  hfile_add_scheme_handler("azb", &genomicsdb_handler);
  hfile_add_scheme_handler("hdfs", &genomicsdb_handler);
  hfile_add_scheme_handler("s3", &genomicsdb_handler);
  hfile_add_scheme_handler("gs", &genomicsdb_handler);
  return 0;
}
