/**
 * The MIT License (MIT)
 * Copyright (c) 2020-2021 Omics Data Automation, Inc.
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
 */

#include "hfile_genomicsdb.h"

#include "genomicsdb_logger.h"
#include "tiledb.h"
#include "tiledb_storage.h"

#include <fcntl.h>
#include <map>
#include <mutex>
#include <regex>

typedef std::map<std::string, TileDB_CTX *> tiledb_ctx_map_t;
tiledb_ctx_map_t *create_tiledb_ctx_map();
void free_tiledb_ctx_map(tiledb_ctx_map_t *tiledb_ctx_map);
thread_local std::unique_ptr<tiledb_ctx_map_t, void(*)(tiledb_ctx_map_t*)> tiledb_ctx_map(create_tiledb_ctx_map(),
                                                                                          free_tiledb_ctx_map);
tiledb_ctx_map_t *create_tiledb_ctx_map() {
  return new tiledb_ctx_map_t;
}

void free_tiledb_ctx_map(tiledb_ctx_map_t *tiledb_ctx_map) {
  logger.debug("Deleting htslib_fs_adapter tiledb_ctx_map");
  for (auto tiledb_ctx_map_entry : *tiledb_ctx_map) {
    TileDB_CTX *tiledb_ctx= tiledb_ctx_map_entry.second;
    tiledb_ctx_finalize(tiledb_ctx);
  }
  delete tiledb_ctx_map;
}

std::once_flag initialize_once;
void genomicsdb_htslib_plugin_initialize() {
  std::call_once(initialize_once, [&](){
      // Initialize htslib before registering ourselves as a plugin
      auto htslib_scheme_handler = find_scheme_handler("az://dummy");
      if (htslib_scheme_handler) {
        // Register genomicsdb supported schemes
        hfile_plugin_init(NULL);
      }
    });
}

std::string get_uri_prefix(std::string uri) {
  std::regex pattern("(.*://.*?/)(.*)");
  std::smatch match;
  if (std::regex_match(uri, match, pattern) && match.size() >= 2) {
    return match[1].str();
  } else {
    return uri;
  }
}

void *genomicsdb_filesystem_init(const char *filename, int mode) {
  if (filename && filename[0] != '\0') {
    TileDB_CTX* tiledb_ctx;
    std::string uri_prefix = get_uri_prefix(filename);
    auto found = tiledb_ctx_map->find(uri_prefix);
    if (found == tiledb_ctx_map->end()) {
      TileDB_Config tiledb_config = {};
      tiledb_config.home_ = filename;
      if (tiledb_ctx_init(&tiledb_ctx, &tiledb_config)) {
        logger.error("htslib_plugin could not open file {} {}", filename,
                     tiledb_errmsg[0]?std::string("\n")+tiledb_errmsg:"");
        errno = EIO;
        return NULL;
      }
      tiledb_ctx_map->insert({uri_prefix, tiledb_ctx});
    } else {
      tiledb_ctx = found->second;
    }
    if ((mode & O_ACCMODE) == O_RDONLY) {
      // Check file exists for O_RDONLY, otherwise return NULL
      if (is_file(tiledb_ctx, filename)) {
	return tiledb_ctx;
      }
    } else {
      return tiledb_ctx;
    }
  }
  return NULL;
}

size_t genomicsdb_filesize(void *context, const char *filename) {
  auto tiledb_ctx = reinterpret_cast<TileDB_CTX *>(context);
  return is_file(tiledb_ctx, filename)?file_size(tiledb_ctx, filename):0;
}

ssize_t genomicsdb_filesystem_read(void *context, const char *filename, off_t offset, void *buffer, size_t length) {
  if (read_file(reinterpret_cast<TileDB_CTX *>(context), filename, offset, buffer, length)) {
    logger.error("hts_plugin read {} error {}", filename, tiledb_errmsg[0]?tiledb_errmsg:"");
    return -1;
  } else {
    return length;
  }
}

ssize_t genomicsdb_filesystem_write(void *context, const char *filename, const void *buffer, size_t nbytes) {
  if (write_file(reinterpret_cast<TileDB_CTX *>(context), filename, buffer, nbytes)) {
    logger.error("hts_plugin write {} error {}", filename, tiledb_errmsg[0]?tiledb_errmsg:"");
    return -1;
  } else {
    return nbytes;
  }
}

int genomicsdb_filesystem_close(void *context, const char *filename) {
  return close_file(reinterpret_cast<TileDB_CTX *>(context), filename);
}


