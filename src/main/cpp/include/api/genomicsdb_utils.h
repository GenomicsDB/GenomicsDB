/**
 * @file genomicsdb_utils.h
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2022-2023 Omics Data Automation, Inc.
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
 *
 * @section DESCRIPTION
 *
 * Header file to genomicsdb utilities
 *
 **/

#pragma once

// Override global project visibility which is set to hidden
#define GENOMICSDB_EXPORT __attribute__((visibility("default")))

#include <string>
#include <vector>

#include "genomicsdb_status.h"

namespace genomicsdb {

/**
 * Get version of the native GenomicsDB Library
 * @return version string of GenomicsDB Library 
 */
GENOMICSDB_EXPORT std::string version();

/**
 * Check if file referenced by filename exists as a file.
 * @param filename Filename could be a file in the local filesystem or could be a cloud URI(s3, gcs or azure)
 * @return true if file exists, false otherwise
 */
GENOMICSDB_EXPORT bool is_file(const std::string& filename);

/**
 * Get size of file referenced by filename
 * @param filename Filename could be a file in the local filesystem or could be a cloud URI(s3, gcs or azure)
 * @return size of file if it exists, -1 otherwise
 */
GENOMICSDB_EXPORT ssize_t file_size(const std::string& filename);

/**
 * Read contents of file referenced by filename
 * @param filename Filename could be a file in the local filesystem or could be a cloud URI(s3, gcs or azure)
 * @param buffer Reference to a buffer pointer. buffer is malloc'ed and has to be freed by calling function
 * @param length Reference to length of buffer.
 * @return GENOMICSDB_OK for success, GENOMICSDB_ERR otherwise
 */
GENOMICSDB_EXPORT int read_entire_file(const std::string& filename, void **buffer, size_t *length);

/** Check whether workspace exists
 * @param workspace storage URL
 * @return true if given workspace exists, false otherwise
 */ 
GENOMICSDB_EXPORT bool workspace_exists(const std::string& workspace);

/** Check whether an array in a given workspace exists
 * @param workspace storage URL
 * @param array_name
 * @return true if array exists in the workspace, false otherwise
 */
GENOMICSDB_EXPORT bool array_exists(const std::string& workspace, const std::string& array_name);

/**
 * Get all array names under a workspace
 * @param workspace storage URL
 * @return vector of array names in the workspace
 */
GENOMICSDB_EXPORT std::vector<std::string> get_array_names(const std::string& workspace);

/**
 * Cache fragment metadata associated with array in a workspace. Used by the underlying storage system only
 * when the workspace is a cloud URL for now.
 * @param workspace storage URL
 * @param array URL under workspace
 * @return GENOMICSDB_OK for success, GENOMICSDB_ERR otherwise
 */
GENOMICSDB_EXPORT int cache_fragment_metadata(const std::string& workspace, const std::string& array_name);

}


