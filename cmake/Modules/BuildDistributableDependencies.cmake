#
# BuildDistributableDependencies.cmake
#
# The MIT License
#
# Copyright (c) 2026 Tim Fennell
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Builds the static, position-independent libraries a distributable library links in: OpenSSL, and on
# Linux libcurl and libuuid too, so that the library then depends at run time only on glibc and zlib.
# macOS provides libcurl and libuuid itself. They are built at configure time, since GenomicsDB, TileDB
# and TileDB's cloud SDKs all look for them when they are configured, once per toolchain in
# GENOMICSDB_TOOLCHAIN_DEPS_DIR, and are found there through OPENSSL_ROOT_DIR, CURL_PREFIX_DIR and
# LIBUUID_DIR.
#
# OpenSSL looks for its configuration and CA certificates under /etc/ssl, where macOS, Amazon Linux,
# Debian, Ubuntu and RHEL 9 keep them, and libcurl uses OpenSSL's. TileDB's cloud storage clients find
# a Linux distribution's CA bundle themselves; elsewhere, e.g. htslib reading https:// URLs on RHEL 8,
# SSL_CERT_FILE must point at it, say /etc/pki/tls/cert.pem.

set(GENOMICSDB_OPENSSL_VERSION "3.5.9" CACHE STRING "Version of OpenSSL built for a distributable library")
set(GENOMICSDB_OPENSSL_URL_HASH "SHA256=603f5602e2eef00d77fbd429d34dcd5822bb301757a1bc9cdb24c670f1eb859a"
    CACHE STRING "Hash of the OpenSSL download")
set(GENOMICSDB_CURL_VERSION "8.23.0" CACHE STRING "Version of libcurl built for a distributable library")
set(GENOMICSDB_CURL_URL_HASH "SHA256=e3de77cfa6367083c776a26647fc9f633b538c06c3baa73d931a252ac8f72b49"
    CACHE STRING "Hash of the curl download")
set(GENOMICSDB_UTIL_LINUX_VERSION "2.42.4"
    CACHE STRING "Version of util-linux, for its libuuid, built for a distributable library")
set(GENOMICSDB_UTIL_LINUX_URL_HASH "SHA256=af3241e7776964dcb6bb9a811ca7b0d93000b563e2ae2a8df8f80a7cd6e04d56"
    CACHE STRING "Hash of the util-linux download")

cmake_host_system_information(RESULT GENOMICSDB_DEPS_BUILD_JOBS QUERY NUMBER_OF_LOGICAL_CORES)

# Downloads, configures, builds and installs a dependency into INSTALL_DIR unless INSTALL_DIR already has
# CHECK_FILE. CONFIGURE_COMMAND runs in the unpacked source; the build and INSTALL_TARGET use make.
function(build_distributable_dependency)
  cmake_parse_arguments(DEP "" "NAME;URL;URL_HASH;INSTALL_DIR;CHECK_FILE;INSTALL_TARGET" "CONFIGURE_COMMAND" ${ARGN})
  if(EXISTS "${DEP_INSTALL_DIR}/${DEP_CHECK_FILE}")
    return()
  endif()
  message(STATUS "Building ${DEP_NAME} for a distributable library in ${DEP_INSTALL_DIR}")
  set(work_dir "${CMAKE_BINARY_DIR}/distributable-dependencies/${DEP_NAME}")
  file(REMOVE_RECURSE "${work_dir}")
  get_filename_component(archive "${DEP_URL}" NAME)
  file(DOWNLOAD "${DEP_URL}" "${work_dir}/${archive}" EXPECTED_HASH ${DEP_URL_HASH} TLS_VERIFY ON STATUS status)
  list(GET status 0 status_code)
  if(NOT status_code EQUAL 0)
    message(FATAL_ERROR "Could not download ${DEP_URL}: ${status}")
  endif()
  file(ARCHIVE_EXTRACT INPUT "${work_dir}/${archive}" DESTINATION "${work_dir}")
  string(REGEX REPLACE "\\.tar\\.gz$" "" source_dir "${work_dir}/${archive}")
  foreach(step configure build install)
    if(step STREQUAL "configure")
      set(command ${DEP_CONFIGURE_COMMAND})
    elseif(step STREQUAL "build")
      set(command make -j ${GENOMICSDB_DEPS_BUILD_JOBS})
    else()
      set(command make ${DEP_INSTALL_TARGET})
    endif()
    execute_process(COMMAND ${command} WORKING_DIRECTORY "${source_dir}" RESULT_VARIABLE result
                    OUTPUT_FILE "${work_dir}/${step}.log" ERROR_FILE "${work_dir}/${step}.log")
    if(NOT result EQUAL 0)
      file(STRINGS "${work_dir}/${step}.log" log_lines)
      list(LENGTH log_lines num_log_lines)
      math(EXPR first_line "${num_log_lines} > 30 ? ${num_log_lines} - 30 : 0")
      list(SUBLIST log_lines ${first_line} -1 log_tail)
      list(JOIN log_tail "\n" log_tail)
      message(FATAL_ERROR "Could not ${step} ${DEP_NAME}; the end of ${work_dir}/${step}.log:\n${log_tail}")
    endif()
  endforeach()
endfunction()

set(OPENSSL_INSTALL_DIR "${GENOMICSDB_TOOLCHAIN_DEPS_DIR}/openssl-install/${GENOMICSDB_OPENSSL_VERSION}")
#OpenSSL otherwise builds for the host's architecture, which on macOS need not be the one the library is built for
if(APPLE)
  set(OPENSSL_TARGET darwin64-${GENOMICSDB_DEPS_ARCH}-cc)
endif()
build_distributable_dependency(NAME openssl
  URL "https://github.com/openssl/openssl/releases/download/openssl-${GENOMICSDB_OPENSSL_VERSION}/openssl-${GENOMICSDB_OPENSSL_VERSION}.tar.gz"
  URL_HASH ${GENOMICSDB_OPENSSL_URL_HASH}
  INSTALL_DIR "${OPENSSL_INSTALL_DIR}" CHECK_FILE lib/libcrypto.a INSTALL_TARGET install_sw
  CONFIGURE_COMMAND perl ./Configure ${OPENSSL_TARGET} --prefix=${OPENSSL_INSTALL_DIR} --libdir=lib --openssldir=/etc/ssl
                    no-shared no-module no-apps no-tests no-docs CC=${CMAKE_C_COMPILER} -fPIC)

set(OPENSSL_ROOT_DIR "${OPENSSL_INSTALL_DIR}")
set(OPENSSL_USE_STATIC_LIBS True)
set(ENV{PKG_CONFIG_PATH} "${OPENSSL_INSTALL_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
if(APPLE)
  return()
endif()

set(CURL_INSTALL_DIR "${GENOMICSDB_TOOLCHAIN_DEPS_DIR}/curl-install/${GENOMICSDB_CURL_VERSION}")
build_distributable_dependency(NAME curl
  URL "https://curl.se/download/curl-${GENOMICSDB_CURL_VERSION}.tar.gz"
  URL_HASH ${GENOMICSDB_CURL_URL_HASH}
  INSTALL_DIR "${CURL_INSTALL_DIR}" CHECK_FILE lib/libcurl.a INSTALL_TARGET install
  CONFIGURE_COMMAND ./configure --prefix=${CURL_INSTALL_DIR} --libdir=${CURL_INSTALL_DIR}/lib
                    --disable-shared --enable-static --with-pic --with-openssl=${OPENSSL_INSTALL_DIR} --with-zlib
                    --with-ca-fallback --without-ca-bundle --without-ca-path --without-libpsl --without-libidn2
                    --without-nghttp2 --without-brotli --without-zstd --without-librtmp --without-libssh2
                    --disable-ldap --disable-docs CC=${CMAKE_C_COMPILER})

set(LIBUUID_INSTALL_DIR "${GENOMICSDB_TOOLCHAIN_DEPS_DIR}/libuuid-install/${GENOMICSDB_UTIL_LINUX_VERSION}")
string(REGEX MATCH "^[0-9]+\\.[0-9]+" _util_linux_series "${GENOMICSDB_UTIL_LINUX_VERSION}")
build_distributable_dependency(NAME libuuid
  URL "https://www.kernel.org/pub/linux/utils/util-linux/v${_util_linux_series}/util-linux-${GENOMICSDB_UTIL_LINUX_VERSION}.tar.gz"
  URL_HASH ${GENOMICSDB_UTIL_LINUX_URL_HASH}
  INSTALL_DIR "${LIBUUID_INSTALL_DIR}" CHECK_FILE lib/libuuid.a INSTALL_TARGET install
  CONFIGURE_COMMAND ./configure --prefix=${LIBUUID_INSTALL_DIR} --libdir=${LIBUUID_INSTALL_DIR}/lib
                    --disable-all-programs --enable-libuuid --disable-shared --enable-static --with-pic
                    --disable-nls --without-python --without-systemd --disable-bash-completion
                    CC=${CMAKE_C_COMPILER})

set(CURL_PREFIX_DIR "${CURL_INSTALL_DIR}")
set(LIBUUID_DIR "${LIBUUID_INSTALL_DIR}")
#For dependencies that look them up with pkg-config, e.g. TileDB's Azure client for libuuid
set(ENV{PKG_CONFIG_PATH} "${CURL_INSTALL_DIR}/lib/pkgconfig:${LIBUUID_INSTALL_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
