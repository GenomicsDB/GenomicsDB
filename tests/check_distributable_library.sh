#!/bin/bash
#
# check_distributable_library.sh
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
# Checks that a distributable native library loads, and finds CA certificates, wherever it is meant to. On Linux it
# may depend only on glibc and zlib, plus any library given with --allow, and must not export the symbols of the
# static libraries linked into it; the newest glibc it needs is reported, and is that of the system it was built on
# at most. On macOS it may depend only on the system's own libraries and must target the minimum macOS version given
# with --macos-minimum. Either way, its OpenSSL must look for CA certificates under /etc/ssl.
#
# usage: check_distributable_library.sh <library> [--allow <library name>]... [--macos-minimum <version>]

set -euo pipefail

LIBRARY=$1
shift
ALLOWED=()
MACOS_MINIMUM=""
while [ $# -gt 0 ]; do
  case $1 in
    --allow) ALLOWED+=("$2"); shift 2 ;;
    --macos-minimum) MACOS_MINIMUM=$2; shift 2 ;;
    *) echo "Unknown argument $1" >&2; exit 2 ;;
  esac
done

failed=0
fail() {
  echo "FAILED: $*" >&2
  failed=1
}

if [ "$(uname)" = Darwin ]; then
  # otool heads its output, and each architecture's in a universal library, with a line ending in a colon
  install_name=$(otool -D "$LIBRARY" | tail -1)
  unexpected=$(otool -L "$LIBRARY" | grep -v ':$' | awk '{print $1}' | grep -vxF "$install_name" \
    | grep -vE '^(/usr/lib/|/System/Library/)' | sort -u || true)
  [ -z "$unexpected" ] || fail "$LIBRARY depends on libraries macOS does not provide:" $unexpected
  minos=$(otool -l "$LIBRARY" | awk '/LC_BUILD_VERSION/ {found = 1} found && $1 == "minos" {print $2; exit}')
  echo "$LIBRARY targets macOS $minos"
  if [ -z "$MACOS_MINIMUM" ]; then
    fail "No --macos-minimum given"
  elif [ "${minos%.0}" != "${MACOS_MINIMUM%.0}" ]; then
    fail "$LIBRARY targets macOS $minos, not $MACOS_MINIMUM"
  fi
else
  allowed_pattern='libc\.so\.6|libm\.so\.6|libdl\.so\.2|libpthread\.so\.0|librt\.so\.1'
  allowed_pattern+='|ld-linux-[-a-z0-9_]+\.so\.[0-9]+|libz\.so\.1'
  for library in "${ALLOWED[@]+"${ALLOWED[@]}"}"; do
    allowed_pattern+="|${library//./\\.}"
  done
  unexpected=$(readelf -d "$LIBRARY" | sed -n 's/.*(NEEDED).*\[\(.*\)\]/\1/p' | grep -vxE "$allowed_pattern" || true)
  [ -z "$unexpected" ] || fail "$LIBRARY depends on libraries other than glibc and zlib:" $unexpected
  # OpenSSL's, libcurl's and libuuid's
  exported=$(nm -D --defined-only "$LIBRARY" | awk '{print $NF}' \
    | grep -E '^(OPENSSL_|SSL_|SSL3_|TLS_|DTLS_|CRYPTO_|EVP_|BIO_|X509_|ERR_|curl_|Curl_|uuid_)' || true)
  [ -z "$exported" ] || fail "$LIBRARY exports $(echo "$exported" | wc -l) symbols of the static libraries" \
    "linked into it, e.g." $(echo "$exported" | head -5)
  glibc=$(objdump -T "$LIBRARY" | grep -oE 'GLIBC_[0-9]+(\.[0-9]+)+' | sort -uV | tail -1)
  echo "$LIBRARY needs ${glibc:-no versioned glibc symbols};" \
    "this system has glibc $(getconf GNU_LIBC_VERSION | awk '{print $2}')"
fi

# OpenSSL looks for CA certificates under the OPENSSLDIR it was built with
openssl_dirs=$(strings -a "$LIBRARY" | sed -n 's/^OPENSSLDIR: "\(.*\)"$/\1/p' | sort -u)
[ "$openssl_dirs" = /etc/ssl ] \
  || fail "$LIBRARY has OpenSSL look for CA certificates in" ${openssl_dirs:-no directory} "rather than /etc/ssl"

[ $failed -eq 0 ] && echo "$LIBRARY is distributable"
exit $failed
