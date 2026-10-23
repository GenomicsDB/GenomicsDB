#!/bin/bash

# The MIT License (MIT)
# Copyright (c) 2026 Tim Fennell
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Makes the bundle the Maven Central Portal takes for a release of org.genomicsdb:genomicsdb and
# org.genomicsdb:genomicsdb-spark: each one's jar, pom, sources and javadoc jars in Maven's repository layout, each
# signed with GPG and with MD5 and SHA-1 checksums. genomicsdb-spark gets a pom of its own (spark_pom.py) and
# genomicsdb's sources and javadoc jars, which cover both. GPG's default key signs; GPG_PASSPHRASE, if set, unlocks
# it, else gpg-agent asks for the passphrase.
#
# usage: central_bundle.sh <version> <dir with genomicsdb-<version>{.jar,.pom,-sources.jar,-javadoc.jar} and
#                          genomicsdb-spark-<version>.jar> <bundle.zip>
set -euo pipefail
VERSION=$1
IN=$(cd "$2" && pwd)
BUNDLE=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
SCRIPTS=$(cd "$(dirname "$0")" && pwd)
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

if ! grep -q "<version>$VERSION</version>" "$IN/genomicsdb-$VERSION.pom"; then
  echo "genomicsdb-$VERSION.pom does not have version $VERSION" >&2
  exit 1
fi

gpg_options=(--batch --yes --armor --detach-sign)
if [ -n "${GPG_PASSPHRASE:-}" ]; then gpg_options+=(--pinentry-mode loopback --passphrase "$GPG_PASSPHRASE"); fi

# Signs each file and writes its checksums
sign() {
  for file in "$@"; do
    gpg "${gpg_options[@]}" --output "$file.asc" "$file"
    gpg --batch --verify "$file.asc" "$file"
    for algorithm in md5 sha1; do
      openssl dgst "-$algorithm" -r "$file" | cut -d' ' -f1 > "$file.$algorithm"
    done
  done
}

CORE="$STAGE/org/genomicsdb/genomicsdb/$VERSION"
SPARK="$STAGE/org/genomicsdb/genomicsdb-spark/$VERSION"
mkdir -p "$CORE" "$SPARK"
for suffix in .jar .pom -sources.jar -javadoc.jar; do
  cp "$IN/genomicsdb-$VERSION$suffix" "$CORE/"
done
cp "$IN/genomicsdb-spark-$VERSION.jar" "$SPARK/"
for suffix in -sources.jar -javadoc.jar; do
  cp "$IN/genomicsdb-$VERSION$suffix" "$SPARK/genomicsdb-spark-$VERSION$suffix"
done
python3 -I "$SCRIPTS/spark_pom.py" "$IN/genomicsdb-$VERSION.pom" > "$SPARK/genomicsdb-spark-$VERSION.pom"
sign "$CORE"/* "$SPARK"/*

rm -f "$BUNDLE"
(cd "$STAGE" && zip -q -r "$BUNDLE" org)
unzip -l "$BUNDLE"
