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

"""Writes THIRD-PARTY-NOTICES.txt, which the jar carries in META-INF: the software compiled into the released native
libraries, each with its version, licence and notices. The notices come from the submodules in the source tree and,
for software the build downloads, from the upstream project at the version built, so the source tree's submodules
must be checked out and the network reachable. When a dependency of the released library is added, removed or
changes version, update COMPONENTS and rerun this.

usage: third_party_notices.py <GenomicsDB source dir> <output file>
"""
import os
import re
import sys
import textwrap
import urllib.request

GITHUB = "https://raw.githubusercontent.com"


def tree(path):
    """A file in the source tree."""
    return ("tree", path)


def url(location):
    """A file fetched from the given URL."""
    return ("url", location)


def first_comment(source):
    """The first /* ... */ comment of a source file given by tree() or url(), without the comment markers."""
    return ("first_comment", source)


# (name, version, project, SPDX licence, where it is compiled in or None for everywhere, [notices])
COMPONENTS = [
    ("TileDB", "0.6.0 (GenomicsDB/TileDB)", "https://github.com/GenomicsDB/TileDB", "MIT", None,
     [tree("dependencies/TileDB/LICENSE")]),
    ("LZ4", "1.9.2 (bundled with TileDB)", "https://github.com/lz4/lz4", "BSD-2-Clause", None,
     [tree("dependencies/TileDB/core/src/codec/external/lz4/LICENSE")]),
    ("Zstandard", "1.5.7", "https://github.com/facebook/zstd", "BSD-3-Clause", None,
     [url(f"{GITHUB}/facebook/zstd/v1.5.7/LICENSE")]),
    ("muparserx", "4.0.12 (bundled with TileDB)", "https://github.com/beltoforion/muparserx", "BSD-2-Clause", None,
     [tree("dependencies/TileDB/deps/muparserx/LICENSE")]),
    ("azure-storage-cpplite", "commit 0611aed (bundled with TileDB)", "https://github.com/Azure/azure-storage-cpplite",
     "MIT", None, [tree("dependencies/TileDB/deps/azure-storage-cpplite/LICENSE")]),
    ("Apache Hadoop libhdfs", "2.7.5, in TileDB's HDFSWrapper", "https://github.com/apache/hadoop",
     "Apache-2.0 AND MIT", "Linux only, where the library is built with HDFS support",
     [url(f"{GITHUB}/apache/hadoop/rel/release-2.7.5/NOTICE.txt"),
      url("https://www.apache.org/licenses/LICENSE-2.0.txt"),
      tree("dependencies/TileDB/deps/HDFSWrapper/LICENSE")]),
    ("HTSlib", "1.10.2 (GenomicsDB/htslib, commit 2cfd94d9)", "https://github.com/GenomicsDB/htslib",
     "MIT AND BSD-3-Clause", None, [tree("dependencies/htslib/LICENSE")]),
    ("libdeflate", "1.26", "https://github.com/ebiggers/libdeflate", "MIT", None,
     [url(f"{GITHUB}/ebiggers/libdeflate/v1.26/COPYING")]),
    ("Protocol Buffers", "3.21.7", "https://github.com/protocolbuffers/protobuf", "BSD-3-Clause", None,
     [url(f"{GITHUB}/protocolbuffers/protobuf/v21.7/LICENSE")]),
    ("spdlog", "1.8.2", "https://github.com/gabime/spdlog", "MIT", None, [tree("dependencies/spdlog/LICENSE")]),
    ("{fmt}", "7.1.3 (bundled with spdlog)", "https://github.com/fmtlib/fmt", "MIT", None,
     [tree("dependencies/spdlog/include/spdlog/fmt/bundled/LICENSE.rst")]),
    ("RapidJSON", "1.1.0 (commit 6534506e)", "https://github.com/Tencent/rapidjson", "MIT", None,
     [tree("dependencies/RapidJSON/license.txt")]),
    ("AWS SDK for C++", "1.8.187", "https://github.com/aws/aws-sdk-cpp", "Apache-2.0", None,
     [url(f"{GITHUB}/aws/aws-sdk-cpp/1.8.187/NOTICE.txt"), url(f"{GITHUB}/aws/aws-sdk-cpp/1.8.187/LICENSE.txt")]),
    ("cJSON", "bundled with the AWS SDK for C++", "https://github.com/DaveGamble/cJSON", "MIT", None,
     [first_comment(url(f"{GITHUB}/aws/aws-sdk-cpp/1.8.187/aws-cpp-sdk-core/source/external/cjson/cJSON.cpp"))]),
    ("TinyXML-2", "bundled with the AWS SDK for C++", "https://github.com/leethomason/tinyxml2", "Zlib", None,
     [first_comment(url(f"{GITHUB}/aws/aws-sdk-cpp/1.8.187/aws-cpp-sdk-core/source/external/tinyxml2/tinyxml2.cpp"))]),
    ("aws-c-common", "0.6.9", "https://github.com/awslabs/aws-c-common", "Apache-2.0", None,
     [url(f"{GITHUB}/awslabs/aws-c-common/v0.6.9/NOTICE"), url(f"{GITHUB}/awslabs/aws-c-common/v0.6.9/LICENSE")]),
    ("aws-checksums", "0.1.12", "https://github.com/awslabs/aws-checksums", "Apache-2.0", None,
     [url(f"{GITHUB}/awslabs/aws-checksums/v0.1.12/LICENSE")]),
    ("aws-c-event-stream", "0.1.5", "https://github.com/awslabs/aws-c-event-stream", "Apache-2.0", None,
     [url(f"{GITHUB}/awslabs/aws-c-event-stream/v0.1.5/NOTICE"),
      url(f"{GITHUB}/awslabs/aws-c-event-stream/v0.1.5/LICENSE")]),
    ("Google Cloud C++ Client Libraries", "1.24.0", "https://github.com/googleapis/google-cloud-cpp", "Apache-2.0",
     None, [url(f"{GITHUB}/googleapis/google-cloud-cpp/v1.24.0/LICENSE")]),
    ("Abseil", "20230802.1", "https://github.com/abseil/abseil-cpp", "Apache-2.0", None,
     [url(f"{GITHUB}/abseil/abseil-cpp/20230802.1/LICENSE")]),
    ("CRC32C", "1.1.2", "https://github.com/google/crc32c", "BSD-3-Clause", None,
     [url(f"{GITHUB}/google/crc32c/1.1.2/LICENSE")]),
    ("JSON for Modern C++", "3.9.1", "https://github.com/nlohmann/json", "MIT", None,
     [url(f"{GITHUB}/nlohmann/json/v3.9.1/LICENSE.MIT")]),
    ("OpenSSL", "3.5.9", "https://github.com/openssl/openssl", "Apache-2.0", None,
     [url(f"{GITHUB}/openssl/openssl/openssl-3.5.9/LICENSE.txt")]),
    ("curl", "8.23.0", "https://github.com/curl/curl", "curl", "Linux only; on macOS the system libcurl is used",
     [url(f"{GITHUB}/curl/curl/curl-8_23_0/COPYING")]),
    ("libuuid (util-linux)", "2.42.4", "https://github.com/util-linux/util-linux", "BSD-3-Clause",
     "Linux only; on macOS the system libuuid is used",
     [first_comment(url(f"{GITHUB}/util-linux/util-linux/v2.42.4/libuuid/src/uuid.h"))]),
    ("GCC runtime libraries (libstdc++, libgcc)", "15", "https://gcc.gnu.org",
     "GPL-3.0-or-later WITH GCC-exception-3.1",
     "Linux only; on macOS the system libc++ is used. The GCC Runtime Library Exception allows them to be linked "
     "into the library and distributed under its licence",
     [url(f"{GITHUB}/gcc-mirror/gcc/releases/gcc-15.2.0/COPYING.RUNTIME")]),
]

HEADER = """\
Third-party software in GenomicsDB's native libraries

The native libraries in this jar (libtiledbgenomicsdb.so for Linux and libtiledbgenomicsdb.dylib for macOS) are
GenomicsDB, licensed under the MIT License in META-INF/LICENSE, with the software below compiled in. Each entry gives
the version built, where it comes from and its licence, followed by its copyright and licence notices as distributed
with it.
"""

RULE = "=" * 100


def text(source_dir, notice):
    """The text of a notice: a whole file, or the first comment of one."""
    kind, location = notice
    if kind == "first_comment":
        block = re.search(r"/\*(.*?)\*/", text(source_dir, location), re.S).group(1)
        lines = [re.sub(r"^ ?\* ?", "", line) for line in block.strip("\n").split("\n")]
        return textwrap.dedent("\n".join(lines)).strip("\n")
    if kind == "tree":
        with open(os.path.join(source_dir, location), encoding="utf-8") as f:
            return f.read().strip("\n")
    with urllib.request.urlopen(location, timeout=60) as response:
        return response.read().decode("utf-8").strip("\n")


def main(source_dir, output):
    out = [HEADER]
    width = max(len(name) for name, *_ in COMPONENTS)
    for name, version, _, licence, _, _ in COMPONENTS:
        out.append(f"  {name.ljust(width)}  {version}; {licence}")
    for name, version, project, licence, platforms, notices in COMPONENTS:
        out.append(f"\n{RULE}\n{name} {version}\n{project}\nLicence: {licence}")
        if platforms:
            out.append(textwrap.fill(platforms, 100))
        for notice in notices:
            out.append(f"\n{text(source_dir, notice)}")
    with open(output, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
