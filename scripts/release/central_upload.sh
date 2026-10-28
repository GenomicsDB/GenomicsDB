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
# Uploads a bundle made by central_bundle.sh to the Maven Central Portal, to be published as soon as it validates, and
# waits until the Portal has validated it and started publishing it. Authenticates with a Portal user token from
# CENTRAL_TOKEN_USERNAME and CENTRAL_TOKEN_PASSWORD. Needs curl and jq.
#
# usage: central_upload.sh <bundle.zip> <deployment name>
set -euo pipefail
BUNDLE=$1
NAME=$2
PORTAL=https://central.sonatype.com
token=$(printf '%s:%s' "$CENTRAL_TOKEN_USERNAME" "$CENTRAL_TOKEN_PASSWORD" | base64 | tr -d '\n')

# POSTs to the Portal and prints the response, or the Portal's error on stderr. The header goes to curl on stdin so
# that the token stays out of the process list.
portal() {
  local response
  if ! response=$(printf 'Authorization: Bearer %s\n' "$token" | curl --fail-with-body -sS -X POST -H @- "$@"); then
    echo "$response" >&2
    return 1
  fi
  echo "$response"
}

id=$(portal -F "bundle=@$BUNDLE" "$PORTAL/api/v1/publisher/upload?name=$NAME&publishingType=AUTOMATIC")
echo "Uploaded $BUNDLE as deployment $id"
# Validation takes minutes. Sonatype asks for gentle polling, and has cut off clients that poll hard.
for _ in $(seq 90); do
  status=$(portal "$PORTAL/api/v1/publisher/status?id=$id")
  state=$(jq -r .deploymentState <<< "$status")
  echo "$(date -u +%T) $state"
  case "$state" in
    PUBLISHING|PUBLISHED)
      echo "Validated; the Portal is publishing it to Maven Central, which can take a while longer"
      exit 0
      ;;
    FAILED)
      jq . <<< "$status" >&2
      exit 1
      ;;
  esac
  sleep 30
done
echo "Still $state after 45 minutes; see $PORTAL/publishing/deployments" >&2
exit 1
