#!/bin/bash

# This script is meant to be run from develop or master to add a version+commit hash tag.
# The tag will trigger the docker release workflow on Github Actions

REPO_ROOT=$(git rev-parse --show-toplevel)

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [[ "$CURRENT_BRANCH" != "develop" && "$CURRENT_BRANCH" != "master" ]]; then
    echo "Release should be run from develop or master branch!"
    exit 1
fi

VERSION=$(grep "GenomicsDB release version" $REPO_ROOT/CMakeLists.txt | cut -d '"' -f 2)
GIT_COMMIT_HASH=$(git log -1 --format=%h)

TAG=v${VERSION}-${GIT_COMMIT_HASH}

echo "Adding tag: $TAG"

git tag -am "Tagging for docker release $TAG" $TAG
git push origin $TAG
