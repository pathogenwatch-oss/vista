#!/usr/bin/env bash
set -euo pipefail

push=false
if [[ "${1:-}" == "--push" ]]; then
    push=true
elif [[ $# -gt 0 ]]; then
    echo "Usage: $0 [--push]" >&2
    exit 2
fi

image="902121496535.dkr.ecr.eu-west-2.amazonaws.com/pathogenwatch-source/vista"
version=$(sed -n 's/^version = "\([^"]*\)"$/\1/p' pyproject.toml)

if [[ -z "$version" ]]; then
    echo "Could not read project version from pyproject.toml" >&2
    exit 1
fi

tagged_image="$image:$version"
docker build --tag "$tagged_image" .

if [[ "$push" == true ]]; then
    docker push "$tagged_image"
else
    echo "Image has not been pushed. To push it, run:"
    echo "docker push $tagged_image"
fi
