#!/usr/bin/env bash
# Build + push multi-arch (amd64/arm64) per-module images for nf-xen.
# Requires: docker buildx, logged in to Docker Hub as saditya88.
set -euo pipefail
TAG=${TAG:-1.0.0}
cd "$(dirname "$0")"
for m in qc clustering spatial celltyping_celltypist celltyping_singler celldex_download report_data; do
    echo ">>> building nf-xen-$m:$TAG"
    cp "../modules/local/$m/environment.yml" "$m/environment.yml"
    docker buildx build \
        --platform linux/amd64,linux/arm64 \
        -t "saditya88/nf-xen-$m:$TAG" \
        --push \
        "$m"
done
echo "all images pushed"
