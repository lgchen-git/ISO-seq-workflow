#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)
cd "$root"
image=iso-seq-bulk:1.1.0
export TMPDIR="$root/tmp/build" DOCKER_CONFIG="$root/tmp/docker-client"
mkdir -p "$TMPDIR" "$DOCKER_CONFIG" "$root/log"

case "${1:-}" in
  "") export_image=false ;;
  --export) export_image=true ;;
  --help|-h)
    echo "Usage: bash pipeline/docker/build_export.sh [--export]"
    exit 0
    ;;
  *) echo "Unknown option: $1" >&2; exit 2 ;;
esac

docker build --network=host -f pipeline/docker/Dockerfile -t "$image" . 2>&1 \
  | tee log/docker-build.log

if [[ "$export_image" == true ]]; then
  mkdir -p outputs/_docker
  archive="outputs/_docker/iso-seq-bulk_1.1.0.tar.gz"
  partial="$archive.partial.$$"
  trap 'rm -f "$partial"' EXIT
  docker save "$image" | gzip -1 > "$partial"
  mv -f "$partial" "$archive"
  trap - EXIT
fi
