#!/usr/bin/env bash
# Dependencies: CMake, C++17 compiler, curl, zlib, zstd, BLAS, LAPACK/LAPACKE.
# Set STRAW_DIR to the directory containing straw.h/straw.cpp if needed.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
args=()
if [[ -n "${STRAW_DIR:-}" ]]; then args+=("-DSTRAW_DIR=$STRAW_DIR"); fi
cmake -S "$script_dir/.." -B "$script_dir/../build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="$script_dir" "${args[@]}" "$@"
cmake --build "$script_dir/../build" --parallel
