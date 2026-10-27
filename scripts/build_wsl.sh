#!/bin/bash
# Build bplex (out-of-tree, MKL off), install it to build/bplex-install and
# build cvrpbcp against that package under WSL.
#   scripts/build_wsl.sh [Release|Debug] [extra cmake args for cvrpbcp]
# REBUILD_BPLEX=1 rebuilds and reinstalls bplex (after a bplex update).
set -e
export PATH="$HOME/.local/bin:$PATH"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BT=${1:-Release}; shift || true
cd "$HERE"
PREFIX="$HERE/build/bplex-install"
if [ ! -f "$PREFIX/lib/cmake/bplex/bplexConfig.cmake" ] || [ "$REBUILD_BPLEX" = 1 ]; then
  cmake -S ../bplex -B build/bplex -G Ninja -DBPLEX_USE_MKL=OFF -DBPLEX_BUILD_TESTS=OFF -DBPLEX_INSTALL=ON \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  cmake --build build/bplex
  cmake --install build/bplex --prefix "$PREFIX" > /dev/null
fi
cmake -S . -B build/cvrpbcp-$BT -G Ninja -DCMAKE_BUILD_TYPE=$BT -DCMAKE_PREFIX_PATH="$PREFIX" "$@"
cmake --build build/cvrpbcp-$BT
echo "built: build/cvrpbcp-$BT/cvrpbcp"
