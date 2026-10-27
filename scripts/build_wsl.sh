#!/bin/bash
# Build bplex (out-of-tree, MKL off) and cvrpbcp under WSL.
#   scripts/build_wsl.sh [Release|Debug] [extra cmake args for cvrpbcp]
set -e
export PATH="$HOME/.local/bin:$PATH"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BT=${1:-Release}; shift || true
cd "$HERE"
if [ ! -f build/bplex/libbplex.a ] || [ "$REBUILD_BPLEX" = 1 ]; then
  cmake -S ../bplex -B build/bplex -G Ninja -DBPLEX_USE_MKL=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  cmake --build build/bplex
fi
cmake -S . -B build/cvrpbcp-$BT -G Ninja -DCMAKE_BUILD_TYPE=$BT -DBPLEX_BUILD="$HERE/build/bplex" "$@"
cmake --build build/cvrpbcp-$BT
echo "built: build/cvrpbcp-$BT/cvrpbcp"
