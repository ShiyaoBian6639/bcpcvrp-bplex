#!/bin/bash
# Benchmark runner (inside WSL):
#   scripts/bench.sh <instance dir> <max n> <time limit s> <output prefix> [extra solve options]
# Runs cvrpbcp solve on every .vrp in the directory with n <= max n, appends
# one JSON line per run to <prefix>.jsonl and the full log to <prefix>.log.
# The binary is copied to <prefix>.bin first so rebuilding during a long
# benchmark does not change the program under test (BIN=... overrides the
# source binary). DEBUGSOL=1 loads each instance's .sol as the known
# solution: the solver then prints "DEBUG ..." lines whenever a node loses it
# (grep DEBUG <prefix>.log afterwards; none expected).
set -u
DIR=${1:?instance dir}; MAXN=${2:?max n}; TL=${3:?time limit}; OUT=${4:?output prefix}; shift 4
HERE="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${BIN:-$HERE/build/cvrpbcp-Release/cvrpbcp}"
[ -x "$SRC" ] || { echo "build first: scripts/build_wsl.sh" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
cp "$SRC" "$OUT.bin"
B="$OUT.bin"
: > "$OUT.log"
for f in "$DIR"/*.vrp; do
  dim=$(grep -m1 -i "DIMENSION" "$f" | tr -dc '0-9')
  n=$((dim - 1))
  [ "$n" -le "$MAXN" ] || continue
  echo "== $f (n=$n)" >> "$OUT.log"
  sol="${f%.vrp}.sol"
  if [ "${DEBUGSOL:-0}" = 1 ] && [ -f "$sol" ]; then
    CVRPBCP_DEBUGSOL="$sol" timeout $((TL + 60)) "$B" solve "$f" --time "$TL" --log 1 --jsonl "$OUT.jsonl" "$@" >> "$OUT.log" 2>&1
  else
    timeout $((TL + 60)) "$B" solve "$f" --time "$TL" --log 1 --jsonl "$OUT.jsonl" "$@" >> "$OUT.log" 2>&1
  fi
  grep "^RESULT" "$OUT.log" | tail -1
done
echo "DEBUG lines: $(grep -c DEBUG "$OUT.log")"
