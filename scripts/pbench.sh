#!/bin/bash
# Parallel benchmark runner (inside WSL):
#   scripts/pbench.sh <instance dir|list file> <time limit s> <output prefix> [extra solve options]
# Runs cvrpbcp solve on every instance (all .vrp of a directory, or the paths
# listed one per line in a file), JOBS at a time (default 4), with the
# instance's .sol loaded as the known solution when present (DEBUGSOL=0
# disables that). One JSON line per run goes to <prefix>.jsonl, per-instance
# logs to <prefix>.d/<name>.log, DEBUG lines are counted at the end.
# The binary is copied to <prefix>.bin first (BIN=... overrides the source).
set -u
SET=${1:?instance dir or list}; TL=${2:?time limit}; OUT=${3:?output prefix}; shift 3
HERE="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${BIN:-$HERE/build/cvrpbcp-Release/cvrpbcp}"
JOBS=${JOBS:-4}
[ -x "$SRC" ] || { echo "build first: scripts/build_wsl.sh" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")" "$OUT.d"
cp "$SRC" "$OUT.bin"
: > "$OUT.jsonl"
if [ -d "$SET" ]; then ls "$SET"/*.vrp; else cat "$SET"; fi > "$OUT.list"
export B="$OUT.bin" TL OUT DEBUGSOL="${DEBUGSOL:-1}"
export OPTS="$*"
run_one() {
  f=$1
  name=$(basename "$f" .vrp)
  sol="${f%.vrp}.sol"
  if [ "$DEBUGSOL" = 1 ] && [ -f "$sol" ]; then export CVRPBCP_DEBUGSOL="$sol"; else unset CVRPBCP_DEBUGSOL; fi
  timeout $((TL + 60)) "$B" solve "$f" --time "$TL" --log 1 --jsonl "$OUT.jsonl" $OPTS > "$OUT.d/$name.log" 2>&1
  grep "^RESULT" "$OUT.d/$name.log" | tail -1 | cut -c1-120
}
export -f run_one
xargs -a "$OUT.list" -P "$JOBS" -I{} bash -c 'run_one {}'
echo "DEBUG lines: $(cat "$OUT.d"/*.log | grep -c DEBUG)"
