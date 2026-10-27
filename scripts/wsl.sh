#!/bin/bash
# Run a command inside WSL from the repo root, from Git Bash on Windows.
#   bash scripts/wsl.sh scripts/build_wsl.sh
#   bash scripts/wsl.sh build/cvrpbcp-Release/cvrpbcp solve instances/cvrplib/A/A-n32-k5.vrp
#   bash scripts/wsl.sh bash -c 'for s in E A P; do echo $s; done'
# wsl.exe re-parses its command line, so the command is written to a temp
# script (reachable from WSL under /mnt) and bash runs that file. Handles the
# WSL PATH (cmake/ninja in ~/.local/bin) and Git Bash path conversion.
HERE="$(cd "$(dirname "$0")/.." && pwd)"                 # e.g. /d/ARTC/AIRouting/exactsolv
DRIVE="$(echo "$HERE" | sed -E 's|^/([a-zA-Z])/.*|\1|' | tr 'A-Z' 'a-z')"
REST="$(echo "$HERE" | sed -E 's|^/[a-zA-Z]/||')"
WSLROOT="/mnt/$DRIVE/$REST"
if [ $# -eq 0 ]; then echo "usage: scripts/wsl.sh <command...>" >&2; exit 2; fi
mkdir -p "$HERE/build/wslcmd"
TMP="$HERE/build/wslcmd/cmd.$$.$RANDOM.sh"
{
  echo 'export PATH="$HOME/.local/bin:$PATH"'
  printf 'cd %q || exit 1\n' "$WSLROOT"
  printf '%q ' "$@"
  echo
} > "$TMP"
WSLTMP="$WSLROOT/build/wslcmd/$(basename "$TMP")"
export MSYS_NO_PATHCONV=1 MSYS2_ARG_CONV_EXCL='*'
wsl -- bash "$WSLTMP"
rc=$?
rm -f "$TMP"
exit $rc
