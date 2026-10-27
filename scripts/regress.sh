#!/bin/bash
# Quick regression (inside WSL): selftest + tests/regression.txt instances.
#   scripts/regress.sh [time limit per instance, default 300]
# Exit status 1 if any instance is not proven to its expected optimum.
TL=${1:-300}
HERE="$(cd "$(dirname "$0")/.." && pwd)"
B="$HERE/build/cvrpbcp-Release/cvrpbcp"
[ -x "$B" ] || { echo "build first: scripts/build_wsl.sh" >&2; exit 1; }
cd "$HERE" || exit 1
fail=0
"$B" selftest 8 | tail -1 | grep -q "8/8 passed" || { echo "FAIL selftest"; fail=1; }
while IFS='|' read -r inst opts expect; do
  inst=$(echo "$inst" | xargs); opts=$(echo "$opts" | xargs); expect=$(echo "$expect" | xargs)
  [ -z "$inst" ] && continue
  case "$inst" in \#*) continue;; esac
  line=$(timeout $((TL + 30)) "$B" solve "$inst" --time "$TL" --log 0 $opts 2>&1 | grep "^RESULT")
  status=$(echo "$line" | sed -E 's/.*status ([a-z_]+).*/\1/')
  ub=$(echo "$line" | sed -E 's/.*\bub (-?[0-9]+).*/\1/')
  t=$(echo "$line" | sed -E 's/.*\btime ([0-9.]+)s.*/\1/')
  if [ "$status" = "optimal" ] && [ "$ub" = "$expect" ]; then
    printf "ok   %-45s %8s  %7ss\n" "$inst" "$ub" "$t"
  else
    printf "FAIL %-45s status=%s ub=%s expected=%s\n" "$inst" "$status" "$ub" "$expect"
    fail=1
  fi
done < tests/regression.txt
[ $fail -eq 0 ] && echo "regression: all passed" || echo "regression: FAILURES"
exit $fail
