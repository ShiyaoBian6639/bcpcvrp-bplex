---
name: cvrpbcp-dev
description: Build, run, benchmark and debug the cvrpbcp branch-cut-and-price solver (exactsolv) under WSL from the Windows Claude Code session. Use for any build/run/test/benchmark step in exactsolv, and before touching bplex-related code.
---

# cvrpbcp development recipes

All compilation and runs happen in WSL Ubuntu (gcc 11, cmake/ninja in
`~/.local/bin`). The Windows Bash tool cannot run WSL commands reliably inline
(Git Bash rewrites `/mnt/...` paths and heredocs containing apostrophes fail),
so always go through `scripts/wsl.sh`, which handles both.

## Rules
- `../bplex/` is read-only. Build it out-of-tree only: `scripts/build_wsl.sh`
  does `cmake -S ../bplex -B build/bplex -DBPLEX_USE_MKL=OFF`, then
  `cmake --install build/bplex --prefix build/bplex-install`; cvrpbcp finds it
  with `find_package(bplex 0.2 CONFIG)` (CMAKE_PREFIX_PATH). Set
  `REBUILD_BPLEX=1` when bplex's sources changed (new API).
- Token discipline (from `../bplex_wt/AGENT_BRIEF.md`): never Read or cat a
  whole source file (grep -n, then a range); grep logs and JSONL, never cat
  them; cap open-ended output with `| tail -n 20`; `git diff --stat` first.
- Write C++ sources with the Write/Edit tools, not shell heredocs.
- Never report an optimum that `cvrpbcp check` or the solver's own checker
  did not verify.

## Build
    bash scripts/wsl.sh scripts/build_wsl.sh               # Release -> build/cvrpbcp-Release/cvrpbcp
    bash scripts/wsl.sh scripts/build_wsl.sh Debug -DCVRPBCP_SANITIZE=ON   # ASan/UBSan -> build/cvrpbcp-Debug
    bash scripts/wsl.sh env REBUILD_BPLEX=1 scripts/build_wsl.sh            # after a bplex update

## Tuning set (medium, parallel)
    bash scripts/wsl.sh env JOBS=4 scripts/pbench.sh instances/xml100/medium 300 results/<tag> [solve options]
12 XML100 instances (two per route-size class, known optima, .sol loaded as
the known solution) run four at a time; one JSON line per run in
`results/<tag>.jsonl`, logs in `results/<tag>.d/`. Compare variants by
solved count and total time; keep the 900 s X-n153-k22 / XML200_2251_68
runs for the final check only (their trees are too noisy to rank variants).

## Run
    B=build/cvrpbcp-Release/cvrpbcp
    bash scripts/wsl.sh $B selftest 12                     # brute force vs bplex MIP vs BCP, tiny instances
    bash scripts/wsl.sh $B info  instances/cvrplib/P/P-n55-k7.vrp
    bash scripts/wsl.sh $B check instances/cvrplib/P/P-n55-k7.vrp instances/cvrplib/P/P-n55-k7.sol
    bash scripts/wsl.sh $B root  instances/cvrplib/P/P-n55-k7.vrp --log 2
    bash scripts/wsl.sh $B solve instances/cvrplib/P/P-n55-k7.vrp --time 600 --log 1 --jsonl results/dev.jsonl
Options: `--K k --ng s --knn k --cols c --alpha a --no-rcc --no-src --no-early
--mono --time T --log L --ub u --no-rmip --rmip-time s --rmip-every n --nodes n
--sb k --purge-age a --no-enum --enum-mip-max r --enum-dump prefix --no-dive
--dive-time s --no-savings --cut-age n --cleanup-every n --plunge d --qbound s
--purge-every n --stage-skip n --out f.sol`.
Defaults tuned 2026-10-27 at n=150: knn 5, 75 cols/iteration, purge every 10 CG
iterations, q-route bound in all stages (`--qbound 3` turns it off).
Fleet: CVRPLIB files state "No of trucks"; by default exactly that many are
used (matches the published optima). `--fleet max` allows fewer,
`--free-fleet` ignores it (pipeline instances: X, XML; use `--free-fleet`).
`solve` prints one `RESULT ...` line (status, lb, ub, published, nodes, time,
verdict), a `times:` breakdown, and exits 1 when a proven value differs from
the published one. `--log 2` adds per-node, cut-round and enumeration lines.

## Regression before and after any solver change
    bash scripts/wsl.sh scripts/regress.sh 300
Runs `selftest` and `tests/regression.txt` (ten instances, ~3 min) and
fails if any optimum differs. Add an instance there when a bug is fixed on
it. A wrong "optimal" is the worst outcome of this project: when one shows
up, run the instance with the known solution loaded,
    bash scripts/wsl.sh env CVRPBCP_DEBUGSOL=path/to/inst.sol build/cvrpbcp-Release/cvrpbcp solve path/to/inst.vrp --log 2
and read the `DEBUG ...` lines: they say at which node the known solution
was feasible but got pruned, declared infeasible, or lost a route in the
enumeration (with a per-prefix trace of the completion bound). Bisect with
flags (`--no-enum --no-src --mono --no-savings --no-dive --cut-age 1000000`)
to find the component. Benchmarks run with `DEBUGSOL=1` so that a lost
known solution is reported even when the final answer is right.

## Benchmark a set
    bash scripts/wsl.sh env DEBUGSOL=1 scripts/bench.sh instances/cvrplib/A 60 600 results/A
Arguments: directory, max n, time limit per instance, output prefix. Writes
`<prefix>.jsonl` (every run) and `<prefix>.log` (the binary is copied to
`<prefix>.bin` first, so rebuilding meanwhile is safe); then
    python -I scripts/summarize.py results/A.jsonl
prints the table and writes `<prefix>.usable.csv` (proven only, cap not
binding, matches published when known).

## Debug
- Crash or wrong answer: rebuild Debug with sanitizers and rerun the same
  command; the trace names the line.
- Master LP failure: `solve_node` writes `cvrpbcp_master_fail.mps`; report it
  to the bplex side with the expected status, do not patch bplex.
- Bound above a published optimum: the CLI prints `BUG: lower bound
  exceeds...`; suspect pricing (missing routes), the Lagrangian bound
  formula, or a wrong cut; run with `--no-rcc` and `--no-early` to isolate.
- Long output: the Bash tool shows only the tail; pipe through
  `grep -v '^Route' | tail -n 40` or redirect to `results/*.log`.

## Where instances are
`instances/cvrplib/{A,B,E,F,P}` (with `.sol`); X: `../xltw/benchmarks/CVRP`;
XML100 + optimal solutions: `../data/xml100/{instances,solutions}`; XML150/200:
`../exactlabel/inst`; Solomon/GH CVRPTW with BKS: `../xltw/benchmarks/VRPTW`.
