# cvrpbcp (exactsolv) — working rules

Branch-cut-and-price solver for CVRP/CVRPTW built on bplex. Plan and status:
`docs/PLAN.md`; bplex API contract and requests: `docs/BPLEX_API_REQUESTS.md`.

## Hard rules
- **Never modify anything under `../bplex/`.** bplex is a separate project under
  active development; we consume only `../bplex/include/bplex.h`. Build it
  out-of-tree into `build/bplex` (MKL off). bplex bugs are reported with an
  MPS dump (`bplex_write_mps`), never patched here.
- Every optimum we report is verified independently (`src/verify`): coverage,
  capacity, time windows, recomputed cost, and comparison with the published
  value. A disagreement is a bug until explained.
- Code quality and solver performance come first; reuse the scripts and
  skills below instead of re-deriving procedures.

## Build and run (Linux/WSL is the platform)
Use the project skill `cvrpbcp-dev` (`.claude/skills/cvrpbcp-dev/SKILL.md`)
for the exact recipes. Short form, from the Windows Bash tool:

    bash scripts/wsl.sh scripts/build_wsl.sh            # bplex (once) + cvrpbcp Release
    bash scripts/wsl.sh build/cvrpbcp-Release/cvrpbcp selftest
    bash scripts/wsl.sh build/cvrpbcp-Release/cvrpbcp solve instances/cvrplib/A/A-n32-k5.vrp
    bash scripts/wsl.sh scripts/bench.sh instances/cvrplib/A 60 600   # set, max n, time limit

`scripts/wsl.sh` runs its arguments inside WSL from the repo root with the
right PATH and path-conversion settings; pass a script or a command line.

## Layout
- `src/instance` readers (CVRPLIB CVRP, PyVRP-style CVRPTW), `src/master`
  bplex master, `src/pricing` ng-route labeling, `src/cuts` CVRPSEP wrapper,
  `src/tree` CG loop and tree search, `src/verify` checker + brute force,
  `src/app/main.cpp` CLI (`selftest | info | check | root | solve`).
- `third_party/cvrpsep` vendored (CPL 1.0), built by our CMake.
- `instances/cvrplib/{A,B,E,F,P}` with `.sol`; other sets live in the
  AIRouting tree (see `instances/cvrplib/README.md`).
- `results/` benchmark outputs (`results.jsonl`, `usable.csv`).

## Conventions
- Costs are integers in cost units (`Instance::scale`: 1 for CVRP, 10 for
  CVRPTW one-decimal truncation). Depot is vertex 0, customers 1..n.
- Master rows: customers, fleet `[klo,khi]`, then robust rows (RCC and
  edge-branching rows, inactive = free bounds). Pricing uses the symmetric
  dual split; route reduced cost = `c_r - sum pi - pi_0 - row terms`.
- Fleet cap `kmax = ffd + 1` by default; a result with the cap binding must
  be re-run with a larger `--K` before it is called proven.
