# bcpcvrp-bplex (exactsolv)

An exact branch-cut-and-price solver for the capacitated vehicle routing problem
(CVRP) and its time-window variant (CVRPTW), written in C++17 on top of the
[bplex](https://github.com/ShiyaoBian6639/bplex) LP/MIP engine.

bplex is consumed only through its public header (`bplex.h`); the two projects
are developed independently and communicate through the API documented in
`docs/BPLEX_API_REQUESTS.md`. Nothing under the bplex tree is modified here.

## What is inside

* `src/instance`  readers for CVRPLIB CVRP and PyVRP-style CVRPTW files
* `src/master`    set-partitioning master on a persistent bplex engine
                  (incremental columns and rows, column pool and purge,
                  robust and non-robust cut rows, restricted MIP)
* `src/pricing`   ng-route labeling (bidirectional for CVRP, time windows,
                  limited-memory subset-row cut states, route enumeration)
* `src/cuts`      rounded capacity cuts (CVRPSEP) and limited-memory
                  subset-row cuts
* `src/tree`      column generation loop and best-bound/plunging tree search
                  with pseudocost edge branching, diving and enumeration + MIP
* `src/heur`      Clarke-Wright seed and the enumeration MIP
* `src/verify`    independent solution checker and brute-force oracles
* `src/app`       the `cvrpbcp` command line (`selftest | info | check | root | solve`)
* `third_party/cvrpsep` vendored CVRPSEP (CPL 1.0)
* `instances/cvrplib/{A,B,E,F,P}` with published solutions
* `scripts/` build, benchmark, regression and summary helpers (WSL/Linux)
* `docs/` plan, status log and the bplex API contract

## Build (Linux or WSL)

bplex (API 0.2 or later) is expected as a sibling checkout (`../bplex`); it
is built out-of-tree, installed to `build/bplex-install` and consumed as the
CMake package `bplex::bplex` (`find_package(bplex 0.2 CONFIG)`):

    scripts/build_wsl.sh                 # cmake+ninja: build/bplex (+install), build/cvrpbcp-Release
    build/cvrpbcp-Release/cvrpbcp selftest
    build/cvrpbcp-Release/cvrpbcp solve instances/cvrplib/A/A-n32-k5.vrp --time 600
    scripts/regress.sh 300               # regression set with known optima

From a Windows Claude Code session use `bash scripts/wsl.sh <command>`.

## Verification policy

Every optimum the solver reports is re-checked independently (coverage,
capacity, time windows, recomputed cost) and compared with the published value.
Benchmarks are run with the known solution loaded (`CVRPBCP_DEBUGSOL`) so a lost
optimal solution is reported, not silently pruned. See `docs/STATUS.md` for
the current results and `docs/PLAN.md` for the milestone plan.
