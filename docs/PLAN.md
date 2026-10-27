# cvrpbcp: a branch-cut-and-price solver for CVRP and CVRPTW on bplex

Status (2026-10-27): milestones 0-2 built and passing; see `docs/STATUS.md`
for the running log of results and design decisions. bplex delivered the
incremental API on 2026-10-26 (`BPLEX_API_REQUESTS.md`, status section).
Working rules and build/run recipes: `CLAUDE.md` and the `cvrpbcp-dev` skill.

## 1. Why this project

The label pipeline (`exactlabel/EXACT_LABEL_GENERATION.md`) needs *proven*
optima for CVRP instances with many short routes (n = 150-400, k = 30-150).
The only exact backend available here is the CLP-only BaPCod binary shipped
with VRPSolverEasy. It closes bounds but has no MIP solver behind it, so the
restricted-master heuristic and the enumeration-to-MIP finish never run and it
rarely finds an incumbent on n >= 150. BaPCod source is not available, so the
decision is to build our own branch-cut-and-price (BCP) code, `cvrpbcp`, that
uses bplex as

* the LP engine for the column-generation master, and
* the MIP engine for primal heuristics (restricted master MIP, enumeration
  MIP), which is what produces incumbents.

Reference design: VRPSolver (Pessoa, Sadykov, Uchoa, Vanderbeck, *A generic
exact solver for vehicle routing and related problems*, Math. Prog. 2020) and
the RouteOpt / VRPSolverEasy feature set. We re-implement the parts needed for
CVRP and CVRPTW, not the generic modelling layer.

## 2. Ground rules

* **bplex is read-only.** Nothing under `D:\ARTC\AIRouting\bplex\` is edited; it
  is under active development. cvrpbcp links the public C API
  (`bplex/include/bplex.h`) from an out-of-tree build
  (`cmake -S ../bplex -B build/bplex -DBPLEX_USE_MKL=OFF`). API changes we
  need are *requested* in `BPLEX_API_REQUESTS.md`, never patched in.
* All code, scripts, results and docs live in `D:\ARTC\AIRouting\exactsolv\`.
* Linux first: WSL Ubuntu 22.04 (gcc 11, 16 threads, 7 GB RAM) is the build
  and run box. A Windows build is a later convenience, not a goal.
* Every reported optimum is verified independently of the solver before it is
  called proven (section 6).

## 3. Problem definitions

**CVRP** (CVRPLIB format): depot 0, customers 1..n with integer coordinates and
demands, capacity Q, cost = integer-rounded Euclidean distance per edge
(`nint(sqrt(dx^2+dy^2))`), objective = total distance, fleet size free unless
the instance fixes it. The `k` in the instance name is a lower bound only; we
add a `<= K` constraint with K = bin-packing bound + slack and re-solve with a
larger K if the cap is binding (the trap described in the handover doc,
section 5.3).

**CVRPTW** (Solomon format): as above plus service times and time windows,
unlimited homogeneous fleet, objective = total distance. The exact-method
convention truncates distances to one decimal; the reader makes rounding a
parameter so published optima can be reproduced.

## 4. Architecture

```
exactsolv/
  src/
    instance/      readers (CVRPLIB, Solomon), distance matrices, rounding
    master/        set-partitioning master on bplex; column pool; basis handling
    pricing/       RCSP labeling: resources (load, time), ng-sets, bidirectional,
                   bucket graph, dominance, completion bounds, arc fixing
    cuts/          rounded capacity cuts (CVRPSEP, vendored), lm-SRC (3-row),
                   their effect on pricing (dual contributions, SRC states)
    tree/          node, branching (edge flow), strong branching, best-bound search,
                   node LP management
    heur/          restricted master MIP (bplex MIP), diving heuristic,
                   enumeration + MIP, optional external UB (HGS-CVRP)
    verify/        independent solution checker (coverage, capacity, time windows,
                   recomputed cost)
    app/           cvrpbcp CLI, results.jsonl / usable.csv writers
  third_party/     CVRPSEP (Lysgaard), nothing else
  tests/           unit tests + instance regression lists with published optima
  docs/            this plan, bplex requests, design notes, results
  build/           out-of-tree builds (bplex, cvrpbcp)
```

Language C++17, CMake, single executable `cvrpbcp` plus a static library so
the pieces can be unit-tested.

### 4.1 Master problem

Set partitioning over routes r with cost c_r and a_ir = number of visits of
customer i on r:

```
min  sum_r c_r lambda_r
s.t. sum_r a_ir lambda_r  = 1     for each customer i   (duals pi_i)
     sum_r lambda_r      <= K     fleet cap             (dual mu <= 0)
     lambda_r >= 0
```

plus cut rows (RCC, SRC). Artificial columns (one per partitioning row, cost =
big-M derived from the instance) keep the master feasible at every node. Column
generation uses an ng-route relaxation (ng-set size 8 by default, enlarged on
cycles), Wentges dual smoothing with automatic alpha, and a Lagrangian bound
each iteration to stop early when the bound cannot improve (`ceil` on integer
costs).

bplex use: one persistent `bplex_t` per node evaluator. The master is loaded
once with the artificial columns; every CG iteration appends the priced
columns with `bplex_add_cols` and re-solves warm (dual simplex, basis, LU and
DSE weights kept by bplex). Cuts go in with `bplex_add_rows`, branching bounds
with `bplex_chg_col_bounds`/`bplex_chg_row_bounds`, strong branching and dives
with `bplex_snap_save/restore`. Column pool: nonbasic columns with large
reduced cost for several iterations are removed with `bplex_del_cols` (only
nonbasic ones, so the factorization survives) but kept in the pool for the
heuristics.

### 4.2 Pricing: resource-constrained shortest path

Generic label `(vertex, cost, resources[], ng-memory bitset, SRC states[])`.
Resources: load for CVRP; load and time for CVRPTW. Stages:

1. Heuristic pricing on a reduced graph (k nearest neighbours per vertex),
   relaxed dominance; stop when enough negative-reduced-cost routes exist.
2. Exact pricing: bidirectional labeling (forward to half the main resource,
   backward from the sink), ng-dominance, completion bounds from a backward
   pass, concatenation with hashing on ng-memory.
3. Bucket graph (resource discretised into buckets per vertex), bucket-arc
   elimination by reduced-cost fixing once a primal bound exists, enumeration
   of all routes with reduced cost below the gap when the gap is small
   (feeds the enumeration MIP).

Milestone 1 ships stages 1-2 with a mono-directional fallback; stage 3 is
milestone 4.

### 4.3 Cuts

* Rounded capacity cuts separated by CVRPSEP on the edge-flow solution
  x_e = sum_r b_er lambda_r. They are robust cuts: their duals enter the arc
  reduced costs only.
* Limited-memory subset-row cuts on triples of customers (Pecin et al. 2017),
  non-robust: each adds a state to the labels; limited memory keeps the damage
  to dominance small. Milestone 4.
* Cut management: add the most violated, remove non-binding cuts after a few
  rounds, keep a pool.

### 4.4 Branching and tree search

* Branching on edge flows x_e (fractional, closest to 0.5 first), evaluated by
  strong branching on a shortlist (phase 1: LP only, phase 2: a few CG
  iterations) with pseudocosts from the start.
* Fleet-size branching when sum lambda is fractional.
* Best-bound node selection with plunging; node LP restored from a stored basis.
* Integrality test on lambda and on the edge flows, with the solution rebuilt
  from the routes and verified by `verify/`.

### 4.5 Primal side (where bplex's MIP earns its place)

* Restricted master MIP: the pool columns of the current node with integrality
  on lambda, solved by bplex (`relax=0`) with a time cap, a cutoff at the
  incumbent and feasibility emphasis. Run at the root after convergence and
  every N nodes, or when the pool changed enough.
* Diving heuristic with limited discrepancy: fix the lambda closest to 1,
  re-run CG with the reduced graph, backtrack a few times.
* Enumeration + MIP: once gap < threshold and the number of routes with
  reduced cost below the gap is below a limit (e.g. 10^5), enumerate them and
  finish the node with one bplex MIP (set partitioning on enumerated routes).
* Optional external upper bound: `--ub` or `--hgs` (HGS-CVRP binaries are in
  this repository) to start with a tight cutoff. Not required for correctness.

## 5. Milestones and gates

| # | Scope | Gate (must pass before the next milestone) |
|---|---|---|
| 0 | WSL setup: out-of-tree bplex build (MKL off), CMake skeleton, CVRPSEP vendored, CVRPLIB A/B/E/F/P/X and Solomon sets downloaded with their published optima | `cvrpbcp --selftest` solves a 5-customer instance and matches brute force |
| 1 | CVRP root: readers, master on bplex with incremental columns and warm solves, labeling pricing (heuristic + exact, ng-routes), smoothing, Lagrangian bound | root bounds of A-n32-k5 and P-n55-k7 equal the known ng-route + RCC root bounds within rounding; the dual simplex re-solves in few iterations after a reload |
| 2 | Exact: RCC via CVRPSEP, edge branching with strong branching, best-bound tree, incumbent verification | all A, E, P instances with n <= 60 proven to their published optima |
| 3 | Incumbents: restricted master MIP on bplex, diving, external UB option, results.jsonl / usable.csv writers | P-n55-k7 and X-n101-k25 have an incumbent right after the root; proven costs equal the published optima; a first n = 150, k >= 35 XML instance is proven |
| 4 | Speed: bidirectional labeling, bucket graph, reduced-cost arc fixing, lm-SRC, enumeration + MIP | B, F, X up to n = 200 and a sample of XML100 proven within the pipeline's time limits; time per instance reported |
| 5 | CVRPTW: time resource, Solomon reader and conventions | Solomon 25- and 50-customer instances proven to published optima, then the easy 100-customer ones |

Milestone 3 is the point where the project pays back: proven large-k labels
start flowing to `data/<set>/` in the shape the pipeline ingests.

Progress (2026-10-27, details in `STATUS.md`):

* Milestones 0, 1, 2: done. 64/64 CVRPLIB A/B/E/F/P instances with n <= 60
  proven within 600 s and equal to the published optima (exact fleet).
* Milestone 3: done except the "incumbent right after the root" wording:
  the root RMIP rarely helps; incumbents come from the Clarke-Wright seed,
  node LPs during the plunge and the enumeration MIP. X-n101-k25 proven in
  13-20 s, XML150_1121_01 (n = 150, k = 37) in 100 s. `results.jsonl` via
  `--jsonl`, `usable.csv` via `scripts/summarize.py`.
* Milestone 4: bidirectional labeling, lm-SRC and enumeration + MIP are in;
  a q-route completion bound prunes labels and arc extensions in every
  pricing stage (2026-10-27), cut rows go in one batch per round, LP
  columns are purged inside the node, and the pricing stages are tuned at
  n = 150 (root of X-n153-k22 92 -> 9 s). Bucket graph, reduced-cost arc
  fixing and proper dual stabilization are open; the n = 200 tree is being
  re-measured (`results/large_v7`).
* Milestone 5 (CVRPTW): time resource in the labeling and enumeration
  (mono-directional exact pass, time-free enumeration bound). Solomon
  100-customer C101, C102, C201, R101, RC101, R201 proven to the published
  optima (0.1-123 s); the C1/R1/RC1 series benchmark is in
  `results/solomon/`. 25/50-customer Solomon files are not in the tree.

## 6. Verification policy

An instance is *proven* only if all of these hold:

1. the global lower bound satisfies `ceil(lb - 1e-6) >= cost` (integer costs);
2. the solution passes the independent checker (every customer exactly once,
   capacity and time windows respected, indices in range);
3. the cost recomputed from coordinates by the checker equals the reported cost
   and, for benchmark instances, the published optimum;
4. the fleet cap was not binding (otherwise re-solve with a larger cap).

Benchmark sets: CVRPLIB A, B, E, F, P, X (n <= 200) for CVRP; Solomon for
CVRPTW. Each run writes `results.jsonl` (every attempt, with status, time,
nodes, root bound) and `usable.csv` (proven only).

Any disagreement between cvrpbcp and a published optimum is a bug until
explained. If the cause is a wrong LP or MIP answer from bplex, the model is
dumped as MPS and reported to the bplex side with the expected result; we do
not patch bplex.

## 7. Risks

* **Pricing speed** dominates everything above n = 100; milestone 4 exists for
  that reason and is where most of the effort goes.
* **Master LP engine**: bplex's incremental API is new (2026-10-26) and
  exercised so far by its own `api_incr` test; cvrpbcp's master keeps a KKT
  self-check (debug builds) and dumps the model with `bplex_write_mps` on any
  disagreement, reported to the bplex side rather than patched.
* **Numerics**: integer data keeps LPs well scaled; tolerances are set
  explicitly (primal/dual 1e-7, reduced-cost acceptance 1e-6, integrality
  1e-6) and the Lagrangian bound is used for the proof, not the LP value.
* **Memory**: 7 GB on WSL bounds label counts and parallel workers; the bucket
  graph and label pruning are sized with that in mind.

## 8. Deliverables

* `cvrpbcp` binary and library, CMake build, WSL build script.
* `docs/`: this plan, `BPLEX_API_REQUESTS.md`, design notes per component,
  results tables per milestone.
* `tests/`: unit tests, regression lists with published optima, a one-command
  benchmark runner producing `results.jsonl` / `usable.csv`.
