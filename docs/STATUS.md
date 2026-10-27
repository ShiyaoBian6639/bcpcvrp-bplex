# cvrpbcp status log

Newest first. Numbers are from WSL Ubuntu 22.04, gcc 11 -O3, single thread,
bplex built with MKL off.

## 2026-10-28: duplicate-column bug, node-level tuning on a medium set

* **Bug (ours): a route listed twice in one `Master::add_routes` batch got
  two LP columns.** The dedupe only looked at routes already in the pool;
  within a batch the second copy of a new route was "found" with
  `col_of_ = -1` and queued again. The tree reads the pool-indexed lambda,
  which keeps only one column's value, so the fleet sum was under-counted:
  on XML100_3313_15 (kmin 12) the root reported 11.089 vehicles, branched
  on the fleet, the "<= 11" child was infeasible and the ">= 12" child was
  the same LP, forever (8000+ nodes, then the per-node enumeration sets ran
  the process out of memory). Edge flows and the integrality test read the
  same lambda. Fix: a `-2` "queued in this batch" mark in `col_of_`; the
  master now throws when a pool route has two LP columns (cheap invariant
  per solve) or when the fleet row's activity leaves its bounds; a fleet
  child with an empty vehicle range is dropped instead of being pushed.
  XML100_3313_15: optimal 22311 in 9 s. A standalone check of
  `bplex_chg_row_bounds` on ranged rows (`tests/bplex_rowbounds.c`) passes:
  bplex was not at fault.
* Node-level changes, measured on the 12-instance XML100 medium set
  (`instances/xml100/medium`, `scripts/pbench.sh`, 300 s, four at a time):
  - CVRPSEP's "existing cuts" list now holds only the RCC rows with a
    nonzero dual at the last solve; duplicates of the others are filtered by
    a set hash. Separation time per call in the tree dropped from ~64 ms to
    ~20 ms on XML200; medium set 1768 -> 1327 s (same 8 solved).
  - lazy pool pricing (a pool route with clearly positive reduced cost is
    skipped for 2/6/15 scans), adaptive heuristic stages (a stage empty
    twice in a row is skipped for the rest of the node, `--stage-skip`),
    CG-iteration and separation counters in the `split:` line.
  - row cleanup every 5 nodes instead of 25 (`--cleanup-every`): medium set
    1327 -> 1214 s and 1344 -> 1255 s in two pairs; adopted.
  - `--stage-skip 1` was worse (1964 s) and hit the LP time limit at the end
    of two runs (reported as "error"; now any node stopped by the clock ends
    the run as `time_limit`). Skipping after 2 empty results vs never:
    medium set 1255 vs 1229 s, CVRPLIB n <= 60 set 321 vs 263 s (both
    64/64 proven, zero debug reports) -> stage skipping is off by default
    (`stage_skip_after` 100).
* Verified state of the defaults (knn 5, 75 cols, purge every 10, q-bound
  all stages, RCC existing = binding rows, lazy pool, cleanup every 5, no
  stage skipping): `results/v9b` 64/64 CVRPLIB n <= 60 proven = published,
  total 263 s, slowest E-n30-k3 40 s; medium set 9/12 in 300 s (class-6
  instances with ~20-customer routes stay open: pricing, not the tree).
  Regression passes.

## 2026-10-27 (night): tree bug fixed, large instances re-measured

* **Bug (ours, not bplex):** `cleanup_rows` marked the branch rows that open
  nodes still reference in an array sized "rows + 4096" but indexed by row
  id; ids are never reused, so after a few hundred nodes an inactive branch
  row of an open node was deleted, and the next switch to that node asked
  bplex to change the bounds of a stale row index (`BPLEX_ERR_ARG`, seen as
  `bplex_chg_row_bounds failed` after 390-455 s on XML200_2251_68 and
  X-n153-k22). The bug predates today; it only surfaced because nodes are
  now fast enough to reach that many ids within the time limit. Fix: marks
  indexed by `Master::row_id_bound()`, `set_row_active` rejects deleted
  ids, failed row edits print rc/row and dump the model, bplex's warnings
  and errors (level 0) reach our log via `bplex_set_log_callback`.
  Regression passes. Earlier 900 s results on these instances were not
  wrong (a relaxed node keeps a valid bound) but the search degraded.
* 900 s, free fleet, with the fix (`results/large_v7/*.fix.log`):

  | instance | lb | ub | nodes | note |
  |---|---|---|---|---|
  | X-n101-k25 | 27591 | 27591 | 9 | proven in 8.4 s (published 27591) |
  | X-n153-k22 | 21136 | 21313 | 268 | published 21220; before today lb 21112 / ub 21415 |
  | XML150_3372_54 | 41702 | 42339 | 126 | before: 41699 / 42821 |
  | XML200_2251_68 | 31355 | 31461 | 1023 | before: 31338 / 31484 |

  Time split X-n153-k22: pricing 423 s, master 253 s, dive 139 s, enum 26 s;
  XML200: pricing 266 s, master 279 s, cuts 222 s, enum 100 s (115 tries,
  18 MIPs, 18 nodes closed). Node throughput, not the root, is now the
  limit; a per-stage split (`split:` line after `times:`) was added and the
  two instances are being re-run to pick the next lever.

## 2026-10-27 (later): repository, q-route completion bound, bplex 0.2

* Code pushed to `github.com/ShiyaoBian6639/bcpcvrp-bplex` (branch `main`;
  cvrpsep vendored as plain files, `build/`, binaries and large MPS dumps
  ignored, `results/solomon/set1` symlinks ignored).
* **q-route completion bound in pricing** (`Pricer::compute_qbound`): the
  cheapest 2-cycle-free walk depot -> v with load <= L, O(Q n^2) per dual
  vector, SRC ignored. Labels and arc extensions with rc + bound > 0 are
  dropped in every stage (`--qbound s` sets the first stage, 3 = off); the
  enumeration's ng bounding pass prunes with threshold gap + vpay[v]. Root
  relaxation, free fleet, same bounds in all runs:

  | instance | before | bound from stage 1 | bound in all stages (default) |
  |---|---|---|---|
  | X-n153-k22 | 92 s (pricing 46) | 58 s (pricing 13) | 48 s (pricing 11) |
  | XML150_3372_54 | 76 s (pricing 54) | 42 s (pricing 18) | 41 s (pricing 16) |

  `results/v6` (DEBUGSOL=1): 64/64 A/B/E/F/P n <= 60 proven = published,
  zero debug reports, total 427 s (v5: 454 s; single instances swing both
  ways with the tree, P-n55-k10 40 -> 96 s, P-n55-k7 83 -> 40 s).
  Regression passes. The remaining root time at n = 150 is master LP
  (~45 ms per solve, ~17 s) and cut separation (~20 s); the ROOT line now
  prints `cuts` and `src` times.
* **bplex API 0.2** (round 2 delivered: package, `bplex_version()`, pool
  getters, log callback, bulk bounds, error contract). cvrpbcp builds with
  `find_package(bplex 0.2 CONFIG)` against `cmake --install` output
  (`build/bplex-install`), no reference to the bplex source tree remains;
  `info` prints the version and every JSONL line carries `"bplex"`.
  Regression and selftest pass on the new library.
* **Heuristic-pricing / master tuning at n = 150** (root, free fleet, two
  runs at a time; `results/pr/s2b,s3,s4`). A grid run concurrently with the
  bench was too noisy (wall times 2x) and was discarded. Findings:
  - sparse stage with 5 nearest successors instead of 10 (`--knn 5`): fewer
    labels and fewer CG iterations; 3 is too sparse (pricing time doubles),
    20 is worse than 10;
  - 75 columns per iteration instead of 150 (`--cols 75`): smaller LP, same
    iteration count; 50 and 100 are no better;
  - label cap 20 vs 40: neutral;
  - **cut rows were added one bplex edit at a time**: at 30k LP columns the
    per-row edits cost 20 s of the 47 s root on X-n153-k22 ("cuts" time);
    now one `add_rows` batch per round: 1.7 s (CVRPSEP itself 0.3 s);
  - purging LP columns every 10 CG iterations inside the node
    (`--purge-every 10`, age 6, floor 20 per customer) halves the master
    time again; every 5 is the same.

  Root time with the new defaults (knn 5, cols 75, purge every 10, batched
  rows): X-n153-k22 92 -> 9 s, XML150_3372_54 76 -> 15 s. The root bound
  varies by a few units between configurations because the SRC rounds and
  the tail-off see different LP solutions.
* `results/v7` (new defaults, DEBUGSOL=1): **64/64 proven = published, zero
  debug reports, total 315 s (v6 427 s, v5 454 s), slowest 35 s** (P-n50-k8;
  P-n55-k10 96 -> 34 s, F-n45-k4 33 -> 8 s, E-n30-k3 10 -> 25 s).
  Solomon (300 s, five concurrent): C101, R101, R105, RC101 proven as before,
  **RC102 now proven (14574, 262 s; open before)**, R104 open (lb 9653 /
  ub 9859, published 9715). Regression passes.

## 2026-10-28: verified benchmark, CVRPTW (milestone 5 started)

* `results/v5` (fixed enumeration, dive, savings, cut aging; run with
  `DEBUGSOL=1`): **64/64 A/B/E/F/P instances with n <= 60 proven, equal to
  the published optima, zero debug reports.** Slowest P-n55-k7 83 s,
  P-n50-k8 45 s, P-n55-k10 40 s, A-n60-k9 27 s, B-n50-k8 27 s, F-n45-k4 26 s.
* Time windows in the labeling: label time = end of service, extension
  waits for the window and rejects late arrivals, dominance adds time,
  return to the depot checked at completion. The exact pass stays
  mono-directional for CVRPTW (the bidirectional join needs symmetric
  halves); the enumeration bound uses a time-free pass (valid either way).
  Fleet for CVRPTW: the instance's VEHICLES is an upper bound, not exact.
  Solomon 100-customer (one-decimal truncation, cost units x10): C101
  0.1 s, C102 0.6 s, C201 9 s, R101 7.7 s, RC101 12 s, all equal to the
  published optima; R201 (wide windows, 2-4 long routes) first failed in
  300 s because each enumeration attempt took 50 s. Enumeration attempts
  now have a 10 s cap (`enum_time_cap`, checked inside the bound pass) and
  a failed attempt blocks retries anywhere until the gap halves: R201
  proven in 123 s (60 nodes).
* Solomon C1/R1/RC1 series, 300 s each, verified (`results/solomon/`):
  20 of 29 proven and equal to the published optima (all C1, R101-R103,
  R105-R107, R109, RC101, RC103, RC105, RC107), zero debug reports; the
  open ones (R104, R108, R110-R112, RC102, RC106, RC108, and RC104 which
  overran the runner timeout on the pre-budget binary) end with valid
  bounds below the published values. Open-window R1/RC1 instances need
  the bidirectional time labeling and faster pricing.
* Sub-steps (dive, RMIP, enumeration MIP) now take the remaining time
  budget, every master LP solve carries bplex's `time_limit` equal to the
  node's remaining budget (a stopped solve ends the node with an error
  instead of hanging), and stdout is line-buffered so a killed run still
  leaves its log. Background: X-n153-k22 twice produced no output in
  1000 s with the old binary; with these changes it runs the tree normally
  (root 130 s, root dive 131 s without a solution, incumbent 21415 from a
  dive at 366 s vs published 21220, bound 21112 at 420 s).
* Frontier for the pipeline: XML150_3372_54 (n = 150, k = 24, routes of 6)
  gets its first incumbent from the dive at 169 s and ends 900 s with a 2.7%
  gap; pricing takes 660 of the 900 s. Faster exact pricing (bucket graph,
  reduced-cost arc fixing) is the next lever for n >= 150 with k < n/5.

## 2026-10-28: second enumeration bug (completion bound with SRC)

* The verified benchmark (`DEBUGSOL=1`) caught A-n61-k9 proven at 1035
  instead of 1034 with routes of reduced cost 0 missing from the
  enumeration. The per-prefix trace showed the completion bound above the
  true suffix cost at one vertex although every arc cost is symmetric.
* Cause: the bound at (v, load) is the cheapest backward ng label at v, and
  that label includes v's own visit. For a subset-row cut containing v, the
  pair (x, v) can be paid by the backward label and again by the prefix
  (v paired with a prefix vertex), so the bound over-counts by up to the
  penalties v itself can incur. Fixed by subtracting, per vertex, the sum of
  |dual| over the active SRCs containing it. The bidirectional join has no
  shared vertex and needs no correction.
* The debug check was then refined: a known route counts as present when
  the set holds any order of the same customers that is no dearer (that is
  what enumeration dominance keeps). After the fix A-n61-k9, P-n50-k8 and
  B-n50-k8 run without a single report.

## 2026-10-28: column-fixing dive, plunge limit

* XML150_3372_54 and X-n153-k22 ran 900 s without any incumbent (depth
  709 and 78 plunges on edge branches, nothing integral), so enumeration
  never started. Added a column-fixing dive: fix the elementary route with
  the largest fractional lambda (column bounds [1,1] through
  `bplex_chg_col_bounds`), remove its customers from pricing, re-run CG
  without cuts, repeat until the LP is integral or infeasible; then undo
  (columns stay in the pool) and re-solve the node. Runs at the root and
  every 20 nodes while no incumbent exists (`--no-dive`, `--dive-time`).
  P-n55-k7: dive gives 597 at 3 s (savings 616). Plunge depth now capped at
  50 (`--plunge`), then best-bound.

## 2026-10-28: enumeration dominance bug (wrong "optimal" on E-n30-k3)

* After the night's changes `scripts/regress.sh` reported E-n30-k3 proven
  at 537 instead of 534. Bisection by flags pointed at enumeration, and a
  new debug hook (`CVRPBCP_DEBUGSOL=<file.sol>`: every node checks whether
  the known solution is feasible there, survives the bound, and has all its
  routes in the enumerated set) showed a route of the optimal solution with
  reduced cost 3.96 missing from a set enumerated with gap 12.15, although no
  prefix was pruned by the completion bound.
* Cause: the elementary enumeration dominated two labels with the same end
  vertex and visited set by reduced cost alone. With arc-dependent duals
  (RCC and branching rows) the order with the smaller reduced cost can be
  the dearer route; the MIP then only sees the dearer one and misses the
  improving solution. Dominance now requires cost <= and reduced cost (plus
  SRC penalty) <=. Pricing labels are unaffected (there only the minimum
  reduced cost matters).
* Also fixed on the way: the enumeration MIP's absolute gap was 0.999 with
  cutoff ub - 0.5, which lets bplex prune a node with bound in
  (ub - 1.5, ub - 1]; now 0.45.
* Lesson recorded in the skill: run `scripts/regress.sh` after every solver
  change, and run benchmarks with `DEBUGSOL=1` (bench.sh) so that a lost
  known solution is reported even when the final answer happens to be right.

## 2026-10-27 (night): pricing stages, savings seed, cut aging

* The bidirectional join was the cost on long-route instances: early
  iterations run on big-M duals, so 25-50k labels all have hugely negative
  reduced cost and nothing prunes the join (F-n45-k4 root 534 s). Now the
  heuristic stages use the cheap mono pass with a per-vertex label cap that
  keeps the best labels by reduced cost (stage 0: sparse k-nearest graph,
  stage 1: full graph, both capped at 40), and only the exact stage is
  bidirectional. F-n45-k4: 534 -> 45 s.
* Clarke-Wright savings seeds the master and gives a first incumbent at
  0 s when its route count fits the fleet (`--no-savings` to disable).
* Routes are de-duplicated in the master (a known route out of the LP is
  re-added instead of copied).
* Cut rows with a zero dual for 30 LP solves are deleted at the periodic
  row cleanup (every 25 nodes); CVRPSEP's "existing cuts" list is rebuilt
  from the master's live RCC rows at each call so a deleted cut can return.
* XML150_3372_54 (n = 150, k = 24) root: 330 CG iterations / 83 s before,
  261 / 68 s now (38 exact passes of ~1 s each remain the cost). Wentges
  smoothing with a fixed alpha (0.5, 0.8) did not converge within 250 s and
  stays off; a proper auto-alpha version is future work.
* `results/v3` (enumeration on, before this round): 63/63 of A/B/E/P
  proven and matching; F-n45-k4 hit the 660 s runner timeout (its root was
  the 534 s case above).

## 2026-10-27 (evening): full A/B/E/F/P benchmark, enumeration

* `results/v2/*.jsonl` (exact fleet, SRC, bidirectional labeling, pseudocost
  branching, no enumeration): **64/64 instances with n <= 60 proven within
  600 s and equal to the published optima.** Slowest: F-n45-k4 534 s (all at
  the root: 4 long routes of ~11 customers make the bidirectional join
  expensive; the mono version took 80 s), P-n55-k10 50 s, B-n50-k8 58 s,
  A-n60-k9 44 s; median well under 10 s.
* Route enumeration + set-partitioning MIP (bplex, fresh object) after a
  node when an incumbent exists: all elementary routes with
  rc <= ub - 1 - LB(pi) + min_rc(pi), completion bound from the mono ng pass.
  Sets above `enum_mip_max_routes` are kept and replace labeling in the
  subtree (pricing by scanning, elementary routes, stronger bounds); they
  shrink at each descendant. On XML150_1121_01 enumeration takes 0.8 s for
  215k routes, but bplex needs 16-60 s per MIP with 5k-15k columns and
  0.2-5 s below 1000 columns, so the MIP threshold is 1000 routes; a failed
  enumeration is retried in the subtree only once the gap shrank by 25%.
  Effect (seconds): XML150_1121_01 280 -> 100, P-n55-k7 37 -> 28, P-n50-k8
  19 -> 29 (12 small MIPs), E-n30-k3 18 -> 20. `--enum-dump prefix` writes
  the MIPs as MPS for the bplex side (`results/enum_mips/`).
* XML200_2251_68 (n = 200, k = 41-44) without enumeration: not closed in
  900 s (lb 31338, ub 31484, 873 nodes); the master LP took 453 s with 42k
  columns and 3700 rows, so cut-row aging is the next item.

## 2026-10-27 (later): bidirectional labeling, B/F sets, first large-k proofs

* Bidirectional labeling on the load resource (one pass to Q/2, then
  concatenation with the ng-memory disjointness test and the SRC junction
  penalty; symmetric data so one label set serves both directions). Effect
  on `scripts/regress.sh` (seconds, mono -> bidirectional): E-n30-k3 with
  exactly 3 trucks >300 -> 18, P-n50-k8 159 -> 19, X-n101-k25 39 -> 13,
  P-n55-k7 30 -> 34 (short routes: no gain). Default on; `--mono` to compare.
* B and F sets (n <= 60, exact fleet, SRC, mono labeling): 18/18 proven and
  equal to the published optima; slowest B-n50-k8 193 s, F-n45-k4 80 s.
* Free-fleet proofs for the pipeline target: X-n101-k25 27591 in 50 s
  (mono; 13 s bidirectional), XML150_1121_01 (n = 150, k = 37) 43972 in
  199 s (mono), 270 nodes, root 9 s. The exactlabel record for a differently
  calibrated copy of that instance (capacity 26) ended at a BaPCod/CLP time
  limit with no solution.
* Milestone 3 gate is therefore met except that incumbents come from node
  LPs during the plunge rather than from the root RMIP (the 5 s RMIP found
  576 on P-n55-k7, nothing on X-n101-k25): diving or a longer RMIP would help.

## 2026-10-27: milestone 2 gate reached (A/E/P, n <= 60)

Benchmark `scripts/bench.sh instances/cvrplib/{E,A,P} 60 600` before SRC
cuts and before the exact-fleet convention: 46/46 proven optimal within 600 s
(max 274 s, P-n55-k10), 41 equal to the published optimum. The 5 others
(P-n22-k8, E-n30-k3, P-n50-k8, P-n55-k15, P-n55-k8) were solved with a free
fleet and found cheaper solutions with a different vehicle count; CVRPLIB's
A/B/E/F/P optima assume exactly the k vehicles of the COMMENT line. The
reader now parses "No of trucks" and `solve` uses that count exactly by
default (`--fleet max` for at most, `--free-fleet` to ignore it; the label
pipeline wants `--free-fleet`).

Changes that mattered, measured on P-n55-k7 (optimum 568):

| version | root bound | nodes | time |
|---|---|---|---|
| plain tree, most-fractional branching, 10 s RMIP every 20 nodes | 558.9 | 448 (not closed) | > 600 s |
| + pseudocost branching, column purge + pool pricing, row artificials | 558.9 | 505 | 107 s |
| + limited-memory SRC on triples (<= 100 active) | 563.1 | 68 | 35 s |

Findings:
* Strong branching by probing (bplex `add_rows` + dual simplex capped at 100
  iterations + `del_rows`) cost 40-90 ms per probe and lost to pseudocosts
  (166 s vs 107 s). Probe rows are now added in one batch as inactive rows and
  switched by bound changes; still off by default (`--sb k` enables).
* LP size dominates node time: purging columns nonbasic with rc > 0 for 6
  solves (floor 20 columns per customer) and re-adding from the pool when
  their reduced cost turns negative cut master time by ~4x despite heavy churn.
* Every ">=" row carries a big-M artificial column so a node LP over the
  current columns is never infeasible (otherwise edge "up" rows made the LP
  infeasible and a wrong `infeasible` prune followed). bigM = 4(n+1)max d(0,i).
* SRC sign: the row is sum coef lambda <= 1 with dual y <= 0, so a route's
  reduced cost *increases* by -y per coefficient unit; the first version had
  the sign reversed and CG looped forever.

## 2026-10-26: milestones 0 and 1

* Milestone 0: CMake skeleton, bplex out-of-tree build, CVRPSEP vendored,
  readers (CVRPLIB CVRP incl. EXPLICIT matrices, PyVRP-style CVRPTW),
  independent checker, brute force (subset DP) and bplex set-partitioning MIP
  oracle. `selftest` agrees on random 5-8 customer instances; all 90 CVRPLIB
  A/B/E/F/P solutions, Solomon C101/R101, X-n101-k25 and XML100_1111_01
  verify to their published costs.
* Milestone 1: master on bplex (incremental columns/rows, warm solves),
  ng-route labeling (ng 8, load buckets), RCC via CVRPSEP. Root: A-n32-k5
  784 (integral), E-n13-k4 247 (integral), P-n55-k7 558.9 in 0.7 s.
