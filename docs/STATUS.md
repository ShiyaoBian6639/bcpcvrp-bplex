# cvrpbcp status log

Newest first. Numbers are from WSL Ubuntu 22.04, gcc 11 -O3, single thread,
bplex built with MKL off.

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
