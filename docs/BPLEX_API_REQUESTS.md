# Requests to bplex from cvrpbcp (branch-cut-and-price on bplex)

Written 2026-10-26 by the cvrpbcp side. bplex is under active development and
is not modified from this project; this document is the list of public-API
additions and behaviour guarantees cvrpbcp needs, in priority order, with the
workaround used until each one exists. Where an internal bplex function already
does the job, it is named so the public wrapper is cheap to add.

Everything refers to `include/bplex.h` and the `bplex_t` object. The LP is
always a minimisation, rows are ranged (`rowlo <= a x <= rowup`), and columns
are compressed sparse column (CSC) as in `bplex_load`.

## Status 2026-10-26: delivered

The bplex side merged the incremental API (`apiedit` branch, commit 7690f7f,
`tests/api_incr.c`). Checked from this project on the same day:

* `include/bplex.h` now has every P0 item (`bplex_add_cols`, `bplex_add_rows`,
  `bplex_del_cols`, `bplex_del_rows`, `bplex_chg_col_bounds`,
  `bplex_chg_row_bounds`, `bplex_chg_obj`, `bplex_chg_col_type`, documented
  warm-start contract on a persistent engine) and every P1 item (MIP solves on
  a model copy that leave the LP engine alone, `bplex_set_mip_start`,
  `bplex_set_incumbent_callback`, cutoff with a valid `bplex_bestbound` on a
  stopped warm solve, `bplex_snap_save/restore/free`, `bplex_get_farkas`,
  `bplex_write_mps`). P2: the header documents row-wise pricing on column-rich
  LPs, no sifting on warm solves, per-object thread safety, a `"work"` info
  counter, the sign convention `dj = c - A^T y` and the infinity rule. The
  "behaviour to confirm" list is answered in the header comment.
* Out-of-tree build in WSL with `-DBPLEX_USE_MKL=OFF` into
  `exactsolv/build/bplex` compiles (gcc 11, Ninja); `api_incr` passes
  628401/628401 checks (2000 random edit/warm-solve sequences, 2.4 iterations
  per solve, snapshots, Farkas, threads).
* `api_incr cg` (m=300, pool 20000, batch 50): incremental edits and
  reload+basis both converge in ~37 solves at ~8.5 ms per solve; a
  0-iteration re-solve costs 0.04 ms incremental vs 0.6 ms reload vs 54 ms
  cold. The workaround table below is therefore historical; cvrpbcp uses the
  incremental calls from milestone 1 on.

Open items for later: none blocking. Minor notes to raise when hit:
`bplex_del_cols` resets DSE weights when a basic column is removed (column
pool cleaning will delete nonbasic columns only); `bplex_add_cols` drops column
names (cvrpbcp keeps its own route index).

## How a branch-cut-and-price code drives an LP solver

Per branch-and-bound node, column generation runs 20-500 iterations. Each
iteration: solve the master LP (dual simplex, warm), read row duals, price,
**add 1-100 columns**, repeat. Every few iterations: **add cut rows**, **remove**
idle columns and non-binding rows. At node boundaries: **change column bounds**
(branching, fixing), **save and restore a basis**, run **strong branching**
(many short re-solves from the same basis). Occasionally: solve the current
master **as a MIP** with integrality on the columns, a cutoff and a time cap, to
find incumbents. The master typically has 100-1,000 rows and 1,000-100,000
columns, integer-valued data, and is heavily degenerate.

The LP solve itself is cheap when warm; the cost model is dominated by how much
state survives between solves. That is what the P0 requests are about.

## P0: incremental editing of a loaded LP with a surviving basis

Needed from milestone 1 (root column generation). Until then cvrpbcp reloads
the whole model each iteration (`bplex_load` + `bplex_set_basis` with the old
basis extended), which is O(nnz) per iteration, re-runs scaling, loses the dual
steepest-edge weights and the LU factors, and skips presolve.

### P0-1 `bplex_add_cols`

```c
/* Append ncols columns. Returns 0 on success.
   New columns become nonbasic at their lower bound (upper if lower is -inf,
   free at zero if both infinite). The current basis and factorization stay
   valid; the next bplex_solve warm-starts from them. */
int bplex_add_cols(bplex_t *bp, int ncols, const double *obj,
                   const double *collo, const double *colup,
                   const int *colstart, const int *rowidx, const double *values);
```

Internal status: nothing equivalent in `spx_api.c` today (`spx_rebuild_csc`
exists for rebuilding the column copy; the MIP never adds columns).

### P0-2 `bplex_add_rows`

```c
/* Append nrows rows. New rows enter the basis through their slack
   (dual simplex continues without refactorization). Returns 0 on success. */
int bplex_add_rows(bplex_t *bp, int nrows, const double *rowlo, const double *rowup,
                   const int *rowstart, const int *colidx, const double *values);
```

Internal status: `bp_spx_add_rows` exists (used by the MIP for cuts); this is a
public wrapper plus bookkeeping in `bplex_t` (model arrays, basis length).

### P0-3 `bplex_del_cols`, `bplex_del_rows`

```c
/* Remove the listed columns / rows (0-based, any order, no duplicates).
   Remaining columns / rows keep their relative order; indices shift down.
   Deleting a basic column or row forces a refactorization; deleting nonbasic
   columns and basic rows keeps the factorization. Returns 0 on success. */
int bplex_del_cols(bplex_t *bp, int ndel, const int *idx);
int bplex_del_rows(bplex_t *bp, int ndel, const int *idx);
```

Internal status: `bp_spx_del_rows(s, mask)` exists; column deletion does not.
cvrpbcp will only delete nonbasic columns and non-binding (basic-slack) rows,
so the simple case is the one that matters; rejecting the other case with an
error code is acceptable.

### P0-4 Bound, cost and right-hand-side changes

```c
int bplex_chg_col_bounds(bplex_t *bp, int j, double lo, double up);
int bplex_chg_row_bounds(bplex_t *bp, int i, double lo, double up);
int bplex_chg_obj(bplex_t *bp, int j, double cost);
```

Internal status: `bp_spx_set_col_bounds`, `bp_spx_set_row_bounds`,
`bp_spx_set_cost` exist (the MIP uses them between nodes). Semantics needed:
a bound change on a nonbasic column moves the column to the new bound and
updates primal values; the basis stays; the next solve is a warm dual simplex.

### P0-5 `bplex_solve` warm-start contract

After any of P0-1 to P0-4, `bplex_solve` must:

* start from the current basis (no crash basis, no presolve unless its
  postsolve returns a basis and the dual steepest-edge weights survive);
* honour `method=1` (dual) and `method=2` (primal) as given, and never race
  threads when `threads=1`;
* leave the object in a state where `bplex_get_solution`, `bplex_get_basis`,
  `bplex_objval` are valid, including after `BPLEX_ITER_LIMIT`,
  `BPLEX_TIME_LIMIT` and `BPLEX_CUTOFF` stops (section P1-4).

Today `bplex_load` discards everything, and `bplex_set_basis` only stores a
start basis that is consumed by the next solve. A cheap intermediate step that
already helps: keep the LU factors and DSE weights when `bplex_set_basis` is
given the basis the object itself produced.

### P0-6 Objective offset and sense are per model; integrality is per column

```c
int bplex_chg_col_type(bplex_t *bp, int j, int is_integer);   /* 0 continuous, 1 integer */
```

so the same object can be solved as an LP (`relax=1`) and, with the lambda
columns marked integer, as a MIP (`relax=0`) without a reload (used by P1-1).

## P1: what the primal heuristics and the tree need

### P1-1 MIP solve on the current model without disturbing the LP state

cvrpbcp solves the restricted master as a MIP from inside column generation.
Requirements:

* `bplex_solve` with `relax=0` on an object whose LP was warm must not destroy
  the LP basis, factors and weights; after the MIP, `relax=1` + `bplex_solve`
  continues the LP from where it was. (If the MIP works on a private copy, as
  `bp_mip_solve(p, par, res, x)` suggests, this is a matter of not touching
  `bp->warm`.)
* Parameters already present and sufficient: `time_limit`, `node_limit`,
  `cutoff`, `mip_rel_gap`, `mip_abs_gap`, `heuristics`, `heur_effort`,
  `threads`.
* Please confirm that `cutoff` is applied by the MIP as "prune nodes whose
  bound is >= cutoff" and that a run that finds nothing below the cutoff
  returns `BPLEX_CUTOFF` or `BPLEX_INFEASIBLE` (either is fine, documented).

### P1-2 MIP start

```c
/* Full-length start vector; infeasible or partial starts are ignored silently. */
int bplex_set_mip_start(bplex_t *bp, const double *x);
```

Internal status: `xstart` exists in `mip_int.h` (presolved space); the public
call needs the map into presolved space. Used to pass the incumbent routes so
the restricted master MIP starts with a cutoff and a feasible point.

### P1-3 Incumbent callback or solution pool

```c
typedef int (*bplex_incumbent_cb)(bplex_t *bp, double obj, const double *x, void *user);
/* cb returns 0 to continue, nonzero to stop the MIP (status BPLEX_INTERRUPTED). */
int bplex_set_incumbent_callback(bplex_t *bp, bplex_incumbent_cb cb, void *user);
```

Lets the heuristic stop as soon as any improving solution exists (we rarely
need the MIP optimum, we need a good solution fast) and lets a time-capped
MIP deliver every improving solution it found. A solution pool getter is an
acceptable substitute:

```c
int bplex_get_pool_size(const bplex_t *bp);
int bplex_get_pool_solution(const bplex_t *bp, int k, double *obj, double *x);
```

### P1-4 Objective cutoff in the LP dual simplex with a usable bound

`cutoff` exists and `bp_spx_set_cutoff` is implemented. Needed guarantees:

* when the dual objective exceeds the cutoff, the solve stops with
  `BPLEX_CUTOFF` and `bplex_bestbound` returns a valid dual bound (the dual
  objective of a dual-feasible iterate);
* `bplex_get_solution` returns the row duals of that iterate (we use them for
  the Lagrangian bound and for pricing one more time).

### P1-5 Basis snapshots for strong branching and diving

```c
typedef struct bplex_snap_s bplex_snap_t;
bplex_snap_t *bplex_snap_save(const bplex_t *bp);     /* basis + weights + factors */
int           bplex_snap_restore(bplex_t *bp, const bplex_snap_t *s);
void          bplex_snap_free(bplex_snap_t *s);
```

Internal status: `bp_spx_snap_save` / `bp_spx_snap_restore` exist for the
MIP's strong branching. Strong branching in cvrpbcp does 10-40 short dual
re-solves per node from the same point; `bplex_get_basis` + `bplex_set_basis`
work but refactorize every time.

### P1-6 Farkas certificate for an infeasible LP

```c
/* Dual ray y with y^T A x <= y^T b violated for all x in bounds. 0 if available. */
int bplex_get_farkas(const bplex_t *bp, double *y);
```

Internal status: `spx_verify_farkas` exists, the vector is not exposed.
cvrpbcp keeps artificial columns, so infeasibility is rare; this is for
Farkas pricing when we drop them.

### P1-7 Model dump for bug reports

```c
int bplex_write_mps(const bplex_t *bp, const char *path);   /* free MPS, full precision */
```

Any LP or MIP answer that disagrees with a published optimum will be reported
with this file and the expected objective. There is a reader but no writer.

## P2: performance and operational guarantees

* **P2-1 Column-rich masters.** n/m of 50-500 is normal. Please confirm that
  the dual simplex prices with partial / hyper-sparse pricing at that shape and
  that sifting is not triggered on a warm solve (it is skipped for `warm` today
  in `lp.c`, which is right).
* **P2-2 Degenerate LPs.** Set-partitioning masters stall; the perturbation
  (`perturb`) and the DSE weights are what prevent long tails. A per-solve
  `iter_limit` is used by strong branching; the status after the limit must
  still give the dual objective (P1-4).
* **P2-3 Thread safety across objects.** cvrpbcp will later evaluate nodes in
  parallel, one `bplex_t` per thread with `threads=1` each. Please confirm no
  shared mutable globals (the `bp_tune` table is read once, which is fine) and
  no shared scratch buffers.
* **P2-4 Deterministic work measure.** `bplex_get_info("iterations")` exists;
  a `"work"` counter (the CLI prints one) would let strong-branching budgets be
  deterministic.
* **P2-5 Reduced costs of nonbasic columns at zero cost.** Already returned by
  `bplex_get_solution`; please keep the sign convention `dj = c - A^T y`
  (as in `lp.c`) and document it in `bplex.h`.
* **P2-6 Infinity.** Document that `|bound| >= BPLEX_INF` is infinite and that
  `rowlo == rowup` is an equality row; cvrpbcp relies on both.

## Behaviour to confirm now (no code change expected)

1. `bplex_load` copies all arrays; the caller may free them after the call.
2. `bplex_set_basis` codes: `BPLEX_BASIC`, `BPLEX_AT_LOWER`, `BPLEX_AT_UPPER`,
   `BPLEX_AT_ZERO`; for rows "at lower" means the row activity is at `rowlo`.
   (`bplex.c` swaps lower/upper for rows internally; the public meaning should
   be stated in the header.)
3. A warm `bplex_solve` skips presolve (`lp.c`: `warm` bypasses it). Good for
   us; please keep it that way or make presolve return a basis.
4. After `BPLEX_OPTIMAL`, `bplex_get_basis` is always available, including when
   the model was presolved (postsolve restores a basis per the README).
5. `relax=1` on a model with integer columns solves the LP relaxation and
   returns LP duals.

## Workarounds in cvrpbcp until the P0 set exists (historical, no longer used)

| need | workaround | cost |
|---|---|---|
| add columns / rows | rebuild CSC arrays in cvrpbcp, `bplex_load`, `bplex_set_basis` with the old basis extended (new columns at lower, new rows basic) | O(nnz) per iteration, rescaling, DSE weights reset |
| delete columns / rows | rebuild without them, map the basis | same |
| change bounds (branching) | rebuild with the new bounds | same |
| strong branching | `bplex_get_basis` once, `bplex_set_basis` before each child solve | one refactorization per child |
| restricted master MIP | a second `bplex_t` loaded with pool columns and integrality, `relax=0`, `cutoff` and `time_limit` set | one extra load; no LP state conflict |
| MIP start / early stop | not available; rely on `cutoff` and `time_limit` | heuristic MIPs run to the time cap |
| Farkas | artificial columns in the master | none |
| model dumps | cvrpbcp writes its own MPS from its arrays | none |

## Round 2 (2026-10-26): what is still asked for a standalone split

**Status 2026-10-27: delivered in bplex API 0.2** (merge 2a6ea8e, history
item 208): R2-1 `bplex_get_pool_size` / `bplex_get_pool_solution` (pool of
the last MIP, `pool_size` parameter), R2-4 `BPLEX_INSTALL` with the
`bplex::bplex` package and `bplexConfigVersion.cmake` (SameMajorVersion),
R2-5 `bplex_version()` plus `BPLEX_VERSION_MINOR` bumps on header changes,
R2-6 `bplex_set_log_callback`, R2-7 `bplex_chg_col_bounds_many`, and an
error contract (`BPLEX_ERR_FATAL`: the object drops its model instead of
aborting the process). cvrpbcp now builds with `find_package(bplex 0.2
CONFIG)` against `cmake --install` output (`scripts/build_wsl.sh`), has no
reference to the bplex source tree any more, and records `bplex_version()`
in every JSONL result. The bullets below are kept as the original request.

Goal agreed with the bplex side: bplex and cvrpbcp are developed
independently and talk only through `include/bplex.h` and the CMake package.
From the dev report (API items 179/180, agents `apipool`, `cliquecoef`) and
our own check of main at 6571ebe:

Already in flight on the bplex side (`apipool`), needed by cvrpbcp:

* **R2-1 Solution pool getters** (P1-3 second half): `bplex_npool`,
  `bplex_get_pool_solution(bp, k, x, &obj)`. The restricted-master MIP yields
  several good partitions; cvrpbcp keeps all of them for the column pool and
  for diving. Until then the incumbent callback is used (it already works).
* **R2-2 Farkas after a cold solve** (P1-6): low priority for us; the master
  always carries artificial columns, so an infeasible LP never occurs.
* **R2-3 API hygiene**: no `abort()` reachable from the API, only `bplex_*`
  exported. Our `nm` of `libbplex.a` shows 285 non-`bplex_` globals, all
  prefixed `bp_`/`basiclu_`, so no clash with CVRPSEP (`CAPSEP_*`, `CMGR_*`)
  is expected either way.

New, small, and what makes the split clean:

* **R2-4 CMake package.** `CMakeLists.txt` has no `install()` or `export()`.
  cvrpbcp wants to build against an *installed* bplex
  (`cmake --install build/bplex --prefix <p>`; then
  `find_package(bplex CONFIG)` and `target_link_libraries(cvrpbcp bplex::bplex)`),
  not against the source tree, so neither project's build reads the other's
  sources. Needed: install of `include/bplex.h`, `libbplex.a`,
  `libbasiclu.a`, an exported target `bplex::bplex` carrying its link
  dependencies (`m`, Threads, MKL when on) and a `bplexConfigVersion.cmake`.
  Until it exists cvrpbcp links `build/bplex/libbplex.a` + `libbasiclu.a` by
  path and includes `../bplex/include`; that is the only cross-tree reference
  and it is confined to one CMake option (`BPLEX_ROOT`).
* **R2-5 Version at run time.** `const char *bplex_version(void)` (or
  `bplex_version_string`) returning what `BPLEX_VERSION_MAJOR.MINOR` plus the
  commit say, and a bump of `BPLEX_VERSION_MINOR` whenever the header
  changes. cvrpbcp prints it in every results line so a result can be tied
  to the bplex build that produced it.
* **R2-6 Log sink.** `bp_log` writes to `stdout` and `log_level` is per
  thread, set by `bplex_solve`. cvrpbcp runs with `log_level = 0` and prints
  its own log, which is enough for now. A `bplex_set_log_callback(bp, fn,
  user)` (or a `FILE *`) per object would let bplex's own log be kept in the
  per-instance log file when debugging a bad LP; low priority.
* **R2-7 Bulk bound change.** Edge branching sets the upper bound of every
  column using an edge to 0, often hundreds of columns at once.
  `bplex_chg_col_bounds` per column is O(1) each and fine; a
  `bplex_chg_col_bounds_many(bp, n, idx, lo, up)` would only save call
  overhead. Not needed.
* **R2-8 Behaviour change acknowledged.** A MIP that finds nothing better
  than `cutoff` returns `BPLEX_INFEASIBLE`; cvrpbcp treats that as "no
  improvement", not as an error, and never passes a cutoff below a known
  feasible value minus the gap tolerances.

Nothing in this list blocks milestones 0-3. R2-4 and R2-5 are the two that
matter for the standalone split; R2-1 is used from milestone 3.

### Performance datapoints from cvrpbcp (2026-10-27), not requests

* **Set-partitioning MIPs from route enumeration.** cvrpbcp enumerates all
  routes that can be in an improving solution and hands bplex a pure set
  partitioning (binary columns, n equality rows, one ranged fleet row, a few
  dozen RCC/SRC rows, objective cutoff). With 1000 or fewer columns bplex
  takes 0.2-5.5 s per solve (P-n50-k8, 49 rows), with 5-15k columns 16-60 s
  (XML150_1121_01, 150 rows), which is why cvrpbcp only sends sets of up to
  1000 routes and keeps branching otherwise. Dumps: `results/enum_mips/
  p50k8_*.mps` (written with `bplex_write_mps`, cutoff not in the file).
  Faster set-partitioning MIPs (clique/odd-set cuts, aggressive diving) would
  let cvrpbcp send larger sets and skip more of its tree.
* **Restricted-master MIP at the root** (5 s, 4-10k columns, same shape)
  rarely finds an incumbent better than the plunge; the first incumbent of
  XML150_1121_01 came after 20 s from a node LP.
* **Warm LP re-solves** behave as documented: 0 iterations after a no-op,
  few pivots after adding columns; the node LP is 20-50 ms with 2-4k columns
  and 300-1500 rows, which is now the dominant tree cost on n = 200.

## Bug report 2026-10-28: warm dual simplex stall on a 2.1k x 254 master LP

Reproducible with the LP trace of XML100_1113_07 (`results/traces/xml1113.trace`,
14 MB, 2856 solves; replay with `build/cvrpbcp-Release/lpreplay <trace> -v`):
solve 2854, an LP with ~2100 columns and 254 rows that the warm dual simplex
reaches after a column deletion/addition step like the 2853 before it (each
taking 15 ms on average), does not finish: recorded 99 s / 51204 iterations
until the time limit, replayed 99 s / 59364 iterations, status TIME_LIMIT.
The LP itself is harmless: `results/traces/xml1113_stall.mps` is that model
(dumped when the stall guard cut the solve after 2 s), a cold solve of it is
instant. The same run with the old build took 4 s in total, so the stall is
specific to this warm-start path (degenerate stalling or cycling of the dual
with the DSE weights it inherited).

Second observation: the first version of our guard retried the cut solve with
`method = BPLEX_METHOD_PRIMAL` and `time_limit` = remaining budget (150 s).
That warm primal solve ran for more than 12 minutes at 100 % CPU without
returning (killed), i.e. the primal warm path does not honour `time_limit`
(or loops). The guard now restarts the dual from the all-slack basis
(`bplex_set_basis`) instead.

Asked from bplex: (1) anti-stalling in the warm dual (perturbation /
bound shifting, or a cycling detector that re-factorizes and perturbs);
(2) `time_limit` and `iter_limit` honoured by the warm primal simplex too;
(3) optionally a `bplex_get_info("ncols")`/`("nrows")`. The trace format is
documented in `src/master/lptrace.h`; `tools/lpreplay.c` also has a CPLEX
backend for a like-for-like comparison once a Windows C compiler is at hand.

## Suggested acceptance tests on the bplex side

* Load a 3x3 LP, solve, add one column with negative reduced cost, solve
  again: iteration count <= 2, objective improves, basis valid.
* Add a violated row to an optimal LP, solve: dual simplex, 1-3 iterations.
* Delete a nonbasic column and a basic-slack row: objective unchanged, 0
  iterations on re-solve.
* Change a nonbasic column's upper bound to 0 (branching down): warm re-solve
  gives the same objective as a cold solve of the modified model.
* Solve a set-partitioning LP, mark columns integer, solve as MIP with
  `cutoff`, then `relax=1` and solve again: 0 iterations.
