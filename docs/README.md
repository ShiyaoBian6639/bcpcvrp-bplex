# exactsolv docs

* `PLAN.md` — cvrpbcp, the branch-cut-and-price CVRP/CVRPTW solver built on
  bplex: goals, architecture, milestones with gates, verification policy.
* `BPLEX_API_REQUESTS.md` — the public-API additions and guarantees cvrpbcp
  asks of bplex (incremental columns/rows, bound changes, warm-start contract,
  MIP start/callback, snapshots, Farkas, MPS writer), with workarounds until
  each exists. bplex itself is never modified from this project.
* `STATUS.md` — running log of milestones, benchmark results and design
  decisions (newest first).
