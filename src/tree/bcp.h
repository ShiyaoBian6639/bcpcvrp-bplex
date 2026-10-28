// Branch-cut-and-price tree: best-bound search with plunging, branching on
// fleet size and on edge flows, incumbents from integral LPs and from the
// restricted master MIP on bplex.
#pragma once
#include "core/types.h"
#include "instance/instance.h"
#include "verify/verify.h"

#include <string>

namespace cvrpbcp {

struct BcpParams : Params {
    int kmax = 0;                 // 0: ffd + 1
    bool free_fleet = false;      // ignore the instance's vehicle count (CLI applies it to the instance)
    // when the instance states a vehicle count: 0 = exactly that many (CVRPLIB
    // A/B/E/F/P convention: P-n55-k8's 588 needs exactly 8), 1 = at most
    int fleet_mode = 0;
    double ub = 1e300;            // external upper bound (cost units), exclusive cutoff is ub
    bool rmip = true;             // restricted master MIP heuristic
    double rmip_time = 5;         // seconds for the root call
    double rmip_time_node = 2;    // seconds for periodic calls in the tree
    int rmip_every = 50;          // nodes between periodic calls
    int sb_cands = 0;             // strong-branching probes per node (0: pseudocost branching only;
                                  // probes with a 100-iteration cap lost to pseudocosts on P-n55-k7)
    int sb_iters = 100;           // dual simplex iterations per probe
    int sb_lookahead = 4;         // stop probing after this many candidates without a better score
    int sb_reliable = 2;          // pseudocost observations per direction before an edge skips probing
    int node_limit = 1000000;
    int plunge_depth = 50;        // keep diving on a child while it is not pruned, up to this depth
    // column-fixing dive heuristic: at the root, and every dive_every nodes while no incumbent exists
    bool dive = true;
    int dive_every = 20;
    double dive_time = 30;        // seconds per dive (120 let a failing root dive eat 130 s of a 300 s run)
    // route enumeration + set-partitioning MIP once an incumbent exists
    bool enumerate = true;
    size_t enum_max_labels = 1000000;
    size_t enum_max_routes = 300000;  // enumeration gives up above this
    size_t enum_mip_max_routes = 1000; // MIP only for sets up to this size; larger sets replace labeling
    double enum_mip_time = 60;    // seconds per MIP
    double enum_time_cap = 10;    // seconds per enumeration attempt (a slow attempt counts as failed)
    std::string enum_dump;        // if set, write each enumeration MIP as <prefix>_<k>.mps (for the bplex side)
};

struct BcpResult {
    std::string status;           // "optimal", "time_limit", "node_limit", "infeasible", "error"
    long long lb = 0;             // global lower bound (cost units)
    long long ub = -1;            // incumbent cost (cost units), -1 if none
    Solution sol;
    double root_lp = 0, root_lb = 0;
    int nodes = 0, depth_max = 0, ncols = 0, nrows = 0;
    double time = 0, time_root = 0;
    double t_pricing = 0, t_master = 0, t_cuts = 0, t_rmip = 0, t_branch = 0, t_switch = 0, t_enum = 0;
    double t_enum_label = 0, t_enum_mip = 0, t_dive = 0;
    int enum_tries = 0, enum_closed = 0, enum_mips = 0;
    double t_sep = 0, t_src = 0;      // parts of t_cuts
    double t_pr[6] = {0, 0, 0, 0, 0, 0}; // pricing split, NodeResult::PR_*
    int n_pr[6] = {0, 0, 0, 0, 0, 0};
    int h_pr[6] = {0, 0, 0, 0, 0, 0}; // calls with columns
    int cg_iters = 0;
    int n_sep = 0, n_sep_hit = 0;
    bool fleet_cap_binding = false; // the incumbent uses kmax vehicles (re-run with larger kmax)
    int kmin = 0, kmax = 0;
    std::string error;
};

BcpResult solve_bcp(const Instance &I, const BcpParams &P);

} // namespace cvrpbcp
