// Column (and cut) generation at one node.
#pragma once
#include "core/types.h"
#include "cuts/rcc.h"
#include "instance/instance.h"
#include "master/master.h"
#include "pricing/labeling.h"
#include "verify/verify.h"

#include <string>
#include <vector>

namespace cvrpbcp {

struct NodeResult {
    bool ok = false;         // finished without error (converged, early-stopped, or infeasible)
    bool infeasible = false; // LP needs artificial columns: node has no feasible solution
    double lp = 0;           // last master objective
    double lb = 0;           // best Lagrangian bound (valid lower bound for the node)
    long long lb_int = 0;    // ceil(lb) in cost units
    int cg_iters = 0;
    int cut_rounds = 0;
    int ncols = 0, ncuts = 0;
    double time = 0, time_pricing = 0, time_master = 0, time_cuts = 0;
    bool integral = false;   // an integer feasible solution was read off the LP
    bool fleet_binding = false; // sum lambda == khi
    bool early_stopped = false; // CG stopped on ceil(lb) == ceil(lp): duals are not LP-optimal
    // Lagrangian bound and minimum reduced cost of the last exact pricing;
    // the pricer still holds those duals (used for route enumeration)
    double lb_last = -1e300, minrc_last = 0;
    Solution sol;            // when integral
    std::vector<double> x;   // edge flows of the final LP (for branching)
    double fleet = 0;        // sum lambda
    std::string error;
};

// Runs CG with RCC rounds on master M. cutoff: stop (node pruned) once the
// bound reaches it (cost units; pass a huge value for none). enum_set: when
// given, pricing scans this route set (all routes of improving solutions in
// the subtree) instead of labeling.
NodeResult solve_node(const Instance &I, const Params &P, Master &M, Pricer &pricer, RccSeparator *rcc,
                      double time_budget, double cutoff = 1e300, const std::vector<Route> *enum_set = nullptr);

// Routes of `set` with reduced cost < -eps under the pricer's duals (and no
// forbidden edge), best first, at most max_routes; min_rc over the set.
std::vector<Route> scan_routes(const std::vector<Route> &set, const Pricer &pr, double eps, int max_routes,
                               double &min_rc);

// Builds the master with singleton routes, runs solve_node at the root.
NodeResult solve_root(const Instance &I, const Params &P, int kmax);

// Extract routes from integral edge flows (each customer degree 2). Returns
// false if x is not integral or not a set of depot cycles.
bool routes_from_flows(const Instance &I, const std::vector<double> &x, Solution &S);

} // namespace cvrpbcp
