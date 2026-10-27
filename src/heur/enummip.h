// Set-partitioning MIP over an enumerated route set, solved by bplex in a
// fresh object. Used after route enumeration: when every route of an
// improving solution is in the set, the MIP finishes the node.
#pragma once
#include "core/types.h"
#include "instance/instance.h"
#include "master/master.h"
#include "verify/verify.h"

#include <vector>

namespace cvrpbcp {

struct EnumMipResult {
    int status = 0;          // bplex status
    bool solved = false;     // OPTIMAL, or INFEASIBLE under the cutoff (no better solution exists)
    bool has_solution = false;
    Solution sol;            // when has_solution
    double bestbound = 0;    // valid bound on improving solutions (min with cutoff)
    double time = 0;
};

// rows: the master's rows; only active ones are added (their coefficients
// are computed per route). cutoff: exclusive upper bound on the cost.
EnumMipResult solve_enum_mip(const Instance &I, const std::vector<Route> &routes, const std::vector<RobustRow> &rows,
                             int klo, int khi, double cutoff, double time_limit, int threads, int log_level,
                             const std::string &dump_path = "");

} // namespace cvrpbcp
