// Brute-force exact solvers for tiny instances (n <= ~9), used by the self
// test and by regression tests as an oracle independent of the BCP code.
#pragma once
#include "instance/instance.h"
#include "verify/verify.h"

#include <vector>

namespace cvrpbcp {

// All customer subsets that fit the capacity, with the cheapest feasible
// route through each (by enumeration of orders; time windows respected when
// the instance has them). Subsets with no feasible order are omitted.
struct EnumRoute {
    unsigned mask;             // bit c-1 set for customer c
    long long cost;
    std::vector<int> order;
};
std::vector<EnumRoute> enumerate_routes(const Instance &I);

// Exact CVRP(TW) optimum by dynamic programming over customer subsets, using
// enumerate_routes. max_vehicles <= 0: fleet free. cost -1 if infeasible.
Solution brute_force_optimum(const Instance &I, int max_vehicles = 0);

// Exact optimum by set partitioning over enumerate_routes solved as a MIP by
// bplex. cost -1 if bplex does not report an optimal solution.
Solution bplex_setpart_optimum(const Instance &I, int max_vehicles = 0, double time_limit = 60);

// Random Euclidean instance for tests (integer coordinates in [0,100]^2,
// integer demands in [1, max_demand]; with_tw adds random windows that keep
// the instance feasible).
Instance random_instance(int n, int capacity, int max_demand, unsigned seed, bool with_tw = false);

} // namespace cvrpbcp
