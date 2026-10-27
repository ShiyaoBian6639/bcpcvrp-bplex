// Independent solution checker. It knows nothing about the solver: it takes
// routes as customer sequences and recomputes everything from the instance.
#pragma once
#include "instance/instance.h"

#include <string>
#include <vector>

namespace cvrpbcp {

struct Solution {
    std::vector<std::vector<int>> routes; // customers 1..n, depot implicit
    long long cost = -1;                  // cost units as reported by the writer, -1 if none
};

struct CheckResult {
    bool feasible = false;
    long long cost = 0;        // recomputed, cost units
    std::string error;         // first violation found
    int vehicles = 0;
};

// Every customer exactly once, loads <= capacity, time windows (if any),
// vehicles <= instance.vehicles when the instance fixes it (enforce_fleet).
CheckResult check_solution(const Instance &I, const Solution &S, bool enforce_fleet = false);

// CVRPLIB .sol format: lines "Route #k: 3 7 2" (customers numbered 1..n,
// depot omitted) and "Cost 784". cost_scale multiplies a decimal cost into
// cost units (10 for the one-decimal CVRPTW convention).
Solution read_sol(const std::string &path, int cost_scale = 1);
void write_sol(const std::string &path, const Solution &S, int cost_scale = 1);
std::string format_sol(const Solution &S, int cost_scale = 1);

} // namespace cvrpbcp
