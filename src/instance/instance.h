// Instance model for CVRP and CVRPTW, with the readers in instance.cpp.
//
// Conventions:
//   - vertex 0 is the depot, customers are 1..n.
//   - all costs are integers in "cost units": for CVRP one unit is one
//     distance unit (nint rounding, CVRPLIB convention); for CVRPTW one unit
//     is 0.1 distance units (DIMACS convention: distances and travel times
//     truncated to one decimal). Instance::scale is the number of units per
//     distance unit (1 or 10), so a published optimum 827.3 is 8273 units.
//   - time windows, service times and travel times are in the same units.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace cvrpbcp {

enum class Rounding { Nearest, TruncOneDecimal, Explicit };

struct Instance {
    std::string name;
    bool has_tw = false;        // CVRPTW
    int n = 0;                  // customers
    int capacity = 0;
    int vehicles = -1;          // -1: fleet size free
    int scale = 1;              // cost units per distance unit
    Rounding rounding = Rounding::Nearest;
    std::vector<double> x, y;   // size n+1, index 0 = depot
    std::vector<int> demand;    // size n+1, demand[0] = 0
    std::vector<int> service;   // cost units, size n+1 (0 unless has_tw)
    std::vector<int> tw_lo, tw_hi; // cost units, size n+1
    std::vector<int> dist;      // (n+1)*(n+1) cost units, row-major

    int nv() const { return n + 1; }
    int d(int i, int j) const { return dist[(size_t)i * nv() + j]; }
    // lower bound on the number of vehicles from total demand
    int min_vehicles() const;
    // first-fit-decreasing bin packing: an upper bound on the vehicles needed
    int ffd_vehicles() const;
    std::string summary() const;
};

// Read a CVRPLIB-style file (.vrp): TYPE CVRP or CVRPTW (PyVRP flavour with
// TIME_WINDOW_SECTION / SERVICE_TIME[_SECTION] / VEHICLES). Throws
// std::runtime_error on a malformed file. rounding_override: if not null,
// forces the distance rounding.
Instance read_instance(const std::string &path, const Rounding *rounding_override = nullptr);

// Published optimum parsed from the COMMENT line ("Optimal value: 784") when
// present, else -1.
long long comment_optimum(const std::string &path);

} // namespace cvrpbcp
