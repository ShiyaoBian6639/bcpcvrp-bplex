// Clarke-Wright parallel savings: a quick feasible CVRP solution used to
// seed the master (sensible first duals) and as a first incumbent.
#pragma once
#include "instance/instance.h"
#include "verify/verify.h"

namespace cvrpbcp {

// Capacity only (time windows ignored: for CVRPTW the routes are still
// valid columns only if the checker accepts them; the caller verifies).
Solution clarke_wright(const Instance &I);

} // namespace cvrpbcp
