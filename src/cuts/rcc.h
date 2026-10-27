// Rounded capacity cuts through CVRPSEP (Lysgaard). CVRPSEP numbers the
// customers 1..n and the depot n+1 and returns x(S:S) <= |S| - k(S); we hand
// back S and k(S) for the crossing form x(delta(S)) >= 2 k(S).
#pragma once
#include "instance/instance.h"

#include <vector>

namespace cvrpbcp {

struct RccCut {
    std::vector<char> in_s; // size n+1, in_s[0] = 0
    int k;
    double violation;
};

class RccSeparator {
public:
    explicit RccSeparator(const Instance &I);
    ~RccSeparator();
    RccSeparator(const RccSeparator &) = delete;
    RccSeparator &operator=(const RccSeparator &) = delete;

    // x: symmetric (n+1)^2 edge flows. existing: the RCC sets currently in
    // the master (in_s vectors), which CVRPSEP must not return again.
    // Returns violated cuts, most violated first.
    std::vector<RccCut> separate(const std::vector<double> &x, const std::vector<const std::vector<char> *> &existing,
                                 int max_cuts, double min_violation, bool &integer_feasible);

private:
    const Instance &I_;
    void *new_cuts_ = nullptr, *old_cuts_ = nullptr; // CnstrMgrPointer (anonymous struct in cnstrmgr.h)
    std::vector<double> demand_;
    std::vector<int> list_;
};

} // namespace cvrpbcp
