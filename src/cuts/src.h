// Limited-memory subset-row cuts on customer triples (Jepsen et al. 2008,
// Pecin et al. 2017): for C = {i,j,k}, sum_r floor(a_rC / 2) lambda_r <= 1.
// With node memory M (C subset of M) the coefficient of a route is the number
// of times its visit counter to C reaches 2 while staying inside M; leaving M
// resets the counter. That weakens the cut (coefficient <= floor(a_rC/2)) but
// keeps the label dominance damage local to M.
#pragma once
#include "instance/instance.h"
#include "master/master.h"

#include <vector>

namespace cvrpbcp {

struct SrcCut {
    std::vector<char> in_c; // size n+1
    std::vector<char> in_m; // size n+1, memory (includes C)
    int c[3];
    double violation;
};

// Separation on the current LP solution of M (routes with lambda > 0).
// `existing` are triples already in the master (encoded c0*nv*nv+c1*nv+c2,
// sorted), skipped. Returns the most violated triples first.
std::vector<SrcCut> separate_src(const Instance &I, const Master &M, const std::vector<long long> &existing,
                                 int max_cuts, double min_violation);

inline long long src_key(int nv, int a, int b, int c)
{
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return ((long long)a * nv + b) * nv + c;
}

} // namespace cvrpbcp
