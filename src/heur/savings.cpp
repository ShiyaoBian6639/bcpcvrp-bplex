#include "heur/savings.h"

#include <algorithm>
#include <numeric>

namespace cvrpbcp {

Solution clarke_wright(const Instance &I)
{
    const int n = I.n;
    struct Sav { long long s; int i, j; };
    std::vector<Sav> sav;
    sav.reserve((size_t)n * (n - 1) / 2);
    for (int i = 1; i <= n; i++)
        for (int j = i + 1; j <= n; j++) sav.push_back({(long long)I.d(0, i) + I.d(j, 0) - I.d(i, j), i, j});
    std::sort(sav.begin(), sav.end(), [](const Sav &a, const Sav &b) {
        if (a.s != b.s) return a.s > b.s;
        return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    // route id per customer; each route kept as a deque-like vector
    std::vector<int> rid(n + 1);
    std::vector<std::vector<int>> routes(n + 1);
    std::vector<int> load(n + 1);
    for (int i = 1; i <= n; i++) { rid[i] = i; routes[i] = {i}; load[i] = I.demand[i]; }
    for (const Sav &s : sav) {
        int ri = rid[s.i], rj = rid[s.j];
        if (ri == rj) continue;
        if (load[ri] + load[rj] > I.capacity) continue;
        auto &A = routes[ri], &B = routes[rj];
        bool i_end = A.back() == s.i, i_start = A.front() == s.i;
        bool j_end = B.back() == s.j, j_start = B.front() == s.j;
        if (!(i_end || i_start) || !(j_end || j_start)) continue; // interior customers
        // orient so that A ends with i and B starts with j
        if (!i_end) std::reverse(A.begin(), A.end());
        if (!j_start) std::reverse(B.begin(), B.end());
        for (int c : B) { A.push_back(c); rid[c] = ri; }
        load[ri] += load[rj];
        B.clear();
        load[rj] = 0;
    }
    Solution S;
    for (int r = 1; r <= n; r++)
        if (!routes[r].empty()) S.routes.push_back(routes[r]);
    CheckResult c = check_solution(I, S);
    S.cost = c.feasible ? c.cost : -1;
    return S;
}

} // namespace cvrpbcp
