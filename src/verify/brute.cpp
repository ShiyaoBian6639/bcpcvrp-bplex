#include "verify/brute.h"

#include "bplex.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace cvrpbcp {

namespace {

// cost of visiting `order` from and back to the depot; -1 if a time window
// is violated
long long route_cost(const Instance &I, const std::vector<int> &order)
{
    long long cost = 0, t = 0;
    int prev = 0;
    for (int c : order) {
        cost += I.d(prev, c);
        if (I.has_tw) {
            t = std::max<long long>(t + I.d(prev, c), I.tw_lo[c]);
            if (t > I.tw_hi[c]) return -1;
            t += I.service[c];
        }
        prev = c;
    }
    cost += I.d(prev, 0);
    if (I.has_tw && t + I.d(prev, 0) > I.tw_hi[0]) return -1;
    return cost;
}

} // namespace

std::vector<EnumRoute> enumerate_routes(const Instance &I)
{
    if (I.n > 12) throw std::runtime_error("enumerate_routes: n too large");
    std::vector<EnumRoute> out;
    const unsigned full = (1u << I.n) - 1;
    for (unsigned mask = 1; mask <= full; mask++) {
        long long load = 0;
        std::vector<int> cust;
        for (int c = 1; c <= I.n; c++)
            if (mask >> (c - 1) & 1) { load += I.demand[c]; cust.push_back(c); }
        if (load > I.capacity) continue;
        EnumRoute best{mask, -1, {}};
        std::sort(cust.begin(), cust.end());
        do {
            long long c = route_cost(I, cust);
            if (c >= 0 && (best.cost < 0 || c < best.cost)) { best.cost = c; best.order = cust; }
        } while (std::next_permutation(cust.begin(), cust.end()));
        if (best.cost >= 0) out.push_back(std::move(best));
    }
    return out;
}

Solution brute_force_optimum(const Instance &I, int max_vehicles)
{
    auto routes = enumerate_routes(I);
    const unsigned full = (1u << I.n) - 1;
    const int K = max_vehicles > 0 ? max_vehicles : I.n;
    const long long INF = (long long)4e18;
    // f[k][S]: cheapest cover of S with exactly k routes
    std::vector<std::vector<long long>> f(K + 1, std::vector<long long>(full + 1, INF));
    std::vector<std::vector<int>> arg(K + 1, std::vector<int>(full + 1, -1));
    f[0][0] = 0;
    // index routes by their lowest customer so each S is split canonically
    std::vector<std::vector<int>> by_low(I.n + 1);
    for (size_t r = 0; r < routes.size(); r++) {
        int low = __builtin_ctz(routes[r].mask) + 1;
        by_low[low].push_back((int)r);
    }
    for (int k = 1; k <= K; k++)
        for (unsigned S = 1; S <= full; S++) {
            int low = __builtin_ctz(S) + 1;
            long long best = INF;
            int bi = -1;
            for (int r : by_low[low]) {
                unsigned m = routes[r].mask;
                if ((m & S) != m) continue;
                long long v = f[k - 1][S ^ m];
                if (v < INF && v + routes[r].cost < best) { best = v + routes[r].cost; bi = r; }
            }
            f[k][S] = best;
            arg[k][S] = bi;
        }
    Solution sol;
    int bestk = -1;
    long long best = INF;
    for (int k = 1; k <= K; k++)
        if (f[k][full] < best) { best = f[k][full]; bestk = k; }
    if (bestk < 0) return sol;
    unsigned S = full;
    for (int k = bestk; k >= 1; k--) {
        int r = arg[k][S];
        sol.routes.push_back(routes[r].order);
        S ^= routes[r].mask;
    }
    sol.cost = best;
    return sol;
}

Solution bplex_setpart_optimum(const Instance &I, int max_vehicles, double time_limit)
{
    auto routes = enumerate_routes(I);
    const int ncols = (int)routes.size();
    const int nrows = I.n + (max_vehicles > 0 ? 1 : 0);
    std::vector<double> obj(ncols), lo(ncols, 0.0), up(ncols, 1.0), rlo(nrows), rup(nrows);
    std::vector<char> integ(ncols, 1);
    std::vector<int> colstart(ncols + 1), rowidx;
    std::vector<double> vals;
    for (int j = 0; j < ncols; j++) {
        colstart[j] = (int)rowidx.size();
        obj[j] = (double)routes[j].cost;
        for (int c = 1; c <= I.n; c++)
            if (routes[j].mask >> (c - 1) & 1) { rowidx.push_back(c - 1); vals.push_back(1.0); }
        if (max_vehicles > 0) { rowidx.push_back(I.n); vals.push_back(1.0); }
    }
    colstart[ncols] = (int)rowidx.size();
    for (int i = 0; i < I.n; i++) rlo[i] = rup[i] = 1.0;
    if (max_vehicles > 0) { rlo[I.n] = -BPLEX_INF; rup[I.n] = max_vehicles; }

    Solution sol;
    bplex_t *bp = bplex_create();
    if (!bp) return sol;
    bplex_set_param(bp, "log_level", 0);
    bplex_set_param(bp, "time_limit", time_limit);
    bplex_set_param(bp, "threads", 1);
    if (bplex_load(bp, ncols, nrows, obj.data(), lo.data(), up.data(), integ.data(), rlo.data(), rup.data(),
                   colstart.data(), rowidx.data(), vals.data(), 0.0, 1) == 0) {
        int st = bplex_solve(bp);
        if (st == BPLEX_OPTIMAL) {
            std::vector<double> x(ncols);
            if (bplex_get_solution(bp, x.data(), nullptr, nullptr, nullptr) == 0) {
                for (int j = 0; j < ncols; j++)
                    if (x[j] > 0.5) sol.routes.push_back(routes[j].order);
                sol.cost = std::llround(bplex_objval(bp));
            }
        }
    }
    bplex_free(bp);
    return sol;
}

Instance random_instance(int n, int capacity, int max_demand, unsigned seed, bool with_tw)
{
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> coord(0, 100), dem(1, max_demand);
    Instance I;
    I.name = "rand-n" + std::to_string(n) + "-s" + std::to_string(seed);
    I.n = n;
    I.capacity = capacity;
    I.has_tw = with_tw;
    I.scale = with_tw ? 10 : 1;
    I.rounding = with_tw ? Rounding::TruncOneDecimal : Rounding::Nearest;
    const int nv = n + 1;
    I.x.resize(nv); I.y.resize(nv); I.demand.assign(nv, 0);
    I.service.assign(nv, 0); I.tw_lo.assign(nv, 0); I.tw_hi.assign(nv, 0);
    for (int v = 0; v < nv; v++) { I.x[v] = coord(rng); I.y[v] = coord(rng); }
    for (int c = 1; c <= n; c++) I.demand[c] = dem(rng);
    I.dist.assign((size_t)nv * nv, 0);
    for (int i = 0; i < nv; i++)
        for (int j = 0; j < nv; j++) {
            double dx = I.x[i] - I.x[j], dy = I.y[i] - I.y[j];
            double e = std::sqrt(dx * dx + dy * dy);
            I.dist[(size_t)i * nv + j] = with_tw ? (int)std::floor(e * 10 + 1e-9) : (int)std::llround(e);
        }
    if (with_tw) {
        // window around the direct arrival time so single-customer routes are
        // always feasible; horizon generous
        std::uniform_int_distribution<int> width(200, 800), svc(0, 100);
        int horizon = 0;
        for (int c = 1; c <= n; c++) {
            int a = I.d(0, c);
            I.tw_lo[c] = std::max(0, a - width(rng) / 2);
            I.tw_hi[c] = a + width(rng);
            I.service[c] = svc(rng);
            horizon = std::max(horizon, I.tw_hi[c] + I.service[c] + I.d(c, 0));
        }
        I.tw_lo[0] = 0;
        I.tw_hi[0] = horizon;
    }
    return I;
}

} // namespace cvrpbcp
