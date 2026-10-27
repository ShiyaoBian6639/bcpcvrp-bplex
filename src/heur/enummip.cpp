#include "heur/enummip.h"

#include "bplex.h"

#include <chrono>
#include <cmath>
#include <cstdio>

namespace cvrpbcp {

EnumMipResult solve_enum_mip(const Instance &I, const std::vector<Route> &routes, const std::vector<RobustRow> &rows,
                             int klo, int khi, double cutoff, double time_limit, int threads, int log_level,
                             const std::string &dump_path)
{
    using namespace std::chrono;
    auto t0 = steady_clock::now();
    EnumMipResult R;
    const int n = I.n;
    const int ncols = (int)routes.size();
    std::vector<const RobustRow *> act;
    for (const RobustRow &r : rows)
        if (r.active) act.push_back(&r);
    const int nrows = n + 1 + (int)act.size();
    std::vector<double> obj(ncols), lo(ncols, 0.0), up(ncols, 1.0), rlo(nrows), rup(nrows), vals;
    std::vector<char> integ(ncols, 1);
    std::vector<int> cs(ncols + 1), ri;
    std::vector<int> cnt(n + 1, 0);
    for (int j = 0; j < ncols; j++) {
        cs[j] = (int)ri.size();
        obj[j] = (double)routes[j].cost;
        for (int v : routes[j].seq) cnt[v]++;
        for (int i = 1; i <= n; i++)
            if (cnt[i]) { ri.push_back(i - 1); vals.push_back(cnt[i]); cnt[i] = 0; }
        ri.push_back(n);
        vals.push_back(1.0);
        for (size_t a = 0; a < act.size(); a++) {
            int c = act[a]->coef(routes[j]);
            if (c) { ri.push_back(n + 1 + (int)a); vals.push_back(c); }
        }
    }
    cs[ncols] = (int)ri.size();
    for (int i = 0; i < n; i++) rlo[i] = rup[i] = 1.0;
    rlo[n] = klo;
    rup[n] = khi;
    for (size_t a = 0; a < act.size(); a++) { rlo[n + 1 + a] = act[a]->lo; rup[n + 1 + a] = act[a]->up; }

    bplex_t *bp = bplex_create();
    if (!bp) { R.status = BPLEX_ERROR; return R; }
    bplex_set_param(bp, "log_level", log_level >= 3 ? 1 : 0);
    bplex_set_param(bp, "threads", threads);
    bplex_set_param(bp, "time_limit", time_limit);
    bplex_set_param(bp, "cutoff", cutoff);
    // Integer costs and cutoff = ub - 0.5: bplex prunes nodes with bound >=
    // cutoff - gap, so the gap must stay below 0.5 or a node holding a
    // solution of cost ub - 1 can be pruned (that happened with 0.999 and
    // produced a wrong "optimal" 537 on E-n30-k3).
    bplex_set_param(bp, "mip_abs_gap", 0.45);
    bplex_set_param(bp, "mip_rel_gap", 0.0);
    if (bplex_load(bp, ncols, nrows, obj.data(), lo.data(), up.data(), integ.data(), rlo.data(), rup.data(), cs.data(),
                   ri.data(), vals.data(), 0.0, 1) == 0) {
        if (!dump_path.empty()) bplex_write_mps(bp, dump_path.c_str());
        R.status = bplex_solve(bp);
        if (R.status == BPLEX_OPTIMAL || R.status == BPLEX_INFEASIBLE) R.solved = true;
        std::vector<double> x(ncols);
        if ((R.status == BPLEX_OPTIMAL || R.status == BPLEX_TIME_LIMIT || R.status == BPLEX_NODE_LIMIT) &&
            bplex_get_solution(bp, x.data(), nullptr, nullptr, nullptr) == 0) {
            for (int j = 0; j < ncols; j++)
                if (x[j] > 0.5) R.sol.routes.push_back(routes[j].seq);
            R.sol.cost = (long long)std::llround(bplex_objval(bp));
            R.has_solution = true;
        }
        if (R.status == BPLEX_OPTIMAL) R.bestbound = R.has_solution ? (double)R.sol.cost : cutoff;
        else if (R.status == BPLEX_INFEASIBLE) R.bestbound = cutoff;
        else {
            double b = bplex_bestbound(bp);
            R.bestbound = (std::fabs(b) < 1e299) ? std::min(b, cutoff) : -1e300;
        }
    } else {
        R.status = BPLEX_ERROR;
    }
    bplex_free(bp);
    R.time = duration<double>(steady_clock::now() - t0).count();
    return R;
}

} // namespace cvrpbcp
