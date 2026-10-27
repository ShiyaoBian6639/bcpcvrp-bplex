#include "tree/colgen.h"

#include "bplex.h"
#include "cuts/src.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <unordered_set>

namespace cvrpbcp {

namespace {
double now()
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}
long long ceil_int(double v) { return (long long)std::ceil(v - 1e-6); }
} // namespace

bool routes_from_flows(const Instance &I, const std::vector<double> &x, Solution &S)
{
    const int nv = I.nv();
    std::vector<std::vector<int>> adj(nv);
    for (int i = 0; i < nv; i++)
        for (int j = i + 1; j < nv; j++) {
            double v = x[(size_t)i * nv + j];
            int k = (int)std::lround(v);
            if (std::fabs(v - k) > 1e-6) return false;
            for (int t = 0; t < k; t++) { adj[i].push_back(j); adj[j].push_back(i); }
        }
    for (int c = 1; c < nv; c++)
        if (adj[c].size() != 2) return false;
    S.routes.clear();
    std::vector<char> used(nv, 0);
    for (int start : adj[0]) {
        if (used[start]) continue;
        std::vector<int> r;
        int prev = 0, cur = start;
        while (cur != 0) {
            if (used[cur]) return false;
            used[cur] = 1;
            r.push_back(cur);
            int next = (adj[cur][0] == prev) ? adj[cur][1] : adj[cur][0];
            if (adj[cur][0] == prev && adj[cur][1] == prev) next = prev; // 0-c-0 route
            prev = cur;
            cur = next;
        }
        S.routes.push_back(r);
    }
    for (int c = 1; c < nv; c++)
        if (!used[c]) return false;
    return true;
}

std::vector<Route> scan_routes(const std::vector<Route> &set, const Pricer &pr, double eps, int max_routes,
                               double &min_rc)
{
    struct C { double rc; int r; };
    std::vector<C> cand;
    min_rc = 0;
    for (size_t r = 0; r < set.size(); r++) {
        const Route &rt = set[r];
        int prev = 0;
        bool bad = false;
        for (int v : rt.seq) {
            if (pr.forbidden(prev, v)) { bad = true; break; }
            prev = v;
        }
        if (bad || pr.forbidden(prev, 0)) continue;
        double rc = pr.route_rc(rt);
        if (rc < min_rc) min_rc = rc;
        if (rc < -eps) cand.push_back({rc, (int)r});
    }
    std::sort(cand.begin(), cand.end(), [](const C &a, const C &b) { return a.rc < b.rc; });
    if ((int)cand.size() > max_routes) cand.resize(max_routes);
    std::vector<Route> out;
    out.reserve(cand.size());
    for (const C &c : cand) out.push_back(set[c.r]);
    return out;
}

NodeResult solve_node(const Instance &I, const Params &P, Master &M, Pricer &pricer, RccSeparator *rcc,
                      double time_budget, double cutoff, const std::vector<Route> *enum_set)
{
    NodeResult R;
    const double t0 = now();
    std::vector<double> pi_in, y_in; // stability centre (best Lagrangian duals)
    std::vector<double> pi_sep, y_sep;
    double alpha = P.smoothing_alpha;
    bool misprice = false;
    R.lb = -1e300;
    bool integer_feasible = false;
    bool pruned = false;
    int src_rounds = 0;
    double lp_before_src = -1, lp_before_src2 = -1;
    int fail_s0 = 0, fail_s1 = 0; // consecutive empty results of the heuristic stages (node-wide)

    while (true) {
        // ---- column generation to convergence for the current rows ----
        while (true) {
            if (now() - t0 > time_budget) { R.error = "time limit"; return R; }
            double tm = now();
            M.set_lp_time_limit(std::max(1.0, time_budget - (now() - t0)));
            if (P.log_level >= 4) std::printf("  lp solve (%d cols, %zu rows)...\n", M.nlpcols(), M.rows().size());
            int st = M.solve();
            R.time_master += now() - tm;
            R.cg_iters++;
            if (st != BPLEX_OPTIMAL) {
                R.error = std::string("master LP status ") + bplex_status_string(st);
                M.write_mps("cvrpbcp_master_fail.mps");
                return R;
            }
            R.lp = M.objval();
            const std::vector<double> &pi_out = M.pi();
            const std::vector<double> &y_out = M.row_duals();
            double a = (alpha > 0 && !pi_in.empty() && !misprice) ? alpha : 0.0;
            pi_sep.resize(pi_out.size());
            y_sep.resize(y_out.size());
            for (size_t i = 0; i < pi_out.size(); i++)
                pi_sep[i] = (a > 0 ? a * pi_in[i] : 0.0) + (1 - a) * pi_out[i];
            for (size_t c = 0; c < y_out.size(); c++)
                y_sep[c] = (a > 0 && c < y_in.size() ? a * y_in[c] : 0.0) + (1 - a) * y_out[c];
            misprice = false;

            double tp = now();
            if (R.cg_iters == 1 || (P.purge_every > 0 && R.cg_iters % P.purge_every == 0))
                M.purge(P.purge_age, P.purge_min_cols_per_cust * I.n);
            pricer.set_duals(pi_sep, M.rows(), y_sep);
            double min_rc = 0;
            std::vector<Route> cols;
            bool exact = false;
            int from_pool = 0;
            double tq = now();
            auto lap = [&](int what) { double t = now(); R.time_pr[what] += t - tq; R.calls_pr[what]++; tq = t; };
            lap(NodeResult::PR_SETUP);
            if (enum_set) {
                // the subtree's improving solutions only use these routes
                cols = scan_routes(*enum_set, pricer, P.rc_eps, P.max_cols_per_iter, min_rc);
                exact = true;
                lap(NodeResult::PR_SCAN);
            } else {
                from_pool = M.price_pool(pricer, P.rc_eps, P.max_cols_per_iter, R.cg_iters == 1);
                lap(NodeResult::PR_POOL);
                if (from_pool == 0) {
                    // adaptive stages: near convergence the sparse and the
                    // capped passes mostly come back empty; after
                    // stage_skip_after empty results in a row a stage is
                    // skipped for the rest of the node
                    if (fail_s0 < P.stage_skip_after) {
                        cols = pricer.price(0, P.max_cols_per_iter, min_rc);
                        lap(NodeResult::PR_S0);
                        fail_s0 = cols.empty() ? fail_s0 + 1 : 0;
                    }
                    if ((int)cols.size() < P.heur_min_cols && fail_s1 < P.stage_skip_after) {
                        std::vector<Route> more = pricer.price(1, P.max_cols_per_iter, min_rc);
                        fail_s1 = more.empty() ? fail_s1 + 1 : 0;
                        if (more.size() > cols.size()) cols.swap(more);
                        lap(NodeResult::PR_S1);
                    }
                    if (cols.empty()) {
                        cols = pricer.price(2, P.max_cols_per_iter, min_rc);
                        exact = true;
                        lap(NodeResult::PR_S2);
                    }
                }
            }
            R.time_pricing += now() - tp;
            if (exact) {
                double lb = M.lagrangian_bound(pi_sep, y_sep, min_rc);
                R.lb_last = lb;
                R.minrc_last = min_rc;
                if (lb > R.lb) {
                    R.lb = lb;
                    pi_in = pi_sep;
                    y_in = y_sep;
                }
                // automatic alpha: compare the in-out direction with a
                // subgradient of the Lagrangian at pi_sep (customer part)
                if (alpha > 0 && P.smoothing_auto && a > 0 && !cols.empty()) {
                    const Route &best = cols.front();
                    std::vector<double> g(pi_out.size(), 1.0);
                    g[0] = 0;
                    for (int v : best.seq) g[v] -= (double)M.khi();
                    double dot = 0, nd = 0, ng = 0;
                    for (size_t i = 1; i < pi_out.size(); i++) {
                        double d = pi_out[i] - pi_in[i];
                        dot += d * g[i];
                        nd += d * d;
                        ng += g[i] * g[i];
                    }
                    if (nd > 0 && ng > 0) {
                        if (dot > 0) alpha = std::min(alpha + 0.1, 0.99);   // in-out direction improves: smooth more
                        else alpha = std::max(alpha - 0.1, 0.0);           // else trust the LP duals more
                    }
                }
            }
            if (P.log_level >= 3 || (P.log_level >= 2 && (R.cg_iters % 10 == 1 || (cols.empty() && !from_pool))))
                std::printf("  it %4d  z %12.4f  lb %12.4f  minrc %10.4f  lp %5d/%6d  rows %4d  %s +%zu  lbl %zu  spx %d  %.1fs\n",
                            R.cg_iters, R.lp, std::max(R.lb, -1e9), min_rc, M.nlpcols(), M.nroutes(),
                            (int)M.rows().size(), from_pool ? "pl" : exact ? "ex" : "hr",
                            from_pool ? (size_t)from_pool : cols.size(), pricer.labels_last(), M.last_iterations(),
                            now() - t0);
            if (exact && ceil_int(R.lb) >= cutoff) { pruned = true; break; }
            if (from_pool) continue;
            if (cols.empty()) {
                if (a > 0) { misprice = true; continue; } // smoothed duals mispriced: retry with LP duals
                break;                                     // converged
            }
            int added = M.add_routes(cols);
            if (added == 0) {
                // every priced route is already in the LP: a misprice under
                // smoothed duals, or a tolerance mismatch with the LP (converged)
                if (a > 0) { misprice = true; continue; }
                break;
            }
            if (P.early_stop && exact && ceil_int(R.lb) >= ceil_int(R.lp)) { R.early_stopped = true; break; }
            R.early_stopped = false;
        }
        R.ncols = M.nroutes();
        R.ncuts = (int)M.rows().size();
        if (pruned) break;
        if (M.artificial_mass() > 1e-6) {
            // only artificials can satisfy the rows: the node is infeasible
            // (the LP value with bigM is huge, so this also prunes by bound)
            R.infeasible = true;
            break;
        }

        // ---- cuts: RCC first; when none is violated, limited-memory SRCs ----
        if (R.cut_rounds >= P.rcc_max_rounds) break;
        double tc = now();
        M.edge_flows(R.x);
        size_t added = 0;
        const char *what = "RCC";
        double maxviol = 0;
        if (rcc && P.rcc) {
            // CVRPSEP's "existing cuts" list costs time per entry; in the tree
            // the master holds many RCC rows, most with a zero dual, so only
            // the binding ones are handed over and duplicates of the rest
            // are filtered by the set hash
            std::vector<const std::vector<char> *> existing;
            std::unordered_set<uint64_t> have;
            auto set_hash = [](const std::vector<char> &s) {
                uint64_t h = 1469598103934665603ull;
                for (char c : s) h = (h ^ (uint64_t)(unsigned char)c) * 1099511628211ull;
                return h;
            };
            for (const RobustRow &r : M.rows())
                if (r.kind == RobustRow::RCC) {
                    if (P.rcc_existing_all || r.zero_streak == 0) existing.push_back(&r.in_s);
                    else have.insert(set_hash(r.in_s));
                }
            const double ts = now();
            auto cuts = rcc->separate(R.x, existing, P.rcc_max_per_round, P.rcc_min_violation, integer_feasible);
            R.time_sep += now() - ts;
            R.sep_calls++;
            if (!have.empty())
                cuts.erase(std::remove_if(cuts.begin(), cuts.end(),
                                          [&](const RccCut &c) { return have.count(set_hash(c.in_s)) > 0; }),
                           cuts.end());
            if (!cuts.empty()) R.sep_hits++;
            std::vector<RobustRow> rows;
            rows.reserve(cuts.size());
            for (const auto &c : cuts) {
                RobustRow r;
                r.kind = RobustRow::RCC;
                r.in_s = c.in_s;
                r.k = c.k;
                r.lo = 2.0 * c.k;
                r.up = BPLEX_INF;
                r.active = true;
                rows.push_back(std::move(r));
            }
            M.add_rows(rows); // one batch: one bplex edit for the round
            added = cuts.size();
            if (added) maxviol = cuts.front().violation;
        }
        if (added == 0 && P.src) {
            // tail-off: two SRC rounds in a row with little progress end the cutting
            if (src_rounds >= 2 && R.lp - lp_before_src < P.cut_tailoff * std::fabs(R.lp) &&
                lp_before_src2 >= 0 && lp_before_src - lp_before_src2 < P.cut_tailoff * std::fabs(R.lp)) {
                R.time_cuts += now() - tc;
                break;
            }
            int active = 0;
            std::vector<long long> existing;
            for (const RobustRow &r : M.rows())
                if (r.kind == RobustRow::SRC) {
                    if (r.active) active++;
                    int c[3], t = 0;
                    for (int v = 1; v <= I.n && t < 3; v++) if (r.in_s[v]) c[t++] = v;
                    existing.push_back(src_key(I.nv(), c[0], c[1], c[2]));
                }
            int room = std::min(P.src_max_per_round, P.src_max_active - active);
            if (room > 0) {
                const double ts = now();
                auto cuts = separate_src(I, M, existing, room, P.src_min_violation);
                R.time_src += now() - ts;
                std::vector<RobustRow> rows;
                rows.reserve(cuts.size());
                for (const auto &c : cuts) {
                    RobustRow r;
                    r.kind = RobustRow::SRC;
                    r.in_s = c.in_c;
                    r.in_m = c.in_m;
                    r.lo = -BPLEX_INF;
                    r.up = 1.0;
                    r.active = true;
                    rows.push_back(std::move(r));
                }
                M.add_rows(rows);
                added = cuts.size();
                what = "SRC";
                if (added) {
                    maxviol = cuts.front().violation;
                    src_rounds++;
                    lp_before_src2 = lp_before_src;
                    lp_before_src = R.lp;
                }
            }
        }
        R.time_cuts += now() - tc;
        if (added == 0) break;
        R.cut_rounds++;
        if (P.log_level >= 2)
            std::printf("  cut round %d: +%zu %s (max viol %.3f), z %.4f\n", R.cut_rounds, added, what, maxviol, R.lp);
    }

    R.lb_int = ceil_int(std::max(R.lb, -1e300));
    R.ok = true;
    R.time = now() - t0;
    if (pruned || R.infeasible) return R;

    // ---- integrality ----
    const auto &lam = M.lambda();
    double sum = 0;
    bool integral = true;
    Solution S;
    for (size_t r = 0; r < lam.size(); r++) {
        sum += lam[r];
        if (lam[r] > 1e-6 && lam[r] < 1 - 1e-6) integral = false;
        if (lam[r] >= 1 - 1e-6) S.routes.push_back(M.route((int)r).seq);
    }
    R.fleet = sum;
    R.fleet_binding = sum >= M.khi() - 1e-6;
    if (R.x.empty()) M.edge_flows(R.x);
    if (!integral) {
        // lambda fractional but the edge flows may still be integral
        Solution T;
        if (routes_from_flows(I, R.x, T)) { S = T; integral = true; }
    }
    if (integral) {
        CheckResult c = check_solution(I, S);
        if (c.feasible) {
            S.cost = c.cost;
            R.sol = S;
            R.integral = true;
        }
    }
    return R;
}

NodeResult solve_root(const Instance &I, const Params &P, int kmax)
{
    int kmin = I.min_vehicles();
    if (kmax <= 0) kmax = std::max(I.ffd_vehicles(), kmin) + 1;
    if (I.vehicles > 0) kmax = std::max(std::min(kmax, I.vehicles), kmin);
    Master M(I, kmin, kmax, P);
    std::vector<Route> init;
    for (int i = 1; i <= I.n; i++) {
        Route r;
        r.seq = {i};
        r.cost = I.d(0, i) + I.d(i, 0);
        init.push_back(r);
    }
    M.add_routes(init);
    Pricer pricer(I, P);
    RccSeparator rcc(I);
    if (P.log_level >= 1)
        std::printf("root: %s  kmin %d kmax %d  ng %d\n", I.summary().c_str(), kmin, kmax, P.ng_size);
    return solve_node(I, P, M, pricer, &rcc, P.time_limit);
}

} // namespace cvrpbcp
