#include "tree/bcp.h"

#include "bplex.h"
#include "cuts/rcc.h"
#include "heur/enummip.h"
#include "heur/savings.h"
#include "master/master.h"
#include "pricing/labeling.h"
#include "tree/colgen.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>

namespace cvrpbcp {

static void add_times(BcpResult &R, const NodeResult &res)
{
    R.t_pricing += res.time_pricing;
    R.t_master += res.time_master;
    R.t_cuts += res.time_cuts;
    R.t_sep += res.time_sep;
    R.t_src += res.time_src;
    R.cg_iters += res.cg_iters;
    R.n_sep += res.sep_calls;
    R.n_sep_hit += res.sep_hits;
    for (int k = 0; k < NodeResult::PR_N; k++) { R.t_pr[k] += res.time_pr[k]; R.n_pr[k] += res.calls_pr[k]; }
}

namespace {

double now()
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct Node {
    int id = 0, parent = -1, depth = 0;
    double lb = 0;                 // parent's bound at creation, then own
    std::vector<int> rows;         // active branch row ids
    int klo = 0, khi = 0;
    // how this node was created (for pseudocost updates)
    int bedge = -1;                // u*nv+v, -1 for fleet branching / root
    int bdir = 0;                  // 0 down, 1 up
    double bfrac = 0;              // distance moved: f (down) or 1-f (up)
    double parent_lp = 0;
    std::shared_ptr<std::vector<Route>> enum_set; // routes of all improving solutions in this subtree
    double enum_fail_gap = 1e300;  // gap at which enumeration last failed in this branch (retry once it halves)
};

struct PseudoCost {
    double sum[2] = {0, 0};
    int cnt[2] = {0, 0};
};

} // namespace

BcpResult solve_bcp(const Instance &I, const BcpParams &P)
{
    BcpResult R;
    const double t0 = now();
    int kmin = I.min_vehicles();
    int kmax = P.kmax > 0 ? P.kmax : std::max(I.ffd_vehicles(), kmin) + 1;
    // time windows can force many more vehicles than capacity does: use the
    // stated fleet (Solomon: VEHICLES) or no cap at all
    if (I.has_tw && P.kmax <= 0) kmax = I.vehicles > 0 ? I.vehicles : I.n;
    if (I.vehicles > 0) {
        kmax = std::max(std::min(kmax, I.vehicles), kmin);
        // CVRPLIB CVRP sets fix the fleet exactly; Solomon's VEHICLES is an upper bound
        if (P.fleet_mode == 0 && !I.has_tw && I.vehicles >= kmin) kmin = kmax = I.vehicles;
    }
    R.kmin = kmin;
    R.kmax = kmax;
    const int nv = I.nv();

    Master M(I, kmin, kmax, P);
    Solution savings;
    {
        std::vector<Route> init;
        for (int i = 1; i <= I.n; i++) {
            Route r;
            r.seq = {i};
            r.cost = I.d(0, i) + I.d(i, 0);
            init.push_back(r);
        }
        if (P.savings_init) {
            savings = clarke_wright(I);
            if (savings.cost >= 0)
                for (const auto &seq : savings.routes) {
                    Route r;
                    r.seq = seq;
                    int prev = 0;
                    for (int v : seq) { r.cost += I.d(prev, v); prev = v; }
                    r.cost += I.d(prev, 0);
                    init.push_back(r);
                }
        }
        M.add_routes(init);
    }
    Pricer pricer(I, P);
    RccSeparator rcc(I);

    double ub = P.ub; // exclusive: we look for cost < ub
    Solution best;
    auto accept = [&](const Solution &S, const char *src) {
        if (S.cost < 0 || (double)S.cost >= ub) return;
        CheckResult c = check_solution(I, S);
        if (!c.feasible || c.cost != S.cost) {
            std::printf("WARNING: %s solution rejected by the checker: %s\n", src, c.error.c_str());
            return;
        }
        ub = (double)S.cost;
        best = S;
        if (P.log_level >= 1) std::printf("* incumbent %lld from %s (%.1fs)\n", S.cost, src, now() - t0);
    };

    if (savings.cost >= 0) {
        int k = (int)savings.routes.size();
        if (k >= kmin && k <= kmax) accept(savings, "Clarke-Wright savings");
    }

    // open nodes: min-heap on lb (then older first)
    std::vector<Node> nodes;
    std::vector<int> open;
    auto cmp = [&](int a, int b) {
        if (nodes[a].lb != nodes[b].lb) return nodes[a].lb > nodes[b].lb;
        return nodes[a].id > nodes[b].id;
    };
    auto push_open = [&](int i) { open.push_back(i); std::push_heap(open.begin(), open.end(), cmp); };
    auto pop_open = [&]() { std::pop_heap(open.begin(), open.end(), cmp); int i = open.back(); open.pop_back(); return i; };

    nodes.push_back(Node{0, -1, 0, 0.0, {}, kmin, kmax});
    push_open(0);
    std::vector<int> active_rows; // branch row ids currently active in the master
    std::vector<std::pair<int, int>> forbidden;

    auto activate = [&](const Node &N) {
        double t = now();
        for (int id : active_rows)
            if (std::find(N.rows.begin(), N.rows.end(), id) == N.rows.end()) M.set_row_active(id, false);
        for (int id : N.rows) M.set_row_active(id, true); // no-op when already active
        active_rows = N.rows;
        forbidden.clear();
        for (int id : N.rows) {
            const RobustRow &row = M.row_by_id(id);
            if (row.up <= 0.5) forbidden.push_back({row.u, row.v});
        }
        pricer.set_forbidden(forbidden);
        M.set_fleet_bounds(N.klo, N.khi);
        R.t_switch += now() - t;
    };

    // Column-fixing dive: fix the route with the largest fractional lambda,
    // drop its customers from pricing, re-run CG (no cuts), repeat until the
    // LP is integral or infeasible. Everything is undone afterwards; the
    // columns generated stay in the pool.
    auto run_dive = [&](const Node &N, double tl) {
        double t = now();
        Params DP = P;
        DP.rcc = false;
        DP.src = false;
        DP.early_stop = false;
        DP.log_level = P.log_level >= 3 ? 3 : 0;
        std::vector<int> fixed;
        std::vector<char> removed(I.nv(), 0);
        bool found = false;
        for (int step = 0; step < I.n; step++) {
            if (now() - t > tl) break;
            NodeResult r = solve_node(I, DP, M, pricer, nullptr, tl - (now() - t), ub < 1e299 ? ub : 1e300);
            if (!r.ok || r.infeasible) break;
            if (r.integral) { accept(r.sol, "dive"); found = true; break; }
            if (ub < 1e299 && std::ceil(r.lb - 1e-6) >= ub) break;
            const auto &lam = M.lambda();
            int best = -1;
            double bl = 0;
            for (size_t k = 0; k < lam.size(); k++) {
                if (lam[k] <= 1e-6 || lam[k] >= 1 - 1e-6) continue;
                if (!M.route((int)k).elementary()) continue;
                if (lam[k] > bl) { bl = lam[k]; best = (int)k; }
            }
            if (best < 0) break;
            M.fix_route(best, true);
            fixed.push_back(best);
            for (int v : M.route(best).seq) removed[v] = 1;
            pricer.set_removed(removed);
            if (P.log_level >= 3)
                std::printf("  dive step %d: fixed route %d (lambda %.3f, %zu customers), lp %.3f, %.1fs\n", step, best, bl,
                            M.route(best).seq.size(), r.lp, now() - t);
        }
        for (int r : fixed) M.fix_route(r, false);
        pricer.clear_removed();
        R.t_dive += now() - t;
        if (P.log_level >= 2)
            std::printf("  dive from node %d: %zu fixings, %s (%.1fs)\n", N.id, fixed.size(),
                        found ? "solution" : "no solution", now() - t);
    };

    int cols_at_last_rmip = 0;
    auto run_rmip = [&](const char *why, double tl) {
        if (!P.rmip) return;
        std::vector<int> used;
        double obj = 0;
        double t = now();
        double cutoff = ub < 1e299 ? ub - 0.5 : BPLEX_INF; // integer costs: strictly better
        int st = M.solve_restricted_mip(cutoff, tl, used, obj);
        cols_at_last_rmip = M.nlpcols();
        R.t_rmip += now() - t;
        if (P.log_level >= 2)
            std::printf("  rmip (%s): status %s obj %.1f in %.1fs\n", why, bplex_status_string(st), obj, now() - t);
        if (!used.empty()) {
            Solution S;
            for (int r : used) S.routes.push_back(M.route(r).seq);
            S.cost = (long long)std::llround(obj);
            accept(S, "restricted master MIP");
        }
    };

    // delete branch rows no open node refers to, and cut rows whose dual has
    // been zero for many solves (they can be separated again later)
    auto cleanup_rows = [&](int cur) {
        // ids are never reused: index the reference marks by the id bound,
        // not by the row count (a branch row of an open node must survive)
        std::vector<char> ref(M.row_id_bound(), 0);
        auto mark = [&](const Node &N) { for (int id : N.rows) ref[id] = 1; };
        for (int i : open) mark(nodes[i]);
        if (cur >= 0) mark(nodes[cur]);
        std::vector<int> del;
        size_t ncut = 0;
        for (const RobustRow &r : M.rows()) {
            if (r.kind == RobustRow::EDGE) {
                if (!r.active && !ref[r.id]) del.push_back(r.id);
            } else if (r.zero_streak >= P.cut_age_limit) {
                del.push_back(r.id);
                ncut++;
            }
        }
        if (!del.empty()) M.delete_rows(del);
        if (P.log_level >= 2)
            std::printf("  cleanup: deleted %zu rows (%zu aged cuts), %zu left\n", del.size(), ncut, M.rows().size());
    };

    // pseudocosts per edge and direction (gain per unit of fractionality)
    std::vector<PseudoCost> pc((size_t)nv * nv);
    double pc_avg_sum[2] = {0, 0};
    int pc_avg_cnt[2] = {0, 0};
    auto pc_update = [&](int e, int dir, double gain, double frac) {
        if (frac < 1e-6 || gain < 0) return;
        double g = std::min(gain, 1e6) / frac;
        pc[e].sum[dir] += g;
        pc[e].cnt[dir]++;
        pc_avg_sum[dir] += g;
        pc_avg_cnt[dir]++;
    };
    auto pc_est = [&](int e, int dir, double frac) {
        double g = pc[e].cnt[dir] ? pc[e].sum[dir] / pc[e].cnt[dir]
                                  : (pc_avg_cnt[dir] ? pc_avg_sum[dir] / pc_avg_cnt[dir] : 1.0);
        return g * frac;
    };

    int next_id = 1, last_rmip_node = 0, last_cleanup = 0, last_dive_node = -1000000;
    double enum_fail_gap_any = 1e300; // smallest gap at which an enumeration attempt failed anywhere
    int cur = -1; // node being plunged into (child kept out of the heap)
    double global_lb = 0;
    std::string stop;

    // Debug: CVRPBCP_DEBUGSOL=<file.sol> loads a known solution; every node
    // reports whether it is feasible there and whether enumeration kept all
    // of its routes.
    Solution dbg;
    bool have_dbg = false;
    if (const char *p = std::getenv("CVRPBCP_DEBUGSOL")) {
        try { dbg = read_sol(p, I.scale); have_dbg = check_solution(I, dbg).feasible; } catch (...) {}
        if (have_dbg) std::printf("debug: loaded solution %s cost %lld\n", p, check_solution(I, dbg).cost);
    }
    auto dbg_feasible_here = [&](const Node &N) {
        for (int id : N.rows) {
            const RobustRow &row = M.row_by_id(id);
            int cnt = 0;
            for (const auto &seq : dbg.routes) {
                Route r;
                r.seq = seq;
                cnt += row.coef(r);
            }
            if (cnt < row.lo - 1e-9 || cnt > row.up + 1e-9) return false;
        }
        int k = (int)dbg.routes.size();
        return k >= N.klo && k <= N.khi;
    };

    while (true) {
        if (cur < 0) {
            if (open.empty()) break;
            cur = pop_open();
        }
        Node N = nodes[cur];
        global_lb = N.lb;
        if (!open.empty()) global_lb = std::min(global_lb, nodes[open.front()].lb);
        if (ub < 1e299 && std::ceil(N.lb - 1e-6) >= ub) { cur = -1; continue; } // pruned by bound
        if (now() - t0 > P.time_limit) { stop = "time_limit"; break; }
        if (R.nodes >= P.node_limit) { stop = "node_limit"; break; }
        if (R.nodes - last_cleanup >= P.cleanup_every) { last_cleanup = R.nodes; cleanup_rows(cur); }

        activate(N);
        double budget = P.time_limit - (now() - t0);
        NodeResult res = solve_node(I, P, M, pricer, &rcc, budget, ub < 1e299 ? ub : 1e300, N.enum_set.get());
        R.nodes++;
        add_times(R, res);
        R.depth_max = std::max(R.depth_max, N.depth);
        if (!res.ok) {
            // a node stopped by the clock (its own check or an LP solve cut
            // off by the remaining budget) ends the run as a time limit
            if (res.error == "time limit" || res.error.find("time limit") != std::string::npos ||
                now() - t0 >= P.time_limit) {
                stop = "time_limit";
                break;
            }
            R.error = res.error;
            stop = "error";
            break;
        }
        if (N.id == 0) {
            R.root_lp = res.lp;
            R.root_lb = res.lb;
            R.time_root = now() - t0;
        }
        nodes[cur].lb = std::max(N.lb, res.lb);
        double nlb = nodes[cur].lb;
        if (N.bedge >= 0) // true pseudocost from the solved child
            pc_update(N.bedge, N.bdir, (res.infeasible ? 1e6 : res.lp) - N.parent_lp, N.bfrac);
        if (P.log_level >= 1 && (N.id == 0 || P.log_level >= 2 || R.nodes % 10 == 0))
            std::printf("node %5d d%-3d lb %10.3f lp %10.3f glb %10.3f ub %s open %zu lp %d/%d rows %d it %d %.1fs\n",
                        N.id, N.depth, res.lb, res.lp, global_lb, ub < 1e299 ? std::to_string((long long)ub).c_str() : "-",
                        open.size(), M.nlpcols(), res.ncols, res.ncuts, res.cg_iters, now() - t0);

        if (res.integral) accept(res.sol, "node LP");
        if (have_dbg && dbg_feasible_here(N)) {
            long long dc = check_solution(I, dbg).cost;
            if (res.infeasible || std::ceil(nlb - 1e-6) > dc)
                std::printf("  DEBUG node %d: known solution (cost %lld) feasible here but node %s (lb %.4f lp %.4f)\n",
                            N.id, dc, res.infeasible ? "declared infeasible" : "bound exceeds it", nlb, res.lp);
        }
        if (res.infeasible) { cur = -1; continue; }
        if (ub < 1e299 && std::ceil(nlb - 1e-6) >= ub) { cur = -1; continue; }
        if (res.integral && (double)res.sol.cost <= std::ceil(nlb - 1e-6) + 1e-9) { cur = -1; continue; }

        if (P.dive && (N.id == 0 || (ub >= 1e299 && R.nodes - last_dive_node >= P.dive_every))) {
            last_dive_node = R.nodes;
            run_dive(N, std::min(P.dive_time, P.time_limit - (now() - t0)));
            if (ub < 1e299 && std::ceil(nlb - 1e-6) >= ub) { cur = -1; continue; }
            // the dive changed the LP (new columns, bounds restored): re-solve
            // the node so the duals, bound and edge flows are current again
            NodeResult again = solve_node(I, P, M, pricer, &rcc, P.time_limit - (now() - t0),
                                          ub < 1e299 ? ub : 1e300, N.enum_set.get());
            add_times(R, again);
            if (again.ok) {
                res = again;
                nodes[cur].lb = std::max(nodes[cur].lb, res.lb);
                nlb = nodes[cur].lb;
                if (res.integral) accept(res.sol, "node LP");
                if (res.infeasible || (ub < 1e299 && std::ceil(nlb - 1e-6) >= ub)) { cur = -1; continue; }
            } else {
                res.lb_last = -1e300; // no fresh duals: skip enumeration at this node
            }
        }
        if (P.rmip) {
            bool grown = M.nlpcols() >= cols_at_last_rmip + std::max(50, cols_at_last_rmip / 20);
            double left = P.time_limit - (now() - t0);
            if (N.id == 0) run_rmip("root", std::min(P.rmip_time, left));
            else if ((ub >= 1e299 || R.nodes - last_rmip_node >= P.rmip_every) && grown) {
                last_rmip_node = R.nodes;
                run_rmip("periodic", std::min(P.rmip_time_node, left));
            }
            if (ub < 1e299 && std::ceil(nlb - 1e-6) >= ub) { cur = -1; continue; }
        }

        // ---- enumeration + MIP: every route of a solution cheaper than ub
        // has reduced cost <= ub - 1 - lb under the node's final duals ----
        // ---- enumeration: with duals pi of the last exact pricing, every route
        // of a solution cheaper than ub has rc(pi) <= ub - 1 - LB(pi) + min_rc(pi).
        // Small sets go to a set-partitioning MIP; larger ones replace labeling
        // in the subtree and shrink at each descendant. ----
        double enum_gap = (ub - 1.0) - res.lb_last + res.minrc_last + 1e-6;
        if (P.enumerate && ub < 1e299 && res.lb_last > -1e299 && now() - t0 < P.time_limit &&
            (N.enum_set || enum_gap < 0.5 * std::min(N.enum_fail_gap, enum_fail_gap_any))) {
            double te = now();
            double gap = enum_gap;
            std::shared_ptr<std::vector<Route>> set;
            bool ok;
            if (N.enum_set) {
                set = std::make_shared<std::vector<Route>>();
                double mr;
                // keep the routes still within the (smaller) gap; scan_routes with a
                // large max returns those with rc < -(-gap), i.e. rc <= gap
                for (const Route &r : *N.enum_set) {
                    int prev = 0;
                    bool bad = false;
                    for (int v : r.seq) { if (pricer.forbidden(prev, v)) { bad = true; break; } prev = v; }
                    if (bad || pricer.forbidden(prev, 0)) continue;
                    if (pricer.route_rc(r) <= gap) set->push_back(r);
                }
                (void)mr;
                ok = true;
            } else {
                set = std::make_shared<std::vector<Route>>();
                ok = pricer.enumerate(gap, P.enum_max_labels, P.enum_max_routes, *set, P.enum_time_cap);
                R.enum_tries++;
                if (!ok) { nodes[cur].enum_fail_gap = gap; enum_fail_gap_any = std::min(enum_fail_gap_any, gap); }
            }
            R.t_enum_label += now() - te;
            if (P.log_level >= 2)
                std::printf("  enumeration: gap %.3f -> %s%zu routes (%.1fs)%s\n", gap, ok ? "" : "FAILED, >",
                            set->size(), now() - te, N.enum_set ? " [filtered]" : "");
            if (have_dbg && ok && dbg_feasible_here(N) && check_solution(I, dbg).cost < ub) {
                for (const auto &seq : dbg.routes) {
                    // the set must hold this route or another order of the same
                    // customers that is no dearer (enumeration dominance keeps those)
                    Route r;
                    r.seq = seq;
                    int prev = 0;
                    for (int v : seq) { r.cost += I.d(prev, v); prev = v; }
                    r.cost += I.d(prev, 0);
                    std::vector<int> key(seq);
                    std::sort(key.begin(), key.end());
                    bool found = false;
                    for (const Route &o : *set) {
                        if (o.seq.size() != seq.size() || o.cost > r.cost) continue;
                        std::vector<int> ok(o.seq);
                        std::sort(ok.begin(), ok.end());
                        if (ok == key) { found = true; break; }
                    }
                    if (!found) {
                        std::printf("  DEBUG node %d: route of the known solution missing from the enumeration: rc %.4f gap %.4f "
                                    "(lb_last %.4f minrc_last %.4f lp %.4f ub %.0f)\n",
                                    N.id, pricer.route_rc(r), gap, res.lb_last, res.minrc_last, res.lp, ub);
                        if (!N.enum_set) pricer.debug_trace_route(r, gap);
                    }
                }
            }
            if (ok) {
                nodes[cur].enum_set = set;
                if (set->size() <= P.enum_mip_max_routes) {
                    std::string dump;
                    if (!P.enum_dump.empty())
                        dump = P.enum_dump + "_" + std::to_string(R.enum_closed + R.enum_tries) + ".mps";
                    double tm = now();
                    EnumMipResult E = solve_enum_mip(I, *set, M.rows(), N.klo, N.khi, ub - 0.5,
                                                     std::max(1.0, std::min(P.enum_mip_time, P.time_limit - (now() - t0))),
                                                     P.threads, P.log_level, dump);
                    R.t_enum_mip += now() - tm;
                    R.enum_mips++;
                    if (P.log_level >= 2)
                        std::printf("  enum MIP: %s, bound %.3f, %s (%.1fs)\n", bplex_status_string(E.status),
                                    E.bestbound, E.has_solution ? "solution" : "no solution", E.time);
                    if (E.has_solution) accept(E.sol, "enumeration MIP");
                    if (have_dbg && E.solved && dbg_feasible_here(N) && check_solution(I, dbg).cost < ub)
                        std::printf("  DEBUG node %d: enumeration MIP closed a node that contains the known solution "
                                    "(cost %lld < ub %.0f, status %s)\n",
                                    N.id, check_solution(I, dbg).cost, ub, bplex_status_string(E.status));
                    if (E.solved) {
                        R.enum_closed++;
                        R.t_enum += now() - te;
                        nodes[cur].lb = std::max(nodes[cur].lb, ub);
                        cur = -1;
                        continue;
                    }
                    if (E.bestbound > nlb) {
                        nodes[cur].lb = std::min(E.bestbound, ub);
                        nlb = nodes[cur].lb;
                    }
                }
            } else {
                nodes[cur].enum_set.reset();
            }
            R.t_enum += now() - te;
            if (ub < 1e299 && std::ceil(nlb - 1e-6) >= ub) { cur = -1; continue; }
        }

        // ---- branching ----
        double tb = now();
        Node c1, c2;
        c1.parent = c2.parent = N.id;
        c1.depth = c2.depth = N.depth + 1;
        c1.lb = c2.lb = nlb;
        c1.rows = c2.rows = N.rows;
        c1.klo = c2.klo = N.klo;
        c1.khi = c2.khi = N.khi;
        c1.parent_lp = c2.parent_lp = res.lp;
        c1.enum_set = c2.enum_set = nodes[cur].enum_set;
        c1.enum_fail_gap = c2.enum_fail_gap = nodes[cur].enum_fail_gap;
        bool branched = false;
        bool plunge_second = true;
        double f = res.fleet;
        if (f - std::floor(f) > 1e-6 && std::ceil(f) - f > 1e-6) {
            c1.khi = (int)std::floor(f);
            c2.klo = (int)std::ceil(f);
            branched = true;
            if (P.log_level >= 2)
                std::printf("  branch on fleet %.3f: <= %d | >= %d  (node [%d,%d] master [%d,%d] art %.6f)\n", f, c1.khi,
                            c2.klo, N.klo, N.khi, M.klo(), M.khi(), M.artificial_mass());
        } else {
            struct Cand { int u, v; double x, score; bool reliable; };
            std::vector<Cand> cands;
            for (int i = 0; i < nv; i++)
                for (int j = i + 1; j < nv; j++) {
                    double v = res.x[(size_t)i * nv + j];
                    double fr = v - std::floor(v);
                    if (fr < 1e-6 || fr > 1 - 1e-6) continue;
                    int e = i * nv + j;
                    double gd = pc_est(e, 0, fr), gu = pc_est(e, 1, 1 - fr);
                    bool rel = pc[e].cnt[0] >= P.sb_reliable && pc[e].cnt[1] >= P.sb_reliable;
                    cands.push_back({i, j, v, std::max(gd, 1e-4) * std::max(gu, 1e-4), rel});
                }
            if (!cands.empty()) {
                std::sort(cands.begin(), cands.end(), [](const Cand &a, const Cand &b) { return a.score > b.score; });
                int bi = 0;
                double bdown = -1, bup = -1;
                if (P.sb_cands > 0) {
                    // probe the best-scored unreliable candidates, with lookahead
                    std::vector<int> probe_idx;
                    for (int c = 0; c < (int)cands.size() && (int)probe_idx.size() < P.sb_cands; c++)
                        if (!cands[c].reliable) probe_idx.push_back(c);
                    std::vector<RobustRow> rs;
                    for (int c : probe_idx) {
                        RobustRow down, up;
                        down.kind = up.kind = RobustRow::EDGE;
                        down.u = up.u = cands[c].u;
                        down.v = up.v = cands[c].v;
                        down.lo = -BPLEX_INF;
                        down.up = std::floor(cands[c].x);
                        up.lo = std::ceil(cands[c].x);
                        up.up = BPLEX_INF;
                        rs.push_back(down);
                        rs.push_back(up);
                    }
                    std::vector<int> ids = M.add_rows(rs);
                    double bestscore = cands.empty() ? -1 : (probe_idx.empty() ? cands[0].score : -1);
                    int since_improve = 0;
                    std::vector<char> used(ids.size(), 0);
                    for (size_t k = 0; k < probe_idx.size(); k++) {
                        int c = probe_idx[k];
                        double d = M.probe(ids[2 * k], P.sb_iters);
                        double u = M.probe(ids[2 * k + 1], P.sb_iters);
                        double fr = cands[c].x - std::floor(cands[c].x);
                        int e = cands[c].u * nv + cands[c].v;
                        // a probe without information (missing columns) falls back to the pseudocost estimate
                        if (d >= 0) pc_update(e, 0, d - res.lp, fr); else d = res.lp + pc_est(e, 0, fr);
                        if (u >= 0) pc_update(e, 1, u - res.lp, 1 - fr); else u = res.lp + pc_est(e, 1, 1 - fr);
                        double gd = std::max(d - res.lp, 1e-4);
                        double gu = std::max(u - res.lp, 1e-4);
                        double s = gd * gu;
                        cands[c].score = s;
                        if (s > bestscore) { bestscore = s; bi = c; bdown = d; bup = u; since_improve = 0; }
                        else since_improve++;
                        if (since_improve >= P.sb_lookahead) break;
                    }
                    // the best may be a reliable (unprobed) candidate with a higher pseudocost score
                    for (int c = 0; c < (int)cands.size(); c++)
                        if (cands[c].reliable && cands[c].score > bestscore) { bestscore = cands[c].score; bi = c; bdown = bup = -1; }
                    // keep the chosen candidate's rows, delete the other probe rows
                    std::vector<int> del;
                    int keep_k = -1;
                    for (size_t k = 0; k < probe_idx.size(); k++)
                        if (probe_idx[k] == bi) keep_k = (int)k;
                    for (size_t k = 0; k < probe_idx.size(); k++)
                        if ((int)k != keep_k) { del.push_back(ids[2 * k]); del.push_back(ids[2 * k + 1]); }
                    M.delete_rows(del);
                    if (keep_k >= 0) {
                        c1.rows.push_back(ids[2 * keep_k]);
                        c2.rows.push_back(ids[2 * keep_k + 1]);
                    }
                }
                if (c1.rows.size() == N.rows.size()) {
                    RobustRow down, up;
                    down.kind = up.kind = RobustRow::EDGE;
                    down.u = up.u = cands[bi].u;
                    down.v = up.v = cands[bi].v;
                    down.lo = -BPLEX_INF;
                    down.up = std::floor(cands[bi].x);
                    up.lo = std::ceil(cands[bi].x);
                    up.up = BPLEX_INF;
                    std::vector<int> ids = M.add_rows({down, up});
                    c1.rows.push_back(ids[0]);
                    c2.rows.push_back(ids[1]);
                }
                double fr = cands[bi].x - std::floor(cands[bi].x);
                c1.bedge = c2.bedge = cands[bi].u * nv + cands[bi].v;
                c1.bdir = 0;
                c2.bdir = 1;
                c1.bfrac = fr;
                c2.bfrac = 1 - fr;
                if (bdown >= 0) plunge_second = bup <= bdown; // probe bounds are over current columns only: not valid node bounds
                branched = true;
                if (P.log_level >= 2)
                    std::printf("  branch on edge (%d,%d) = %.3f  probes %.3f | %.3f  score %.3f\n", cands[bi].u,
                                cands[bi].v, cands[bi].x, std::min(bdown, 1e9), std::min(bup, 1e9), cands[bi].score);
            }
        }
        R.t_branch += now() - tb;
        if (!branched) {
            std::printf("WARNING: node %d: nothing to branch on (lp %.3f)\n", N.id, res.lp);
            cur = -1;
            continue;
        }
        c1.id = next_id++;
        c2.id = next_id++;
        nodes.push_back(c1);
        nodes.push_back(c2);
        int i1 = (int)nodes.size() - 2, i2 = (int)nodes.size() - 1;
        int dive = plunge_second ? i2 : i1, other = plunge_second ? i1 : i2;
        // a fleet child with an empty vehicle range is infeasible: drop it
        auto feasible = [&](int i) { return nodes[i].klo <= nodes[i].khi; };
        if (feasible(other)) push_open(other);
        bool dive_ok = feasible(dive) && (N.depth + 1 <= P.plunge_depth) &&
                       !(ub < 1e299 && std::ceil(nodes[dive].lb - 1e-6) >= ub);
        if (dive_ok) cur = dive;
        else { if (feasible(dive)) push_open(dive); cur = -1; }
    }

    R.time = now() - t0;
    R.ncols = M.nroutes();
    R.nrows = (int)M.rows().size();
    if (stop.empty()) {
        if (best.cost >= 0) {
            R.status = "optimal";
            R.lb = best.cost;
        } else {
            R.status = "infeasible";
            R.lb = (long long)std::ceil(global_lb - 1e-6);
        }
    } else {
        R.status = stop;
        double lb = 1e300;
        if (cur >= 0) lb = std::min(lb, nodes[cur].lb);
        for (int i : open) lb = std::min(lb, nodes[i].lb);
        if (lb > 1e299) lb = best.cost >= 0 ? (double)best.cost : 0;
        R.lb = (long long)std::ceil(lb - 1e-6);
        if (best.cost >= 0 && R.lb >= best.cost) { R.status = "optimal"; R.lb = best.cost; }
    }
    if (best.cost >= 0) {
        R.ub = best.cost;
        R.sol = best;
        int used = 0;
        for (auto &r : best.routes) used += !r.empty();
        R.fleet_cap_binding = (used >= kmax) && (I.vehicles <= 0);
    }
    if (P.log_level >= 1) {
        std::printf("times: pricing %.1f master %.1f cuts %.1f rmip %.1f dive %.1f branch %.1f switch %.1f enum %.1f "
                    "(label %.1f in %d tries, mip %.1f in %d solves, %d nodes closed) total %.1f\n",
                    R.t_pricing, R.t_master, R.t_cuts, R.t_rmip, R.t_dive, R.t_branch, R.t_switch, R.t_enum,
                    R.t_enum_label, R.enum_tries, R.t_enum_mip, R.enum_mips, R.enum_closed, R.time);
        std::printf("split: cg iters %d; pricing setup %.1f pool %.1f/%d s0 %.1f/%d s1 %.1f/%d exact %.1f/%d scan %.1f/%d; "
                    "cuts rccsep %.1f/%d (%d with cuts) src %.1f rows %.1f\n",
                    R.cg_iters, R.t_pr[0], R.t_pr[1], R.n_pr[1], R.t_pr[2], R.n_pr[2], R.t_pr[3], R.n_pr[3], R.t_pr[4],
                    R.n_pr[4], R.t_pr[5], R.n_pr[5], R.t_sep, R.n_sep, R.n_sep_hit, R.t_src,
                    R.t_cuts - R.t_sep - R.t_src);
        const auto &ps = M.probe_stats;
        std::printf("probes: %d  solve %.1f (%.0f iters avg)  purged %ld repriced %ld lp cols %d pool %d\n", ps.n,
                    ps.t_solve, ps.n ? ps.iters / ps.n : 0.0, M.purged_total, M.repriced_total, M.nlpcols(),
                    M.nroutes());
    }
    return R;
}

} // namespace cvrpbcp
