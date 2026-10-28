#include "master/master.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "bplex.h"
#include "pricing/labeling.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

namespace cvrpbcp {

// ---- LP solver wrappers (bplex + optional trace) ----
int Master::bp_load(int ncols, int nrows, const double *obj, const double *lo, const double *up, const double *rlo,
                    const double *rup, const int *cs, const int *ri, const double *v)
{
    if (trace_) trace_->load(ncols, nrows, obj, lo, up, rlo, rup, cs, ri, v);
    return bplex_load(bp_, ncols, nrows, obj, lo, up, nullptr, rlo, rup, cs, ri, v, 0.0, 1);
}
int Master::bp_add_cols(int n, const double *obj, const double *lo, const double *up, const int *cs, const int *ri,
                        const double *v)
{
    if (trace_) trace_->add_cols(n, obj, lo, up, cs, ri, v);
    return bplex_add_cols(bp_, n, obj, lo, up, cs, ri, v);
}
int Master::bp_add_rows(int n, const double *lo, const double *up, const int *rs, const int *ci, const double *v)
{
    if (trace_) trace_->add_rows(n, lo, up, rs, ci, v);
    return bplex_add_rows(bp_, n, lo, up, rs, ci, v);
}
int Master::bp_del_cols(int n, const int *idx)
{
    if (trace_) trace_->del_cols(n, idx);
    return bplex_del_cols(bp_, n, idx);
}
int Master::bp_del_rows(int n, const int *idx)
{
    if (trace_) trace_->del_rows(n, idx);
    return bplex_del_rows(bp_, n, idx);
}
int Master::bp_chg_col(int j, double lo, double up)
{
    if (trace_) trace_->chg_col(j, lo, up);
    return bplex_chg_col_bounds(bp_, j, lo, up);
}
int Master::bp_chg_row(int i, double lo, double up)
{
    if (trace_) trace_->chg_row(i, lo, up);
    return bplex_chg_row_bounds(bp_, i, lo, up);
}
int Master::bp_chg_type(int j, int is_integer)
{
    if (trace_) trace_->chg_type(j, is_integer);
    return bplex_chg_col_type(bp_, j, is_integer);
}
void Master::bp_slack_basis()
{
    if (trace_) trace_->slack_basis();
    std::vector<int> cst(nfix_ + (int)lp_route_.size(), BPLEX_AT_LOWER), rst(nrows_, BPLEX_BASIC);
    bplex_set_basis(bp_, cst.data(), rst.data());
}
void Master::bp_param(const char *name, double value)
{
    if (trace_) trace_->param(name, value);
    bplex_set_param(bp_, name, value);
}
int Master::bp_solve(int kind)
{
    if (!trace_) return bplex_solve(bp_);
    using namespace std::chrono;
    trace_->solve(kind);
    double t = duration<double>(steady_clock::now().time_since_epoch()).count();
    int st = bplex_solve(bp_);
    t = duration<double>(steady_clock::now().time_since_epoch()).count() - t;
    double obj = (st == BPLEX_OPTIMAL || st == BPLEX_TIME_LIMIT || st == BPLEX_NODE_LIMIT) ? bplex_objval(bp_) : 0.0;
    trace_->result(st, obj, bplex_get_info(bp_, "iterations"), t);
    return st;
}

bool Route::elementary() const
{
    std::vector<int> s(seq);
    std::sort(s.begin(), s.end());
    return std::adjacent_find(s.begin(), s.end()) == s.end();
}

int RobustRow::arc_coef(int i, int j) const
{
    if (kind == RCC) return in_s[i] != in_s[j];
    if (kind == EDGE) return (std::min(i, j) == u && std::max(i, j) == v);
    return 0;
}

int RobustRow::coef(const Route &r) const
{
    if (kind == SRC) {
        int c = 0, state = 0;
        for (int x : r.seq) {
            if (!in_m[x]) { state = 0; continue; }
            if (in_s[x] && ++state == 2) { c++; state = 0; }
        }
        return c;
    }
    int c = 0, prev = 0;
    for (int x : r.seq) { c += arc_coef(prev, x); prev = x; }
    c += arc_coef(prev, 0);
    return c;
}

Master::Master(const Instance &I, int kmin, int kmax, const Params &P) : I_(I), P_(P), klo_(kmin), khi_(kmax)
{
    const int n = I.n;
    nfix_ = n + 1;
    int maxd = 0;
    for (int i = 1; i <= n; i++) maxd = std::max(maxd, I.d(0, i));
    bigM_ = 4.0 * (n + 1) * maxd + 1.0; // above the cost of any route solution

    bp_ = bplex_create();
    if (!bp_) throw std::runtime_error("bplex_create failed");
    if (const char *tp = std::getenv("CVRPBCP_LPTRACE")) trace_.reset(new LpTrace(tp));
    bp_param("log_level", 0);
    // bplex's own warnings and errors (level 0) go to our log with a prefix;
    // everything else only at log level >= 4
    bplex_set_log_callback(
        bp_,
        [](bplex_s *, int level, const char *msg, void *user) {
            const Params *p = static_cast<const Params *>(user);
            if (level == 0 || p->log_level >= 4) std::printf("bplex[%d]: %s%s", level, msg,
                                                            (msg[0] && msg[std::strlen(msg) - 1] == '\n') ? "" : "\n");
        },
        const_cast<Params *>(&P_));
    bp_param("threads", P.threads);
    bp_param("presolve", 0);

    nrows_ = n + 1;
    std::vector<double> obj(nfix_, bigM_), lo(nfix_, 0.0), up(nfix_, BPLEX_INF), rlo(nrows_, 1.0), rup(nrows_, 1.0);
    std::vector<int> cs(nfix_ + 1), ri(nfix_);
    std::vector<double> v(nfix_, 1.0);
    for (int j = 0; j < nfix_; j++) { cs[j] = j; ri[j] = j; } // column n -> fleet row n
    cs[nfix_] = nfix_;
    rlo[n] = kmin;
    rup[n] = kmax;
    if (bp_load(nfix_, nrows_, obj.data(), lo.data(), up.data(), rlo.data(), rup.data(), cs.data(), ri.data(),
                v.data()) != 0)
        throw std::runtime_error("bplex_load failed");
    pi_.assign(n + 1, 0.0);
}

Master::~Master()
{
    if (bp_) bplex_free(bp_);
}

void Master::append_lp_columns(const std::vector<int> &route_idx)
{
    if (route_idx.empty()) return;
    const int n = I_.n;
    std::vector<double> obj, vals;
    std::vector<int> cs, ri;
    std::vector<int> cnt(n + 1, 0);
    cs.push_back(0);
    for (int ridx : route_idx) {
        const Route &r = routes_[ridx];
        obj.push_back((double)r.cost);
        for (int v : r.seq) cnt[v]++;
        for (int i = 1; i <= n; i++)
            if (cnt[i]) { ri.push_back(i - 1); vals.push_back(cnt[i]); cnt[i] = 0; }
        ri.push_back(n);
        vals.push_back(1.0);
        for (const RobustRow &c : rows_) {
            int x = c.coef(r);
            if (x) { ri.push_back(c.row); vals.push_back(x); }
        }
        cs.push_back((int)ri.size());
        col_of_[ridx] = (int)lp_route_.size();
        lp_route_.push_back(ridx);
        lp_art_row_.push_back(-1);
        age_.push_back(0);
        nlp_routes_++;
    }
    int rc = bp_add_cols((int)route_idx.size(), obj.data(), nullptr, nullptr, cs.data(), ri.data(),
                            vals.data());
    if (rc != 0) throw std::runtime_error("bplex_add_cols failed: " + std::to_string(rc));
}

int Master::add_routes(const std::vector<Route> &rs)
{
    if (rs.empty()) return 0;
    std::vector<int> idx;
    for (const Route &r : rs) {
        // no duplicates: a route already in the pool is re-added to the LP
        // only if it is out of it
        uint64_t h = 1469598103934665603ull;
        for (int v : r.seq) h = (h ^ (uint64_t)v) * 1099511628211ull;
        auto range = route_index_.equal_range(h);
        int found = -1;
        for (auto it = range.first; it != range.second; ++it)
            if (routes_[it->second].seq == r.seq) { found = it->second; break; }
        if (found >= 0) {
            // -2: queued in this batch already (a route listed twice must
            // get one column, or the pool-indexed lambda loses its value)
            if (col_of_[found] == -1) { idx.push_back(found); col_of_[found] = -2; }
            continue;
        }
        route_index_.emplace(h, (int)routes_.size());
        idx.push_back((int)routes_.size());
        routes_.push_back(r);
        pool_seq_.insert(pool_seq_.end(), r.seq.begin(), r.seq.end());
        pool_off_.push_back((int)pool_seq_.size());
        col_of_.push_back(-2);
        lambda_.push_back(0.0);
        pool_skip_.push_back(0);
    }
    append_lp_columns(idx);
    return (int)idx.size();
}

int Master::price_pool(const Pricer &pr, double eps, int max_cols, bool fresh)
{
    struct C { double rc; int r; };
    std::vector<C> cand;
    if (fresh) std::fill(pool_skip_.begin(), pool_skip_.end(), 0);
    // dual scale for the skip rule: the average customer dual
    double dscale = 0;
    for (int i = 1; i <= I_.n; i++) dscale += std::fabs(pi_[i]);
    dscale = std::max(dscale / I_.n, 1.0);
    // active SRC rows with negative duals only add (-y * coef >= 0) to a
    // route's reduced cost, so the arc part alone decides whether to look
    std::vector<int> srcrows;
    for (size_t c = 0; c < rows_.size(); c++)
        if (rows_[c].kind == RobustRow::SRC && rows_[c].active && yrow_[c] < -1e-9) srcrows.push_back((int)c);
    // flat sequences, one route per iteration, no shared writes besides the
    // route's own skip counter: the scan splits over threads trivially
    const int nr = (int)routes_.size();
    const int T = std::max(1, scan_threads);
    std::vector<std::vector<C>> cands(T);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T) if (T > 1)
#endif
    for (int r = 0; r < nr; r++) {
        if (col_of_[r] >= 0) continue;
        if (pool_skip_[r] > 0) { pool_skip_[r]--; continue; }
        const int *sq = &pool_seq_[pool_off_[r]];
        const int len = pool_off_[r + 1] - pool_off_[r];
        double rc = 0;
        int prev = 0;
        bool bad = false;
        for (int k = 0; k < len; k++) {
            const int v = sq[k];
            if (pr.forbidden(prev, v)) { bad = true; break; }
            rc += pr.arc_rc(prev, v);
            prev = v;
        }
        if (bad || pr.forbidden(prev, 0)) continue;
        rc += pr.arc_rc(prev, 0);
        if (rc >= -eps) {
            pool_skip_[r] = rc < dscale ? 2 : rc < 4 * dscale ? 6 : 15;
            continue;
        }
        for (int c : srcrows) rc -= yrow_[c] * rows_[c].coef(routes_[r]);
        if (rc < -eps) {
#ifdef _OPENMP
            cands[omp_get_thread_num()].push_back({rc, r});
#else
            cands[0].push_back({rc, r});
#endif
        }
    }
    for (int t = 0; t < T; t++) cand.insert(cand.end(), cands[t].begin(), cands[t].end());
    if (cand.empty()) return 0;
    std::sort(cand.begin(), cand.end(), [](const C &a, const C &b) { return a.rc < b.rc; });
    if ((int)cand.size() > max_cols) cand.resize(max_cols);
    std::vector<int> idx;
    for (const C &c : cand) idx.push_back(c.r);
    append_lp_columns(idx);
    repriced_total += (long)idx.size();
    return (int)idx.size();
}

void Master::delete_dynamic_columns(const std::vector<int> &js)
{
    if (js.empty()) return;
    std::vector<int> cols;
    for (int j : js) cols.push_back(nfix_ + j);
    int rc = bp_del_cols((int)cols.size(), cols.data());
    if (rc != 0) throw std::runtime_error("bplex_del_cols failed: " + std::to_string(rc));
    std::vector<char> gone(lp_route_.size(), 0);
    for (int j : js) gone[j] = 1;
    std::vector<int> lp2, art2, age2;
    nlp_routes_ = 0;
    for (size_t j = 0; j < lp_route_.size(); j++) {
        if (gone[j]) {
            if (lp_route_[j] >= 0) col_of_[lp_route_[j]] = -1;
            continue;
        }
        if (lp_route_[j] >= 0) { col_of_[lp_route_[j]] = (int)lp2.size(); nlp_routes_++; }
        lp2.push_back(lp_route_[j]);
        art2.push_back(lp_art_row_[j]);
        age2.push_back(age_[j]);
    }
    lp_route_.swap(lp2);
    lp_art_row_.swap(art2);
    age_.swap(age2);
    redcost_.clear();
}

int Master::purge(int age, int min_cols)
{
    if (nlp_routes_ <= min_cols || redcost_.empty()) return 0;
    std::vector<int> del;
    for (size_t j = 0; j < lp_route_.size(); j++)
        if (lp_route_[j] >= 0 && age_[j] >= age && lambda_[lp_route_[j]] <= 1e-9) del.push_back((int)j);
    if (nlp_routes_ - (int)del.size() < min_cols) del.resize(nlp_routes_ - min_cols);
    if (del.empty()) return 0;
    delete_dynamic_columns(del);
    purged_total += (long)del.size();
    return (int)del.size();
}

std::vector<int> Master::add_rows(const std::vector<RobustRow> &rs)
{
    std::vector<int> ids;
    if (rs.empty()) return ids;
    std::vector<int> rstart = {0}, ci;
    std::vector<double> vals, lo, up;
    for (RobustRow r : rs) {
        r.id = (int)pos_.size();
        r.row = nrows_ + (int)ids.size();
        for (size_t j = 0; j < lp_route_.size(); j++) {
            if (lp_route_[j] < 0) continue;
            int x = r.coef(routes_[lp_route_[j]]);
            if (x) { ci.push_back(nfix_ + (int)j); vals.push_back(x); }
        }
        rstart.push_back((int)ci.size());
        lo.push_back(r.active ? r.lo : -BPLEX_INF);
        up.push_back(r.active ? r.up : BPLEX_INF);
        pos_.push_back((int)rows_.size());
        rows_.push_back(std::move(r));
        yrow_.push_back(0.0);
        ids.push_back((int)pos_.size() - 1);
    }
    int rc = bp_add_rows((int)rs.size(), lo.data(), up.data(), rstart.data(), ci.data(), vals.data());
    if (rc != 0) throw std::runtime_error("bplex_add_rows failed: " + std::to_string(rc));
    nrows_ += (int)rs.size();
    // one artificial column per ">=" row
    std::vector<double> aobj, avals;
    std::vector<int> acs = {0}, ari;
    for (size_t k = 0; k < rs.size(); k++) {
        const RobustRow &r = rows_[pos_[ids[k]]];
        if (r.lo <= -1e20) continue;
        aobj.push_back(bigM_);
        ari.push_back(r.row);
        avals.push_back(1.0);
        acs.push_back((int)ari.size());
        lp_route_.push_back(-1);
        lp_art_row_.push_back(r.id);
        age_.push_back(INT_MIN);
    }
    if (!aobj.empty()) {
        rc = bp_add_cols((int)aobj.size(), aobj.data(), nullptr, nullptr, acs.data(), ari.data(), avals.data());
        if (rc != 0) throw std::runtime_error("bplex_add_cols (artificial) failed: " + std::to_string(rc));
    }
    return ids;
}

void Master::set_row_active(int id, bool active)
{
    if (!row_exists(id)) throw std::runtime_error("set_row_active: row id " + std::to_string(id) + " was deleted");
    RobustRow &r = rows_[pos_[id]];
    if (r.active == active) return;
    r.active = active;
    int rc = active ? bp_chg_row(r.row, r.lo, r.up)
                    : bp_chg_row(r.row, -BPLEX_INF, BPLEX_INF);
    if (rc != 0) row_edit_failed("set_row_active", rc, r.row);
}

void Master::row_edit_failed(const char *what, int rc, int row) const
{
    std::printf("bplex row edit failed: %s rc %d row %d of %d (ours) status %d; dumping cvrpbcp_rowedit_fail.mps\n",
                what, rc, row, nrows_, bplex_status(bp_));
    bplex_write_mps(bp_, "cvrpbcp_rowedit_fail.mps");
    throw std::runtime_error(std::string("bplex_chg_row_bounds failed in ") + what + " rc " + std::to_string(rc));
}

void Master::delete_rows(const std::vector<int> &ids)
{
    if (ids.empty()) return;
    std::vector<int> bprows;
    std::vector<char> del(rows_.size(), 0), delid(pos_.size(), 0);
    for (int id : ids) {
        if (!row_exists(id)) continue;
        int p = pos_[id];
        if (del[p]) continue;
        del[p] = 1;
        delid[id] = 1;
        bprows.push_back(rows_[p].row);
        pos_[id] = -1;
    }
    if (bprows.empty()) return;
    // their artificial columns first (independent indices)
    std::vector<int> js;
    for (size_t j = 0; j < lp_route_.size(); j++)
        if (lp_art_row_[j] >= 0 && delid[lp_art_row_[j]]) js.push_back((int)j);
    delete_dynamic_columns(js);
    int rc = bp_del_rows((int)bprows.size(), bprows.data());
    if (rc != 0) throw std::runtime_error("bplex_del_rows failed: " + std::to_string(rc));
    std::vector<RobustRow> kept;
    std::vector<double> ykept;
    kept.reserve(rows_.size());
    for (size_t p = 0; p < rows_.size(); p++) {
        if (del[p]) continue;
        RobustRow r = std::move(rows_[p]);
        r.row = I_.n + 1 + (int)kept.size();
        pos_[r.id] = (int)kept.size();
        kept.push_back(std::move(r));
        ykept.push_back(yrow_[p]);
    }
    rows_.swap(kept);
    yrow_.swap(ykept);
    nrows_ = I_.n + 1 + (int)rows_.size();
}

void Master::fix_route(int r, bool fixed)
{
    if (col_of_[r] < 0) append_lp_columns({r});
    int j = nfix_ + col_of_[r];
    int rc = fixed ? bp_chg_col(j, 1.0, 1.0) : bp_chg_col(j, 0.0, BPLEX_INF);
    if (rc != 0) throw std::runtime_error("bplex_chg_col_bounds failed");
    if (fixed) age_[col_of_[r]] = INT_MIN; // never purged while fixed
    else age_[col_of_[r]] = 0;
}

void Master::set_fleet_bounds(int klo, int khi)
{
    if (klo == klo_ && khi == khi_) return;
    if (klo > khi) throw std::runtime_error("set_fleet_bounds: klo " + std::to_string(klo) + " > khi " + std::to_string(khi));
    klo_ = klo;
    khi_ = khi;
    int rc = bp_chg_row(I_.n, klo, khi);
    if (rc != 0) row_edit_failed("set_fleet_bounds", rc, I_.n);
}

int Master::solve()
{
    using namespace std::chrono;
    auto tnow = []() { return duration<double>(steady_clock::now().time_since_epoch()).count(); };
    // A warm solve that takes far longer than the running average is a
    // simplex stall (seen on XML100_1113_07: one 2100 x 254 LP ran for the
    // rest of the budget). It is cut off and retried with the primal
    // simplex; the first two cases are dumped for the bplex team (the LP
    // trace, CVRPBCP_LPTRACE, reproduces the warm start itself).
    const double avg = lp_solves_ > 20 ? lp_time_sum_ / (double)lp_solves_ : 0.0;
    const double cap = (avg > 0 && P_.lp_stall_factor > 0) ? std::max(P_.lp_stall_min, P_.lp_stall_factor * avg)
                                                            : lp_time_limit_;
    const double t0 = tnow();
    bp_param("time_limit", std::min(lp_time_limit_, cap));
    int st = bp_solve(0);
    if (st == BPLEX_TIME_LIMIT && cap < lp_time_limit_) {
        // Restart the dual simplex from the slack basis with the remaining
        // budget. (A retry with the primal simplex was tried first and ran
        // past its time limit, see docs/BPLEX_API_REQUESTS.md.)
        lp_stalls++;
        const double used = tnow() - t0;
        if (lp_stalls <= 2) write_mps("cvrpbcp_lp_stall" + std::to_string(lp_stalls) + ".mps");
        std::printf("  LP stall: warm solve cut after %.1fs (average %.3fs, %d cols, %d rows); restarting from the slack basis\n",
                    used, avg, nfix_ + (int)lp_route_.size(), nrows_);
        std::fflush(stdout);
        bp_slack_basis();
        bp_param("time_limit", std::max(1.0, lp_time_limit_ - used));
        st = bp_solve(0);
    }
    lp_time_sum_ += tnow() - t0;
    lp_solves_++;
    last_iters_ = (int)bplex_get_info(bp_, "iterations");
    if (st != BPLEX_OPTIMAL) return st;
    const int n = I_.n;
    const int ncols = nfix_ + (int)lp_route_.size();
    std::vector<double> x(ncols), y(nrows_), act(nrows_);
    redcost_.assign(ncols, 0.0);
    if (bplex_get_solution(bp_, x.data(), act.data(), y.data(), redcost_.data()) != 0)
        throw std::runtime_error("bplex_get_solution failed after OPTIMAL");
    // invariants: one LP column per pool route (the tree reads the
    // pool-indexed lambda, so a duplicate column would lose its value) and
    // the fleet row satisfied by the LP optimum
    double sl = 0;
    for (size_t j = 0; j < lp_route_.size(); j++) if (lp_route_[j] >= 0) sl += x[nfix_ + j];
    {
        std::vector<char> seen(routes_.size(), 0);
        for (int r : lp_route_)
            if (r >= 0 && seen[r]++) throw std::runtime_error("master: pool route " + std::to_string(r) + " has two LP columns");
    }
    if (act[n] < klo_ - 1e-6 || act[n] > khi_ + 1e-6 || std::fabs(act[n] - sl - x[n]) > 1e-6) {
        // the fleet row is not satisfied by the reported optimum: our bounds
        // or bplex's warm start are wrong; dump and stop
        std::printf("FLEET ROW VIOLATED: activity %.6f bounds [%d, %d] sum lambda %.6f artificial x[n] %.6f iters %d; "
                    "dumping cvrpbcp_fleet_fail.mps\n", act[n], klo_, khi_, sl, x[n], last_iters_);
        bplex_write_mps(bp_, "cvrpbcp_fleet_fail.mps");
        throw std::runtime_error("fleet row violated after an OPTIMAL solve");
    }
    obj_ = bplex_objval(bp_);
    art_mass_ = 0;
    for (int j = 0; j < nfix_; j++) art_mass_ += x[j];
    std::fill(lambda_.begin(), lambda_.end(), 0.0);
    for (size_t j = 0; j < lp_route_.size(); j++) {
        double xj = x[nfix_ + j];
        if (lp_route_[j] < 0) { art_mass_ += xj; continue; }
        lambda_[lp_route_[j]] = xj;
        if (redcost_[nfix_ + j] > 1e-6 && xj <= 1e-9) age_[j]++;
        else age_[j] = 0;
    }
    for (int i = 1; i <= n; i++) pi_[i] = y[i - 1];
    pi_[0] = y[n];
    for (size_t c = 0; c < rows_.size(); c++) {
        yrow_[c] = rows_[c].active ? y[rows_[c].row] : 0.0;
        if (std::fabs(yrow_[c]) < 1e-9) rows_[c].zero_streak++;
        else rows_[c].zero_streak = 0;
    }
    return st;
}

double Master::lagrangian_bound(const std::vector<double> &pi, const std::vector<double> &yrow, double min_rc) const
{
    double lb = 0;
    for (int i = 1; i <= I_.n; i++) lb += pi[i];
    for (size_t c = 0; c < rows_.size(); c++) {
        if (!rows_[c].active || yrow[c] == 0) continue;
        double b = yrow[c] > 0 ? rows_[c].lo : rows_[c].up;
        if (std::fabs(b) >= 1e20) continue;
        lb += yrow[c] * b;
    }
    double m = min_rc + pi[0]; // reduced cost without the fleet dual
    lb += (m < 0) ? khi_ * m : klo_ * m;
    return lb;
}

void Master::edge_flows(std::vector<double> &x) const
{
    const int nv = I_.nv();
    x.assign((size_t)nv * nv, 0.0);
    for (int r : lp_route_) {
        if (r < 0) continue;
        double l = lambda_[r];
        if (l < 1e-9) continue;
        int prev = 0;
        for (int v : routes_[r].seq) {
            x[(size_t)std::min(prev, v) * nv + std::max(prev, v)] += l;
            prev = v;
        }
        x[(size_t)0 * nv + prev] += l;
    }
    for (int i = 0; i < nv; i++)
        for (int j = i + 1; j < nv; j++) x[(size_t)j * nv + i] = x[(size_t)i * nv + j];
}

int Master::solve_restricted_mip(double cutoff, double time_limit, std::vector<int> &routes_used, double &obj)
{
    const int ncols = nfix_ + (int)lp_route_.size();
    routes_used.clear();
    for (int j = 0; j < nfix_; j++) bp_chg_col(j, 0.0, 0.0); // no artificials in the MIP
    for (size_t j = 0; j < lp_route_.size(); j++) {
        if (lp_route_[j] < 0) bp_chg_col(nfix_ + (int)j, 0.0, 0.0);
        else bp_chg_type(nfix_ + (int)j, 1);
    }
    bp_param("relax", 0);
    bp_param("time_limit", time_limit);
    bp_param("cutoff", cutoff);
    bp_param("mip_abs_gap", 0.999); // integer costs
    bp_param("mip_rel_gap", 0.0);
    int st = bp_solve(1);
    std::vector<double> x(ncols);
    if ((st == BPLEX_OPTIMAL || st == BPLEX_TIME_LIMIT || st == BPLEX_NODE_LIMIT) &&
        bplex_get_solution(bp_, x.data(), nullptr, nullptr, nullptr) == 0) {
        obj = bplex_objval(bp_);
        for (size_t j = 0; j < lp_route_.size(); j++)
            if (lp_route_[j] >= 0 && x[nfix_ + j] > 0.5) routes_used.push_back(lp_route_[j]);
    }
    bp_param("relax", 1);
    bp_param("cutoff", BPLEX_INF);
    bp_param("time_limit", lp_time_limit_);
    for (int j = 0; j < nfix_; j++) bp_chg_col(j, 0.0, BPLEX_INF);
    for (size_t j = 0; j < lp_route_.size(); j++) {
        if (lp_route_[j] < 0) bp_chg_col(nfix_ + (int)j, 0.0, BPLEX_INF);
        else bp_chg_type(nfix_ + (int)j, 0);
    }
    return st;
}

double Master::probe(int id, int iter_limit)
{
    using namespace std::chrono;
    auto tnow = []() { return duration<double>(steady_clock::now().time_since_epoch()).count(); };
    set_row_active(id, true);
    double old_limit = bplex_get_param(bp_, "iter_limit");
    bp_param("iter_limit", iter_limit);
    bp_param("method", BPLEX_METHOD_DUAL);
    double t = tnow();
    int st = bp_solve(2);
    probe_stats.n++;
    probe_stats.t_solve += tnow() - t;
    probe_stats.iters += bplex_get_info(bp_, "iterations");
    double val = -1;
    if (st == BPLEX_OPTIMAL || st == BPLEX_ITER_LIMIT || st == BPLEX_CUTOFF) {
        val = (st == BPLEX_OPTIMAL) ? bplex_objval(bp_) : bplex_bestbound(bp_);
        // artificials in the iterate: the bound is bigM-inflated, no information
        std::vector<double> x(nfix_ + lp_route_.size());
        if (bplex_get_solution(bp_, x.data(), nullptr, nullptr, nullptr) == 0) {
            double a = 0;
            for (int j = 0; j < nfix_; j++) a += x[j];
            for (size_t j = 0; j < lp_route_.size(); j++)
                if (lp_route_[j] < 0) a += x[nfix_ + j];
            if (a > 1e-6 || val > 0.5 * bigM_) val = -1;
        }
    }
    bp_param("iter_limit", old_limit);
    bp_param("method", BPLEX_METHOD_AUTO);
    set_row_active(id, false);
    return val;
}

void Master::write_mps(const std::string &path) const { bplex_write_mps(bp_, path.c_str()); }

} // namespace cvrpbcp
