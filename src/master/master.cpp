#include "master/master.h"

#include "bplex.h"
#include "pricing/labeling.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace cvrpbcp {

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
    bplex_set_param(bp_, "log_level", 0);
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
    bplex_set_param(bp_, "threads", P.threads);
    bplex_set_param(bp_, "presolve", 0);

    nrows_ = n + 1;
    std::vector<double> obj(nfix_, bigM_), lo(nfix_, 0.0), up(nfix_, BPLEX_INF), rlo(nrows_, 1.0), rup(nrows_, 1.0);
    std::vector<int> cs(nfix_ + 1), ri(nfix_);
    std::vector<double> v(nfix_, 1.0);
    for (int j = 0; j < nfix_; j++) { cs[j] = j; ri[j] = j; } // column n -> fleet row n
    cs[nfix_] = nfix_;
    rlo[n] = kmin;
    rup[n] = kmax;
    if (bplex_load(bp_, nfix_, nrows_, obj.data(), lo.data(), up.data(), nullptr, rlo.data(), rup.data(), cs.data(),
                   ri.data(), v.data(), 0.0, 1) != 0)
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
    int rc = bplex_add_cols(bp_, (int)route_idx.size(), obj.data(), nullptr, nullptr, cs.data(), ri.data(),
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
            if (col_of_[found] < 0) idx.push_back(found);
            continue;
        }
        route_index_.emplace(h, (int)routes_.size());
        idx.push_back((int)routes_.size());
        routes_.push_back(r);
        col_of_.push_back(-1);
        lambda_.push_back(0.0);
    }
    append_lp_columns(idx);
    return (int)idx.size();
}

int Master::price_pool(const Pricer &pr, double eps, int max_cols)
{
    struct C { double rc; int r; };
    std::vector<C> cand;
    // active SRC rows with negative duals only add (-y * coef >= 0) to a
    // route's reduced cost, so the arc part alone decides whether to look
    std::vector<int> srcrows;
    for (size_t c = 0; c < rows_.size(); c++)
        if (rows_[c].kind == RobustRow::SRC && rows_[c].active && yrow_[c] < -1e-9) srcrows.push_back((int)c);
    for (size_t r = 0; r < routes_.size(); r++) {
        if (col_of_[r] >= 0) continue;
        const Route &rt = routes_[r];
        double rc = 0;
        int prev = 0;
        bool bad = false;
        for (int v : rt.seq) {
            if (pr.forbidden(prev, v)) { bad = true; break; }
            rc += pr.arc_rc(prev, v);
            prev = v;
        }
        if (bad || pr.forbidden(prev, 0)) continue;
        rc += pr.arc_rc(prev, 0);
        if (rc >= -eps) continue;
        for (int c : srcrows) rc -= yrow_[c] * rows_[c].coef(rt);
        if (rc < -eps) cand.push_back({rc, (int)r});
    }
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
    int rc = bplex_del_cols(bp_, (int)cols.size(), cols.data());
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
    int rc = bplex_add_rows(bp_, (int)rs.size(), lo.data(), up.data(), rstart.data(), ci.data(), vals.data());
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
        rc = bplex_add_cols(bp_, (int)aobj.size(), aobj.data(), nullptr, nullptr, acs.data(), ari.data(), avals.data());
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
    int rc = active ? bplex_chg_row_bounds(bp_, r.row, r.lo, r.up)
                    : bplex_chg_row_bounds(bp_, r.row, -BPLEX_INF, BPLEX_INF);
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
    int rc = bplex_del_rows(bp_, (int)bprows.size(), bprows.data());
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
    int rc = fixed ? bplex_chg_col_bounds(bp_, j, 1.0, 1.0) : bplex_chg_col_bounds(bp_, j, 0.0, BPLEX_INF);
    if (rc != 0) throw std::runtime_error("bplex_chg_col_bounds failed");
    if (fixed) age_[col_of_[r]] = INT_MIN; // never purged while fixed
    else age_[col_of_[r]] = 0;
}

void Master::set_fleet_bounds(int klo, int khi)
{
    if (klo == klo_ && khi == khi_) return;
    klo_ = klo;
    khi_ = khi;
    int rc = bplex_chg_row_bounds(bp_, I_.n, klo, khi);
    if (rc != 0) row_edit_failed("set_fleet_bounds", rc, I_.n);
}

int Master::solve()
{
    bplex_set_param(bp_, "time_limit", lp_time_limit_);
    int st = bplex_solve(bp_);
    last_iters_ = (int)bplex_get_info(bp_, "iterations");
    if (st != BPLEX_OPTIMAL) return st;
    const int n = I_.n;
    const int ncols = nfix_ + (int)lp_route_.size();
    std::vector<double> x(ncols), y(nrows_);
    redcost_.assign(ncols, 0.0);
    if (bplex_get_solution(bp_, x.data(), nullptr, y.data(), redcost_.data()) != 0)
        throw std::runtime_error("bplex_get_solution failed after OPTIMAL");
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
    for (int j = 0; j < nfix_; j++) bplex_chg_col_bounds(bp_, j, 0.0, 0.0); // no artificials in the MIP
    for (size_t j = 0; j < lp_route_.size(); j++) {
        if (lp_route_[j] < 0) bplex_chg_col_bounds(bp_, nfix_ + (int)j, 0.0, 0.0);
        else bplex_chg_col_type(bp_, nfix_ + (int)j, 1);
    }
    bplex_set_param(bp_, "relax", 0);
    bplex_set_param(bp_, "time_limit", time_limit);
    bplex_set_param(bp_, "cutoff", cutoff);
    bplex_set_param(bp_, "mip_abs_gap", 0.999); // integer costs
    bplex_set_param(bp_, "mip_rel_gap", 0.0);
    int st = bplex_solve(bp_);
    std::vector<double> x(ncols);
    if ((st == BPLEX_OPTIMAL || st == BPLEX_TIME_LIMIT || st == BPLEX_NODE_LIMIT) &&
        bplex_get_solution(bp_, x.data(), nullptr, nullptr, nullptr) == 0) {
        obj = bplex_objval(bp_);
        for (size_t j = 0; j < lp_route_.size(); j++)
            if (lp_route_[j] >= 0 && x[nfix_ + j] > 0.5) routes_used.push_back(lp_route_[j]);
    }
    bplex_set_param(bp_, "relax", 1);
    bplex_set_param(bp_, "cutoff", BPLEX_INF);
    bplex_set_param(bp_, "time_limit", lp_time_limit_);
    for (int j = 0; j < nfix_; j++) bplex_chg_col_bounds(bp_, j, 0.0, BPLEX_INF);
    for (size_t j = 0; j < lp_route_.size(); j++) {
        if (lp_route_[j] < 0) bplex_chg_col_bounds(bp_, nfix_ + (int)j, 0.0, BPLEX_INF);
        else bplex_chg_col_type(bp_, nfix_ + (int)j, 0);
    }
    return st;
}

double Master::probe(int id, int iter_limit)
{
    using namespace std::chrono;
    auto tnow = []() { return duration<double>(steady_clock::now().time_since_epoch()).count(); };
    set_row_active(id, true);
    double old_limit = bplex_get_param(bp_, "iter_limit");
    bplex_set_param(bp_, "iter_limit", iter_limit);
    bplex_set_param(bp_, "method", BPLEX_METHOD_DUAL);
    double t = tnow();
    int st = bplex_solve(bp_);
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
    bplex_set_param(bp_, "iter_limit", old_limit);
    bplex_set_param(bp_, "method", BPLEX_METHOD_AUTO);
    set_row_active(id, false);
    return val;
}

void Master::write_mps(const std::string &path) const { bplex_write_mps(bp_, path.c_str()); }

} // namespace cvrpbcp
