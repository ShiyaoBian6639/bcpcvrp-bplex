#include "pricing/labeling.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace cvrpbcp {

Pricer::Pricer(const Instance &I, const Params &P) : I_(I), P_(P), nv_(I.nv())
{
    const int n = I.n;
    ngsz_ = std::min(std::max(P.ng_size, 1), 16);
    ngsz_ = std::min(ngsz_, n);
    ng_.assign(nv_, {});
    ngpos_.assign((size_t)nv_ * nv_, -1);
    knn_.assign(nv_, {});
    std::vector<int> order(n);
    for (int v = 1; v <= n; v++) {
        std::iota(order.begin(), order.end(), 1);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            if (I.d(v, a) != I.d(v, b)) return I.d(v, a) < I.d(v, b);
            return a < b;
        });
        std::vector<int> s(order.begin(), order.begin() + ngsz_);
        if (std::find(s.begin(), s.end(), v) == s.end()) s.back() = v;
        std::sort(s.begin(), s.end());
        ng_[v] = s;
        for (int a = 0; a < (int)s.size(); a++) ngpos_[(size_t)v * nv_ + s[a]] = (int8_t)a;
        int k = std::min(P.heur_knn, n - 1);
        for (int a = 0, c = 0; a < n && c < k; a++)
            if (order[a] != v) { knn_[v].push_back(order[a]); c++; }
    }
    transfer_.assign((size_t)nv_ * nv_ * 16, -1);
    for (int v = 1; v <= n; v++)
        for (int w = 1; w <= n; w++)
            for (int a = 0; a < (int)ng_[v].size(); a++)
                transfer_[((size_t)v * nv_ + w) * 16 + a] = ngpos_[(size_t)w * nv_ + ng_[v][a]];
    rc_.assign((size_t)nv_ * nv_, 0.0);
    forbid_.assign((size_t)nv_ * nv_, 0);
    removed_.assign(nv_, 0);
    inC_.assign(nv_, {});
    memmask_.assign((size_t)nv_ * 2, 0);
    at_.assign(nv_, {});
}

void Pricer::set_duals(const std::vector<double> &pi, const std::vector<RobustRow> &rows,
                       const std::vector<double> &yrow)
{
    qb_dirty_ = true;
    for (int i = 0; i < nv_; i++)
        for (int j = 0; j < nv_; j++) rc_[(size_t)i * nv_ + j] = I_.d(i, j) - 0.5 * (pi[i] + pi[j]);
    nsrc_ = 0;
    sigma_.clear();
    for (auto &l : inC_) l.clear();
    std::fill(memmask_.begin(), memmask_.end(), 0);
    for (size_t c = 0; c < rows.size(); c++) {
        double y = yrow[c];
        const RobustRow &r = rows[c];
        if (std::fabs(y) <= 1e-12 || !r.active) continue;
        if (r.kind == RobustRow::EDGE) {
            rc_[(size_t)r.u * nv_ + r.v] -= y;
            rc_[(size_t)r.v * nv_ + r.u] -= y;
        } else if (r.kind == RobustRow::RCC) {
            for (int i = 0; i < nv_; i++)
                for (int j = 0; j < nv_; j++)
                    if (r.in_s[i] != r.in_s[j]) rc_[(size_t)i * nv_ + j] -= y;
        } else { // SRC
            if (y >= -1e-9) continue;
            if (nsrc_ >= MAX_SRC) throw std::runtime_error("more active SRC rows than the pricer supports");
            int k = nsrc_++;
            sigma_.push_back(y);
            for (int v = 1; v < nv_; v++) {
                if (r.in_s[v]) inC_[v].push_back(k);
                if (r.in_m[v]) memmask_[(size_t)v * 2 + (k >> 6)] |= (uint64_t)1 << (k & 63);
            }
        }
    }
}

void Pricer::set_forbidden(const std::vector<std::pair<int, int>> &edges)
{
    qb_dirty_ = true;
    std::fill(forbid_.begin(), forbid_.end(), 0);
    for (auto [u, v] : edges) {
        forbid_[(size_t)u * nv_ + v] = 1;
        forbid_[(size_t)v * nv_ + u] = 1;
    }
}

// q-route DP over the load: for every load q and vertex v the cheapest and
// the second cheapest (different predecessor) walk depot -> ... -> v with
// load exactly q and no 2-cycle; then prefix minima over the load.
void Pricer::compute_qbound()
{
    const int n = I_.n, Q = I_.capacity;
    const double INF = 1e300;
    const size_t stride = (size_t)Q + 1;
    qb_.assign((size_t)nv_ * stride, INF);
    std::vector<double> b1((size_t)(Q + 1) * nv_, INF), b2((size_t)(Q + 1) * nv_, INF);
    std::vector<int> p1((size_t)(Q + 1) * nv_, -1);
    std::vector<int> all;
    for (int v = 1; v <= n; v++)
        if (!removed_[v] && I_.demand[v] <= Q) all.push_back(v);
    for (int v : all)
        if (!forbidden(0, v)) {
            b1[(size_t)I_.demand[v] * nv_ + v] = arc_rc(0, v);
            p1[(size_t)I_.demand[v] * nv_ + v] = 0;
        }
    for (int q = 1; q <= Q; q++) {
        const size_t row = (size_t)q * nv_;
        for (int v : all) {
            const double c1 = b1[row + v];
            if (c1 >= INF) continue;
            const double c2 = b2[row + v];
            const int pv = p1[row + v];
            const char *fb = &forbid_[(size_t)v * nv_];
            const double *rcv = &rc_[(size_t)v * nv_];
            for (int w : all) {
                if (w == v || fb[w]) continue;
                const int nq = q + I_.demand[w];
                if (nq > Q) continue;
                const double c = (pv == w) ? c2 : c1;
                if (c >= INF) continue;
                const double val = c + rcv[w];
                const size_t t = (size_t)nq * nv_ + w;
                if (val < b1[t]) {
                    if (p1[t] != v) b2[t] = b1[t];
                    b1[t] = val;
                    p1[t] = v;
                } else if (val < b2[t] && p1[t] != v) {
                    b2[t] = val;
                }
            }
        }
    }
    for (int v : all) {
        double m = INF;
        for (int q = 0; q <= Q; q++) {
            m = std::min(m, b1[(size_t)q * nv_ + v]);
            qb_[(size_t)v * stride + q] = m;
        }
    }
    qb_dirty_ = false;
}

double Pricer::route_rc(const Route &r) const
{
    double rc = 0;
    int prev = 0;
    Label L{0, 0, 0, -1, 0, false, {0, 0}};
    for (int v : r.seq) {
        rc += arc_rc(prev, v);
        if (nsrc_) { L.rc = 0; src_enter(L, v); rc += L.rc; }
        prev = v;
    }
    return rc + arc_rc(prev, 0);
}

double Pricer::penalty(uint64_t d0, uint64_t d1) const
{
    double p = 0;
    while (d0) { int k = __builtin_ctzll(d0); p -= sigma_[k]; d0 &= d0 - 1; }
    while (d1) { int k = __builtin_ctzll(d1); p -= sigma_[64 + k]; d1 &= d1 - 1; }
    return p;
}

void Pricer::src_enter(Label &N, int w) const
{
    N.s[0] &= memmask_[(size_t)w * 2];
    N.s[1] &= memmask_[(size_t)w * 2 + 1];
    for (int k : inC_[w]) {
        uint64_t &word = N.s[k >> 6];
        uint64_t bit = (uint64_t)1 << (k & 63);
        if (word & bit) { word &= ~bit; N.rc -= sigma_[k]; } // coefficient +1: rc grows by |sigma|
        else word |= bit;
    }
}

bool Pricer::insert(Label L)
{
    auto &list = at_[L.v];
    for (int idx : list) {
        const Label &o = arena_[idx];
        if (o.dead) continue;
        if (o.load <= L.load && o.time <= L.time && (o.mem & ~L.mem) == 0) {
            double orc = o.rc;
            if (nsrc_) orc += penalty(o.s[0] & ~L.s[0], o.s[1] & ~L.s[1]);
            if (orc <= L.rc + 1e-9) return false;
        }
    }
    size_t w = 0;
    for (size_t k = 0; k < list.size(); k++) {
        Label &o = arena_[list[k]];
        if (o.dead) continue;
        if (L.load <= o.load && L.time <= o.time && (L.mem & ~o.mem) == 0) {
            double lrc = L.rc;
            if (nsrc_) lrc += penalty(L.s[0] & ~o.s[0], L.s[1] & ~o.s[1]);
            if (lrc <= o.rc) { o.dead = true; continue; }
        }
        list[w++] = list[k];
    }
    list.resize(w);
    if (cap_ > 0 && (int)list.size() >= cap_) {
        // heuristic pricing: bounded label count per vertex, keep the best by rc
        size_t worst = 0;
        for (size_t k = 1; k < list.size(); k++)
            if (arena_[list[k]].rc > arena_[list[worst]].rc) worst = k;
        if (arena_[list[worst]].rc <= L.rc) return false;
        arena_[list[worst]].dead = true;
        list[worst] = list.back();
        list.pop_back();
    }
    arena_.push_back(L);
    list.push_back((int)arena_.size() - 1);
    return true;
}

Route Pricer::build_route(int label) const
{
    Route r;
    for (int l = label; l >= 0; l = arena_[l].parent) r.seq.push_back(arena_[l].v);
    std::reverse(r.seq.begin(), r.seq.end());
    int prev = 0;
    for (int v : r.seq) { r.cost += I_.d(prev, v); prev = v; }
    r.cost += I_.d(prev, 0);
    return r;
}

void Pricer::label_pass(bool heuristic, int extend_limit, int cap)
{
    const int n = I_.n, Q = I_.capacity;
    const bool prune = prune_thr_ < 1e299 && !qb_dirty_;
    cap_ = cap;
    arena_.clear();
    for (auto &l : at_) l.clear();
    std::vector<std::vector<int>> by_load(Q + 1);
    for (int j = 1; j <= n; j++) {
        if (I_.demand[j] > Q || forbidden(0, j) || removed_[j]) continue;
        if (I_.has_tw && !ignore_tw_ && I_.tw_lo[0] + I_.d(0, j) > I_.tw_hi[j]) continue;
        if (prune && pruned(j, I_.demand[j], arc_rc(0, j))) continue;
        Label L{arc_rc(0, j), I_.demand[j], j, -1, (uint16_t)(1u << ngpos_[(size_t)j * nv_ + j]), false, {0, 0},
                start_time(j)};
        if (nsrc_) src_enter(L, j);
        if (insert(L)) by_load[L.load].push_back((int)arena_.size() - 1);
    }
    std::vector<int> all;
    for (int v = 1; v <= n; v++) if (!removed_[v]) all.push_back(v);
    aborted_ = false;
    size_t processed = 0;
    for (int q = 1; q <= extend_limit; q++) {
        for (size_t k = 0; k < by_load[q].size(); k++) {
            const int idx = by_load[q][k];
            if (arena_[idx].dead) continue;
            if (deadline_ > 0 && (++processed & 255) == 0 &&
                std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count() > deadline_) {
                aborted_ = true;
                labels_last_ = arena_.size();
                return;
            }
            const Label L = arena_[idx]; // copy: arena_ may reallocate
            const std::vector<int> &succ = heuristic ? knn_[L.v] : all;
            const char *fb = &forbid_[(size_t)L.v * nv_];
            const double *rcv = &rc_[(size_t)L.v * nv_];
            for (int w : succ) {
                if (w == L.v || fb[w] || removed_[w]) continue;
                // arc-level bound check before the label is built (SRC
                // penalties only add, so the arc reduced cost suffices)
                if (prune && L.load + I_.demand[w] <= Q && pruned(w, L.load + I_.demand[w], L.rc + rcv[w])) continue;
                Label N;
                if (!extend(L, idx, w, N)) continue;
                if (insert(N)) by_load[N.load].push_back((int)arena_.size() - 1);
            }
        }
    }
    labels_last_ = arena_.size();
}

bool Pricer::extend(const Label &L, int idx, int w, Label &N) const
{
    int8_t pw = ngpos_[(size_t)L.v * nv_ + w];
    if (pw >= 0 && (L.mem >> pw & 1)) return false; // w in memory: forbidden
    int load = L.load + I_.demand[w];
    if (load > I_.capacity) return false;
    int time = 0;
    if (I_.has_tw && !ignore_tw_) {
        int start = std::max(L.time + I_.d(L.v, w), I_.tw_lo[w]);
        if (start > I_.tw_hi[w]) return false;
        time = start + I_.service[w];
    }
    const int8_t *t = &transfer_[((size_t)L.v * nv_ + w) * 16];
    uint16_t mem = 0;
    for (int a = 0; a < ngsz_; a++)
        if ((L.mem >> a & 1) && t[a] >= 0) mem |= (uint16_t)(1u << t[a]);
    mem |= (uint16_t)(1u << ngpos_[(size_t)w * nv_ + w]);
    N = Label{L.rc + arc_rc(L.v, w), load, w, idx, mem, false, {L.s[0], L.s[1]}, time};
    if (nsrc_) src_enter(N, w);
    return true;
}

int Pricer::start_time(int j) const
{
    if (!I_.has_tw || ignore_tw_) return 0;
    return std::max(I_.tw_lo[0] + I_.d(0, j), I_.tw_lo[j]) + I_.service[j];
}

std::vector<Route> Pricer::price(int stage, int max_routes, double &min_rc)
{
    // completion-bound pruning from stage qbound_stage on; the bound is
    // computed once per set of duals (O(Q n^2))
    if (stage >= P_.qbound_stage && qb_dirty_) compute_qbound();
    prune_thr_ = (!qb_dirty_ && stage >= P_.qbound_stage) ? -1e-9 : 1e300;
    prune_vpay_ = nullptr;
    std::vector<Route> out;
    if (stage == 0) { stage_cap_ = P_.heur_label_cap; out = price_mono(true, max_routes, min_rc); }
    else if (stage == 1) { stage_cap_ = P_.heur_label_cap; out = price_mono(false, max_routes, min_rc); }
    else {
        stage_cap_ = 0;
        // the bidirectional join relies on symmetric half-paths: not with time windows
        out = (P_.bidirectional && !I_.has_tw) ? price_bidir(false, max_routes, min_rc)
                                                : price_mono(false, max_routes, min_rc);
    }
    prune_thr_ = 1e300;
    return out;
}

std::vector<Route> Pricer::price_mono(bool heuristic, int max_routes, double &min_rc)
{
    label_pass(heuristic, I_.capacity, stage_cap_);
    struct Cand { double rc; int label; };
    std::vector<Cand> cands;
    min_rc = 0;
    const double eps = P_.rc_eps;
    for (int v = 1; v <= I_.n; v++) {
        if (forbidden(v, 0)) continue;
        for (int idx : at_[v]) {
            const Label &L = arena_[idx];
            if (L.dead || !can_return(v, L.time)) continue;
            double total = L.rc + arc_rc(v, 0);
            if (total < min_rc) min_rc = total;
            if (total < -eps) cands.push_back({total, idx});
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand &a, const Cand &b) { return a.rc < b.rc; });
    std::vector<Route> out;
    for (const Cand &c : cands) {
        if ((int)out.size() >= max_routes) break;
        out.push_back(build_route(c.label));
    }
    return out;
}

// Bidirectional labeling on the load resource. The graph, the reduced costs
// and the row coefficients are symmetric, so one labeling pass gives both
// the forward half-paths and (read backwards) the backward half-paths.
// Labels are extended only while load <= Q/2; labels created beyond that are
// kept as possible backward ends. A route with load > Q/2 is generated once,
// at the unique split where the forward prefix has load <= Q/2 and adding
// the next customer exceeds Q/2; routes with load <= Q/2 are complete
// forward paths. Joining f (at v) and g (at w) is ng-feasible iff the two
// ng-memories are disjoint as customer sets (Baldacci et al. 2011); an SRC
// state set on both sides completes one more pair at the junction.
std::vector<Route> Pricer::price_bidir(bool heuristic, int max_routes, double &min_rc)
{
    const int n = I_.n, Q = I_.capacity, half = Q / 2;
    label_pass(heuristic, half);

    // per vertex: alive labels sorted by load, with prefix minima of rc
    std::vector<std::vector<int>> lab(nv_);
    std::vector<std::vector<double>> pmin(nv_);
    for (int v = 1; v <= n; v++) {
        for (int idx : at_[v])
            if (!arena_[idx].dead) lab[v].push_back(idx);
        std::sort(lab[v].begin(), lab[v].end(), [&](int a, int b) { return arena_[a].load < arena_[b].load; });
        pmin[v].resize(lab[v].size());
        double m = 1e300;
        for (size_t i = 0; i < lab[v].size(); i++) { m = std::min(m, arena_[lab[v][i]].rc); pmin[v][i] = m; }
    }

    struct Cand { double rc; int f, g; }; // g = -1: complete forward route
    std::vector<Cand> cands;
    min_rc = 0;
    const double eps = P_.rc_eps;
    std::vector<int> fmem; // forward memory as customers
    for (int v = 1; v <= n; v++) {
        for (int fi : lab[v]) {
            const Label &F = arena_[fi];
            if (F.load > half) continue;
            if (!forbidden(v, 0)) {
                double total = F.rc + arc_rc(v, 0);
                if (total < min_rc) min_rc = total;
                if (total < -eps) cands.push_back({total, fi, -1});
            }
            fmem.clear();
            for (int a = 0; a < ngsz_; a++)
                if (F.mem >> a & 1) fmem.push_back(ng_[v][a]);
            const char *fb = &forbid_[(size_t)v * nv_];
            for (int w = 1; w <= n; w++) {
                if (w == v || fb[w] || lab[w].empty()) continue;
                if (F.load + I_.demand[w] <= half) continue; // split rule
                int8_t pw = ngpos_[(size_t)v * nv_ + w];
                if (pw >= 0 && (F.mem >> pw & 1)) continue;   // w in forward memory
                int maxload = Q - F.load;
                const auto &lw = lab[w];
                int hi = (int)(std::upper_bound(lw.begin(), lw.end(), maxload,
                                                [&](int ld, int idx) { return ld < arena_[idx].load; }) -
                               lw.begin());
                if (hi == 0) continue;
                double base = F.rc + arc_rc(v, w);
                if (base + pmin[w][hi - 1] >= -eps && base + pmin[w][hi - 1] >= min_rc) continue;
                for (int gi = 0; gi < hi; gi++) {
                    const Label &G = arena_[lw[gi]];
                    double total = base + G.rc;
                    if (total >= -eps && total >= min_rc) continue;
                    bool ok = true;
                    for (int u : fmem) {
                        int8_t p = ngpos_[(size_t)w * nv_ + u];
                        if (p >= 0 && (G.mem >> p & 1)) { ok = false; break; }
                    }
                    if (!ok) continue;
                    if (nsrc_) total += penalty(F.s[0] & G.s[0], F.s[1] & G.s[1]);
                    if (total < min_rc) min_rc = total;
                    if (total < -eps) cands.push_back({total, fi, lw[gi]});
                }
            }
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand &a, const Cand &b) { return a.rc < b.rc; });
    std::vector<Route> out;
    for (const Cand &c : cands) {
        if ((int)out.size() >= max_routes) break;
        Route r = build_route(c.f);
        if (c.g >= 0) {
            for (int l = c.g; l >= 0; l = arena_[l].parent) r.seq.push_back(arena_[l].v);
            r.cost = 0;
            int prev = 0;
            for (int v : r.seq) { r.cost += I_.d(prev, v); prev = v; }
            r.cost += I_.d(prev, 0);
        }
        out.push_back(std::move(r));
    }
    return out;
}

// Enumeration of elementary routes with reduced cost <= gap. Completion
// bound from the mono ng pass: by symmetry, the cheapest completion from
// (v, load q) costs at least min{G.rc : G an ng label at v, G.load <= Q - q +
// d_v} (G read backwards is a path v -> depot; SRC penalties only add).
bool Pricer::enumerate(double gap, size_t max_labels, size_t max_routes, std::vector<Route> &out, double time_cap)
{
    using namespace std::chrono;
    const auto t_start = steady_clock::now();
    auto over_time = [&]() { return duration<double>(steady_clock::now() - t_start).count() > time_cap; };
    out.clear();
    deadline_ = time_cap < 1e299 ? duration<double>(t_start.time_since_epoch()).count() + time_cap : 0;
    const int n = I_.n, Q = I_.capacity;
    // The backward label at v includes v's own visit: an SRC pair (x, v)
    // paid by that label may already be paid by the prefix (v paired with a
    // prefix vertex), so subtract what v itself can pay. Elsewhere the
    // junction only adds penalties.
    std::vector<double> vpay(nv_, 0.0);
    for (int v = 1; v <= n; v++)
        for (int k : inC_[v]) vpay[v] -= sigma_[k];
    // ng pass pruned by the q-route bound: a label G at v with
    // G.rc + qcompletion(G) > gap + vpay[v] cannot be the reversed suffix of
    // a route with rc <= gap (route rc >= prefix rc + G.rc - vpay[v], and the
    // reversed prefix is a q-walk, so prefix rc >= qcompletion(G))
    if (P_.qbound_stage <= 2 && qb_dirty_) compute_qbound();
    prune_thr_ = qb_dirty_ ? 1e300 : gap + 1e-6;
    prune_vpay_ = &vpay;
    ignore_tw_ = true; // the bound must cover completions in either direction: time-free labels
    label_pass(false, Q);
    ignore_tw_ = false;
    prune_thr_ = 1e300;
    prune_vpay_ = nullptr;
    deadline_ = 0;
    if (aborted_ || over_time()) return false;
    // bound[v][L] = min rc over ng labels at v with load <= L
    std::vector<std::vector<double>> bound(nv_);
    for (int v = 1; v <= n; v++) {
        bound[v].assign(Q + 1, 1e300);
        for (int idx : at_[v]) {
            const Label &L = arena_[idx];
            if (!L.dead && L.rc < bound[v][L.load]) bound[v][L.load] = L.rc;
        }
        for (int L = 1; L <= Q; L++) bound[v][L] = std::min(bound[v][L], bound[v][L - 1]);
    }
    auto completion = [&](int v, int q) {
        int L = Q - q + I_.demand[v];
        if (L > Q) L = Q;
        return bound[v][L] - vpay[v];
    };

    // elementary labels: visited set in words
    const int W = (n + 64) / 64;
    // cost: true distance so far. With arc-dependent duals (RCC, EDGE rows)
    // the reduced cost alone cannot decide dominance between two orders of
    // the same customer set: the cheaper-rc order may be the dearer route
    // and the MIP needs the cheap one. Both must be no worse.
    struct EL { double rc; long long cost; int load, v, parent; uint64_t s[2]; bool dead; int time; };
    std::vector<EL> el;
    std::vector<uint64_t> vis; // W words per label
    std::unordered_map<uint64_t, std::vector<int>> buckets; // (v, visited set) -> labels
    std::vector<std::vector<int>> by_load(Q + 1);
    const double tol = 1e-6;
    auto key_of = [&](int v, const uint64_t *set) {
        uint64_t h = (uint64_t)v * 0x9E3779B97F4A7C15ull;
        for (int w = 0; w < W; w++) h = (h ^ set[w]) * 0xBF58476D1CE4E5B9ull + (uint64_t)w;
        return h;
    };
    auto push = [&](EL L, const uint64_t *set) -> bool {
        // dominance among labels with the same visited set at v
        auto &b = buckets[key_of(L.v, set)];
        for (int idx : b) {
            EL &o = el[idx];
            if (o.dead || o.v != L.v) continue;
            if (!std::equal(set, set + W, &vis[(size_t)idx * W])) continue;
            double orc = o.rc + (nsrc_ ? penalty(o.s[0] & ~L.s[0], o.s[1] & ~L.s[1]) : 0.0);
            if (o.cost <= L.cost && orc <= L.rc + 1e-9) return false;
            double lrc = L.rc + (nsrc_ ? penalty(L.s[0] & ~o.s[0], L.s[1] & ~o.s[1]) : 0.0);
            if (L.cost <= o.cost && lrc <= o.rc) o.dead = true;
        }
        el.push_back(L);
        vis.insert(vis.end(), set, set + W);
        b.push_back((int)el.size() - 1);
        by_load[L.load].push_back((int)el.size() - 1);
        return true;
    };
    std::vector<uint64_t> tmp(W);
    Label scratch{0, 0, 0, -1, 0, false, {0, 0}};
    for (int j = 1; j <= n; j++) {
        if (I_.demand[j] > Q || forbidden(0, j)) continue;
        if (I_.has_tw && I_.tw_lo[0] + I_.d(0, j) > I_.tw_hi[j]) continue;
        EL L{arc_rc(0, j), (long long)I_.d(0, j), I_.demand[j], j, -1, {0, 0}, false, start_time(j)};
        if (nsrc_) { scratch.rc = 0; scratch.s[0] = scratch.s[1] = 0; src_enter(scratch, j); L.rc += scratch.rc; L.s[0] = scratch.s[0]; L.s[1] = scratch.s[1]; }
        if (L.rc + completion(j, L.load) > gap + tol) continue;
        std::fill(tmp.begin(), tmp.end(), 0);
        tmp[j >> 6] |= (uint64_t)1 << (j & 63);
        push(L, tmp.data());
    }
    std::vector<int> routes;
    for (int q = 1; q <= Q; q++) {
        for (size_t k = 0; k < by_load[q].size(); k++) {
            const int idx = by_load[q][k];
            if (el[idx].dead) continue;
            const EL L = el[idx];
            if (!forbidden(L.v, 0) && can_return(L.v, L.time) && L.rc + arc_rc(L.v, 0) <= gap + tol) {
                routes.push_back(idx);
                if (routes.size() > max_routes) return false;
            }
            const char *fb = &forbid_[(size_t)L.v * nv_];
            for (int w = 1; w <= n; w++) {
                if (w == L.v || fb[w]) continue;
                if (vis[(size_t)idx * W + (w >> 6)] >> (w & 63) & 1) continue;
                int load = L.load + I_.demand[w];
                if (load > Q) continue;
                int time = 0;
                if (I_.has_tw) {
                    int start = std::max(L.time + I_.d(L.v, w), I_.tw_lo[w]);
                    if (start > I_.tw_hi[w]) continue;
                    time = start + I_.service[w];
                }
                EL N{L.rc + arc_rc(L.v, w), L.cost + I_.d(L.v, w), load, w, idx, {L.s[0], L.s[1]}, false, time};
                if (nsrc_) { scratch.rc = 0; scratch.s[0] = N.s[0]; scratch.s[1] = N.s[1]; src_enter(scratch, w); N.rc += scratch.rc; N.s[0] = scratch.s[0]; N.s[1] = scratch.s[1]; }
                if (N.rc + completion(w, load) > gap + tol) continue;
                std::copy(&vis[(size_t)idx * W], &vis[(size_t)idx * W] + W, tmp.begin());
                tmp[w >> 6] |= (uint64_t)1 << (w & 63);
                push(N, tmp.data());
                if (el.size() > max_labels) return false;
                if ((el.size() & 8191) == 0 && over_time()) return false;
            }
        }
    }
    out.reserve(routes.size());
    for (int idx : routes) {
        Route r;
        for (int l = idx; l >= 0; l = el[l].parent) r.seq.push_back(el[l].v);
        std::reverse(r.seq.begin(), r.seq.end());
        int prev = 0;
        for (int v : r.seq) { r.cost += I_.d(prev, v); prev = v; }
        r.cost += I_.d(prev, 0);
        out.push_back(std::move(r));
    }
    return true;
}

} // namespace cvrpbcp

namespace cvrpbcp {

void Pricer::debug_trace_route(const Route &r, double gap)
{
    const int n = I_.n, Q = I_.capacity;
    label_pass(false, Q);
    std::vector<std::vector<double>> bound(nv_);
    for (int v = 1; v <= n; v++) {
        bound[v].assign(Q + 1, 1e300);
        for (int idx : at_[v]) {
            const Label &L = arena_[idx];
            if (!L.dead && L.rc < bound[v][L.load]) bound[v][L.load] = L.rc;
        }
        for (int L = 1; L <= Q; L++) bound[v][L] = std::min(bound[v][L], bound[v][L - 1]);
    }
    // prefix reduced costs (with SRC) and suffix costs (route order)
    std::vector<double> prefix(r.seq.size());
    Label L{0, 0, 0, -1, 0, false, {0, 0}};
    double rc = 0;
    int prev = 0, load = 0;
    std::printf("  TRACE route:");
    for (int v : r.seq) std::printf(" %d", v);
    std::printf("  total rc %.4f gap %.4f\n", route_rc(r), gap);
    for (size_t k = 0; k < r.seq.size(); k++) {
        int v = r.seq[k];
        rc += arc_rc(prev, v);
        if (nsrc_) { L.rc = 0; src_enter(L, v); rc += L.rc; }
        load += I_.demand[v];
        prefix[k] = rc;
        // true suffix: from v to the depot along the route (SRC states continue)
        double suf = 0;
        Label S = L;
        int p = v;
        for (size_t j = k + 1; j < r.seq.size(); j++) {
            suf += arc_rc(p, r.seq[j]);
            if (nsrc_) { S.rc = 0; src_enter(S, r.seq[j]); suf += S.rc; }
            p = r.seq[j];
        }
        suf += arc_rc(p, 0);
        int Lc = std::min(Q, Q - load + I_.demand[v]);
        double vpay = 0;
        for (int kk : inC_[v]) vpay -= sigma_[kk];
        double b = bound[v][Lc] - vpay;
        std::printf("  TRACE   prefix to %3d: load %4d rc %9.4f  bound %9.4f  suffix %9.4f  %s%s\n", v, load, rc, b, suf,
                    rc + b > gap + 1e-6 ? "PRUNED " : "", b > suf + 1e-6 ? "BOUND>SUFFIX" : "");
        if (b > suf + 1e-6) {
            // arc-by-arc symmetry check of the suffix, with fresh SRC states both ways
            {
                std::vector<int> fw(r.seq.begin() + k, r.seq.end());
                double f = 0, g = 0;
                Label A{0, 0, 0, -1, 0, false, {0, 0}}, B{0, 0, 0, -1, 0, false, {0, 0}};
                std::printf("  TRACE     forward standalone:");
                int p = fw[0];
                if (nsrc_) { A.rc = 0; src_enter(A, p); f += A.rc; std::printf(" [%d pay %.2f]", p, A.rc); }
                for (size_t j = 1; j < fw.size(); j++) {
                    double a = arc_rc(p, fw[j]);
                    f += a;
                    double pay = 0;
                    if (nsrc_) { A.rc = 0; src_enter(A, fw[j]); pay = A.rc; f += pay; }
                    std::printf(" %d->%d %.2f [pay %.2f]", p, fw[j], a, pay);
                    p = fw[j];
                }
                f += arc_rc(p, 0);
                std::printf(" %d->0 %.2f  = %.4f\n", p, arc_rc(p, 0), f);
                std::printf("  TRACE     reversed standalone:");
                p = 0;
                for (size_t j = fw.size(); j-- > 0;) {
                    double a = arc_rc(p, fw[j]);
                    g += a;
                    double pay = 0;
                    if (nsrc_) { B.rc = 0; src_enter(B, fw[j]); pay = B.rc; g += pay; }
                    std::printf(" %d->%d %.2f [pay %.2f]", p, fw[j], a, pay);
                    p = fw[j];
                }
                std::printf("  = %.4f\n", g);
            }
            // rebuild the reversed suffix depot -> ... -> v as forward labels and
            // report, at each step, the alive label that dominates it
            std::vector<int> rev(r.seq.begin() + k, r.seq.end());
            std::reverse(rev.begin(), rev.end());
            Label G{arc_rc(0, rev[0]), I_.demand[rev[0]], rev[0], -1,
                    (uint16_t)(1u << ngpos_[(size_t)rev[0] * nv_ + rev[0]]), false, {0, 0}};
            if (nsrc_) src_enter(G, rev[0]);
            for (size_t j = 0; j < rev.size(); j++) {
                if (j > 0) {
                    Label N;
                    if (!extend(G, -1, rev[j], N)) {
                        std::printf("  TRACE     reversed suffix: ng or capacity forbids %d after %d\n", rev[j], G.v);
                        break;
                    }
                    G = N;
                }
                int dom = -1;
                for (int idx : at_[G.v]) {
                    const Label &o = arena_[idx];
                    if (o.dead) continue;
                    if (o.load <= G.load && (o.mem & ~G.mem) == 0) {
                        double orc = o.rc + (nsrc_ ? penalty(o.s[0] & ~G.s[0], o.s[1] & ~G.s[1]) : 0.0);
                        if (orc <= G.rc + 1e-9) { dom = idx; break; }
                    }
                }
                int alive = 0;
                for (int idx : at_[G.v]) alive += !arena_[idx].dead;
                std::printf("  TRACE     reversed suffix at %3d: load %4d rc %9.4f mem %04x s %llx; alive here %d, min rc (load<=%d) %.4f; ",
                            G.v, G.load, G.rc, G.mem, (unsigned long long)G.s[0], alive, G.load, bound[G.v][G.load]);
                if (dom >= 0)
                    std::printf("dominated by rc %.4f load %d mem %04x s %llx\n", arena_[dom].rc, arena_[dom].load,
                                arena_[dom].mem, (unsigned long long)arena_[dom].s[0]);
                else std::printf("no dominator\n");
            }
        }
        prev = v;
    }
}

} // namespace cvrpbcp
