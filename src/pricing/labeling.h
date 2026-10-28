// Pricing by labeling with the ng-route relaxation (load resource; time
// added in milestone 5). Labels are processed in increasing load, which is a
// topological order because every demand is >= 1.
//
// Arc reduced costs use the symmetric split: rc(i,j) = d(i,j) - (pi_i +
// pi_j)/2 - sum_c y_c coef_c(i,j) over the robust rows, with pi_0 the fleet
// dual, so a route's reduced cost is c_r - sum pi - pi_0 - row terms,
// matching the master's columns. Forbidden edges (branching x_e = 0) are
// never extended.
//
// Limited-memory subset-row cuts (SRC rows with negative duals) add one
// state bit per cut to the labels: visiting a member of C toggles the bit,
// and the toggle from 1 to 0 charges -dual; leaving the cut's memory resets
// the bit. Dominance charges the dominating label for every bit it has set
// and the other has not (Jepsen et al. 2008).
//
// Two pricing modes: mono-directional (labels to the full capacity) and
// bidirectional (labels to Q/2, then concatenation). Route enumeration
// (all elementary routes with reduced cost <= gap) reuses the mono pass as
// a completion bound.
#pragma once
#include "core/types.h"
#include "instance/instance.h"
#include "master/master.h"

#include <cstdint>
#include <vector>

namespace cvrpbcp {

class Pricer {
public:
    static constexpr int MAX_SRC = 128;

    Pricer(const Instance &I, const Params &P);

    // throws if more than MAX_SRC active SRC rows carry a negative dual
    void set_duals(const std::vector<double> &pi, const std::vector<RobustRow> &rows,
                   const std::vector<double> &yrow);
    void set_forbidden(const std::vector<std::pair<int, int>> &edges);
    // customers that pricing must not visit (diving: covered by fixed columns)
    void set_removed(const std::vector<char> &removed) { removed_ = removed; qb_dirty_ = true; }
    void clear_removed() { std::fill(removed_.begin(), removed_.end(), 0); qb_dirty_ = true; }

    // Routes with reduced cost < -rc_eps, best first, at most max_routes.
    // stage 0: sparse graph (heur_knn successors) with a label cap;
    // stage 1: full graph with a label cap; stage 2: exact (min_rc is then
    // the minimum reduced cost over all ng-routes, 0 when none is negative).
    std::vector<Route> price(int stage, int max_routes, double &min_rc);

    // All elementary routes with reduced cost <= gap under the current duals
    // and forbidden edges. Returns false (out empty) when the label or route
    // caps are exceeded. Elementary dominance: same (vertex, visited set),
    // reduced cost plus SRC penalty.
    bool enumerate(double gap, size_t max_labels, size_t max_routes, std::vector<Route> &out,
                   double time_cap = 1e300);

    double arc_rc(int i, int j) const { return rc_[(size_t)i * nv_ + j]; }
    bool forbidden(int i, int j) const { return forbid_[(size_t)i * nv_ + j] != 0; }
    // CVRPTW: can a label at v with service ended at t return to the depot in time?
    bool can_return(int v, int t) const { return !I_.has_tw || ignore_tw_ || t + I_.d(v, 0) <= I_.tw_hi[0]; }
    size_t labels_last() const { return labels_last_; }
    int nsrc() const { return nsrc_; }
    // reduced cost of a given route under the current duals (incl. SRC)
    double route_rc(const Route &r) const;
    // debug: for each prefix of r print its reduced cost, the enumeration's
    // completion bound and the true suffix cost (must satisfy bound <= suffix)
    void debug_trace_route(const Route &r, double gap);

private:
    struct Label {
        double rc;
        int load;
        int v;
        int parent;
        uint16_t mem;
        bool dead;
        uint64_t s[2]; // SRC states
        int time = 0;  // CVRPTW: time when service at v ends (cost units)
    };
    bool ignore_tw_ = false; // enumeration bound pass: time-free labels
    double deadline_ = 0;    // label_pass aborts past this wall-clock time (0: none)
    bool aborted_ = false;

    // Completion bound from q-routes (walks without 2-cycles, load resource,
    // SRC ignored): qb_[v*(Q+1) + L] = min reduced cost of a q-walk
    // depot -> ... -> v with load <= L (arc into v included). By symmetry it
    // bounds the cheapest completion of a label at v. Labels and extensions
    // whose rc + bound exceeds prune_thr_ (+ prune_vpay_[v], see enumerate)
    // are dropped; valid because SRC penalties only add and the ng, time and
    // elementarity restrictions only remove completions.
    std::vector<double> qb_;
    bool qb_dirty_ = true;
    double prune_thr_ = 1e300;            // off
    const std::vector<double> *prune_vpay_ = nullptr;
    void compute_qbound();
    double qcompletion(int v, int load) const
    {
        int L = I_.capacity - load + I_.demand[v];
        if (L > I_.capacity) L = I_.capacity;
        return qb_[(size_t)v * (I_.capacity + 1) + L];
    }
    bool pruned(int v, int load, double rc) const
    {
        return rc + qcompletion(v, load) > prune_thr_ + (prune_vpay_ ? (*prune_vpay_)[v] : 0.0);
    }
    const Instance &I_;
    const Params &P_;
    int nv_, ngsz_;
    std::vector<std::vector<int>> ng_;      // NG(v) as a sorted customer list, v included
    std::vector<int8_t> ngpos_;             // ngpos_[v*nv+u] = index of u in NG(v) or -1
    std::vector<int8_t> transfer_;          // [(v*nv+w)*16 + a] = index in NG(w) of NG(v)[a], or -1
    std::vector<std::vector<int>> knn_;     // heuristic successor lists
    std::vector<double> rc_;
    std::vector<char> forbid_;
    std::vector<char> removed_;
    // SRC data for the current duals
    int nsrc_ = 0;
    std::vector<double> sigma_;             // dual (<= 0) per SRC index
    std::vector<std::vector<int>> inC_;     // per vertex: SRC indices whose triple contains it
    std::vector<uint64_t> memmask_;         // per vertex: 2 words, bit set if the vertex is in the cut's memory
    std::vector<Label> arena_;
    // alive label indices per (vertex, load bucket), each list sorted by rc:
    // at_[v*nb_ + load/bstep_]. A label can only be dominated by labels in
    // its own or a lower bucket with rc <= its own, and can only dominate
    // labels in its own or a higher bucket with rc >= its own.
    std::vector<std::vector<int>> at_;
    int nb_ = 1, bstep_ = 1;
    std::vector<int> nalive_; // alive labels per vertex (heuristic cap)
    template <class F> void each_label(int v, F f) const
    {
        for (int b = 0; b < nb_; b++)
            for (int idx : at_[(size_t)v * nb_ + b]) f(idx);
    }
    size_t labels_last_ = 0;

    // ng labeling from the depot; labels with load <= extend_limit are
    // extended, those created beyond it are kept (bidirectional ends).
    // cap > 0: at most cap alive labels per vertex (heuristic pricing)
    void label_pass(bool heuristic, int extend_limit, int cap = 0);
    int cap_ = 0, stage_cap_ = 0;
    bool relax_dom_ = false; // heuristic stages: dominance ignores ng memory and SRC states
    void src_enter(Label &N, int w) const;
    // ng extension of L (arena index idx, -1 if none) to w; false if forbidden by ng or capacity
    bool extend(const Label &L, int idx, int w, Label &N) const;
    // CVRPTW: service end time of a first visit to j (0 without time windows)
    int start_time(int j) const;
    bool insert(Label L);
    Route build_route(int label) const;
    double penalty(uint64_t d0, uint64_t d1) const;
    std::vector<Route> price_mono(bool heuristic, int max_routes, double &min_rc);
    std::vector<Route> price_bidir(bool heuristic, int max_routes, double &min_rc);
};

} // namespace cvrpbcp
