// Set-partitioning master on bplex (persistent engine, incremental edits).
//
// Rows:   0..n-1   customer i -> row i-1, sum_r a_ir lambda_r = 1
//         n        fleet: klo <= sum_r lambda_r <= khi (node-specific bounds)
//         n+1+p    robust row at position p: sum_r coef(r) lambda_r in [lo, up]
//                  RCC   coef = crossings of S         (lo = 2k, up = inf)
//                  EDGE  coef = traversals of edge {u,v} (branching: >= 1 or <= 0)
//                  an inactive row has free bounds (non-binding, slack basic)
// Columns: 0..n-1  artificials of the customer rows, n: artificial of the
//          fleet lower bound; then dynamic columns: LP routes and one
//          artificial per ">=" robust row. Every artificial costs bigM, so
//          the LP over the current columns is always feasible and, once
//          pricing converges with an artificial still positive, the node is
//          infeasible (bigM exceeds the cost of any route solution).
//
// Routes live in a pool (index r); a subset is in the LP (dynamic column
// j, lp_route_[j] = r, or -1 for a row artificial). Columns with positive
// reduced cost for several solves are purged from the LP (bplex_del_cols on
// nonbasic columns keeps the factorization) and come back through pool
// pricing when their reduced cost turns negative again.
//
// Robust rows have stable ids (what the tree stores); positions shift when
// rows are deleted. Rows are structural: every LP column carries its
// coefficient in every row, active or not.
#pragma once
#include "core/types.h"
#include "instance/instance.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

struct bplex_s;

namespace cvrpbcp {

// A master row beyond the partitioning and fleet rows. RCC and EDGE rows are
// robust (their duals enter arc costs); SRC rows are non-robust (limited-
// memory subset-row cuts: a state per cut in the pricing labels).
struct RobustRow {
    enum Kind { RCC, EDGE, SRC } kind = RCC;
    std::vector<char> in_s; // RCC: set S (size n+1, in_s[0] = 0); SRC: the triple C
    std::vector<char> in_m; // SRC: memory M (contains C)
    int k = 0;              // RCC: rhs = 2k
    int u = 0, v = 0;       // EDGE: endpoints, u < v
    double lo = 0, up = 0;  // bounds when active
    bool active = false;
    int id = -1;            // stable id
    int row = -1;           // current bplex row index
    int zero_streak = 0;    // consecutive LP solves with a zero dual (cut aging)
    int coef(const Route &r) const;
    int arc_coef(int i, int j) const; // 0 for SRC
};

class Pricer;

class Master {
public:
    Master(const Instance &I, int kmin, int kmax, const Params &P);
    ~Master();
    Master(const Master &) = delete;
    Master &operator=(const Master &) = delete;

    int n() const { return I_.n; }
    int klo() const { return klo_; }
    int khi() const { return khi_; }
    int nroutes() const { return (int)routes_.size(); } // pool size
    int nlpcols() const { return nlp_routes_; }          // routes in the LP
    const Route &route(int r) const { return routes_[r]; }
    const std::vector<RobustRow> &rows() const { return rows_; } // position order
    const RobustRow &row_by_id(int id) const { return rows_[pos_[id]]; }
    bool row_exists(int id) const { return id >= 0 && id < (int)pos_.size() && pos_[id] >= 0; }
    int row_id_bound() const { return (int)pos_.size(); } // every id ever issued is < this

    // returns the number of columns that entered the LP (duplicates of
    // routes already in the LP are skipped)
    int add_routes(const std::vector<Route> &rs);
    int price_pool(const Pricer &pr, double eps, int max_cols);
    int purge(int age, int min_cols);

    std::vector<int> add_rows(const std::vector<RobustRow> &rs);
    int add_row(const RobustRow &r) { return add_rows({r})[0]; }
    void set_row_active(int id, bool active);
    void delete_rows(const std::vector<int> &ids);
    void set_fleet_bounds(int klo, int khi);
    // diving: fix a pool route's column to 1 (it must be in the LP) / release it
    void fix_route(int r, bool fixed);

    int solve();
    // wall-clock cap for the next LP solves (bplex time_limit); a stopped
    // solve returns its status and the node treats it as an error
    void set_lp_time_limit(double seconds) { lp_time_limit_ = seconds; }
    double objval() const { return obj_; }
    const std::vector<double> &pi() const { return pi_; }          // pi[0] = fleet dual
    const std::vector<double> &row_duals() const { return yrow_; } // position order
    const std::vector<double> &lambda() const { return lambda_; }  // per pool route (0 if not in LP)
    double artificial_mass() const { return art_mass_; }

    double lagrangian_bound(const std::vector<double> &pi, const std::vector<double> &yrow, double min_rc) const;
    void edge_flows(std::vector<double> &x) const;

    int solve_restricted_mip(double cutoff, double time_limit, std::vector<int> &routes_used, double &obj);
    // Strong-branching probe on an existing inactive row: activate, run at
    // most iter_limit dual simplex iterations, deactivate. Returns a bound on
    // the child's LP over the current columns, or -1 when the probe gives no
    // information (artificials needed, i.e. columns are missing).
    double probe(int id, int iter_limit);

    int last_iterations() const { return last_iters_; }
    struct ProbeStats { int n = 0; double t_solve = 0, iters = 0; } probe_stats;
    long purged_total = 0, repriced_total = 0;
    bplex_s *raw() { return bp_; }
    void write_mps(const std::string &path) const;

private:
    const Instance &I_;
    const Params &P_;
    int klo_, khi_;
    int nfix_;                     // fixed artificial columns (n + 1)
    double bigM_;
    bplex_s *bp_ = nullptr;
    std::vector<Route> routes_;    // pool
    std::unordered_multimap<uint64_t, int> route_index_; // sequence hash -> pool index
    std::vector<int> lp_route_;    // dynamic column j -> pool route, -1 for a row artificial
    std::vector<int> lp_art_row_;  // dynamic column j -> row id of its artificial, -1 for routes
    std::vector<int> col_of_;      // pool route -> dynamic column j, -1 if out
    std::vector<int> age_;         // per dynamic column: consecutive solves with rc > 0
    int nlp_routes_ = 0;
    std::vector<RobustRow> rows_;  // position order
    std::vector<int> pos_;         // id -> position, -1 if deleted
    int nrows_ = 0;
    double obj_ = 0, art_mass_ = 0;
    std::vector<double> pi_, yrow_, lambda_, redcost_;
    int last_iters_ = 0;
    double lp_time_limit_ = 1e30;
    void append_lp_columns(const std::vector<int> &route_idx);
    [[noreturn]] void row_edit_failed(const char *what, int rc, int row) const;
    void delete_dynamic_columns(const std::vector<int> &js);
};

} // namespace cvrpbcp
