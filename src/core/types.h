// Shared small types and the parameter block.
#pragma once
#include <string>
#include <vector>

namespace cvrpbcp {

// A route: customer sequence (depot implicit at both ends) and its cost in
// cost units. ng-routes may visit a customer more than once.
struct Route {
    std::vector<int> seq;
    long long cost = 0;
    bool elementary() const;
};

struct Params {
    // pricing
    int ng_size = 8;              // ng-neighbourhood size (<= 16)
    bool bidirectional = true;    // labeling to Q/2 from both ends, then concatenation
    int heur_knn = 10;            // successors per vertex in heuristic pricing
    int heur_label_cap = 40;      // alive labels per vertex in heuristic pricing
    int heur_min_cols = 1;        // below this many columns from the sparse pass, try the full graph (capped)
    bool savings_init = true;     // seed the master (and the incumbent) with Clarke-Wright routes
    int max_cols_per_iter = 150;  // columns added per CG iteration
    double rc_eps = 1e-6;         // a route prices out when rc < -rc_eps
    int qbound_stage = 0;         // q-route completion bound prunes labels from this pricing stage on (3 = off)
    // master
    double smoothing_alpha = 0.0; // Wentges dual smoothing (0 = off); with smoothing_auto this is the start value
    bool smoothing_auto = true;   // adjust alpha by the subgradient angle rule (Pessoa et al. 2018)
    bool early_stop = true;       // stop CG once ceil(LB) == ceil(z)
    int purge_age = 6;            // LP columns nonbasic with rc > 0 for this many solves are purged (node start)
    int purge_min_cols_per_cust = 20; // never purge below this many LP columns per customer
    // cuts
    bool rcc = true;
    int rcc_max_per_round = 50;
    double rcc_min_violation = 0.05;
    int rcc_max_rounds = 100;         // total cut rounds per node (RCC and SRC)
    bool src = true;                  // limited-memory subset-row cuts on customer triples
    int src_max_per_round = 30;
    int src_max_active = 100;         // active SRC rows (each adds a state to the labels)
    double src_min_violation = 0.1;
    double cut_tailoff = 0.0005;      // stop SRC rounds when z improves less than this fraction twice
    int cut_age_limit = 30;           // delete a cut row after this many LP solves with a zero dual
    int cleanup_every = 25;           // nodes between row cleanups
    // limits
    double time_limit = 3600;
    int log_level = 1;
    int threads = 1;
};

} // namespace cvrpbcp
