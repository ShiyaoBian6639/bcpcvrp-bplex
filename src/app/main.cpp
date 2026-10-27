// cvrpbcp command line.
//   cvrpbcp selftest [cases]          brute force vs bplex set partitioning vs BCP on random tiny instances
//   cvrpbcp info <inst.vrp>           print instance summary
//   cvrpbcp check <inst.vrp> <sol>    verify a .sol against the instance (and the COMMENT optimum)
//   cvrpbcp root <inst.vrp> [opts]    root column-and-cut generation bound
//   cvrpbcp solve <inst.vrp> [opts]   branch-cut-and-price to optimality
#include "instance/instance.h"
#include "tree/bcp.h"
#include "tree/colgen.h"
#include "verify/brute.h"
#include "verify/verify.h"

#include "bplex.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>

using namespace cvrpbcp;

static int usage()
{
    std::fprintf(stderr,
                 "usage:\n"
                 "  cvrpbcp selftest [cases]\n"
                 "  cvrpbcp info <inst.vrp>\n"
                 "  cvrpbcp check <inst.vrp> <sol.sol>\n"
                 "  cvrpbcp root  <inst.vrp> [options]\n"
                 "  cvrpbcp solve <inst.vrp> [options] [--out sol.sol] [--jsonl results.jsonl]\n"
                 "options: --K k --ng s --knn k --cols c --alpha a --no-rcc --no-early --time T --log L\n"
                 "         --ub u --no-rmip --rmip-time s --rmip-every n --nodes n\n"
                 "         --purge-age a --purge-min c --sb k --sb-iters i --sb-rel r\n"
                 "         --free-fleet (ignore the instance's vehicle count) --fleet exact|max\n"
                 "         --no-src --src-max m\n");
    return 2;
}

static void parse_opts(int argc, char **argv, int from, BcpParams &P, std::string *out, std::string *jsonl)
{
    for (int a = from; a < argc; a++) {
        std::string o = argv[a];
        auto val = [&]() {
            if (a + 1 >= argc) throw std::runtime_error("missing value for " + o);
            return std::string(argv[++a]);
        };
        if (o == "--K") P.kmax = std::stoi(val());
        else if (o == "--ng") P.ng_size = std::stoi(val());
        else if (o == "--knn") P.heur_knn = std::stoi(val());
        else if (o == "--cols") P.max_cols_per_iter = std::stoi(val());
        else if (o == "--alpha") P.smoothing_alpha = std::stod(val());
        else if (o == "--alpha-fixed") P.smoothing_auto = false;
        else if (o == "--no-rcc") P.rcc = false;
        else if (o == "--no-early") P.early_stop = false;
        else if (o == "--time") P.time_limit = std::stod(val());
        else if (o == "--log") P.log_level = std::stoi(val());
        else if (o == "--ub") P.ub = std::stod(val());
        else if (o == "--no-rmip") P.rmip = false;
        else if (o == "--rmip-time") P.rmip_time = std::stod(val());
        else if (o == "--rmip-every") P.rmip_every = std::stoi(val());
        else if (o == "--nodes") P.node_limit = std::stoi(val());
        else if (o == "--purge-age") P.purge_age = std::stoi(val());
        else if (o == "--purge-min") P.purge_min_cols_per_cust = std::stoi(val());
        else if (o == "--sb") P.sb_cands = std::stoi(val());
        else if (o == "--sb-iters") P.sb_iters = std::stoi(val());
        else if (o == "--sb-rel") P.sb_reliable = std::stoi(val());
        else if (o == "--threads") P.threads = std::stoi(val());
        else if (o == "--free-fleet") P.free_fleet = true;
        else if (o == "--fleet") { std::string m = val(); P.fleet_mode = (m == "max") ? 1 : 0; }
        else if (o == "--no-src") P.src = false;
        else if (o == "--mono") P.bidirectional = false;
        else if (o == "--no-savings") P.savings_init = false;
        else if (o == "--cut-age") P.cut_age_limit = std::stoi(val());
        else if (o == "--no-dive") P.dive = false;
        else if (o == "--dive-time") P.dive_time = std::stod(val());
        else if (o == "--plunge") P.plunge_depth = std::stoi(val());
        else if (o == "--heur-cap") P.heur_label_cap = std::stoi(val());
        else if (o == "--qbound") P.qbound_stage = std::stoi(val());
        else if (o == "--no-enum") P.enumerate = false;
        else if (o == "--enum-time") P.enum_mip_time = std::stod(val());
        else if (o == "--enum-mip-max") P.enum_mip_max_routes = (size_t)std::stoll(val());
        else if (o == "--enum-dump") P.enum_dump = val();
        else if (o == "--src-max") P.src_max_active = std::stoi(val());
        else if (o == "--out" && out) *out = val();
        else if (o == "--jsonl" && jsonl) *jsonl = val();
        else throw std::runtime_error("unknown option " + o);
    }
}

static int cmd_selftest(int cases)
{
    int fails = 0, done = 0;
    for (int s = 1; s <= cases; s++) {
        bool tw = (s % 3 == 0);
        int n = 5 + (s % 4);
        Instance I = random_instance(n, 10, 4, 1000 + s, tw);
        int K = (s % 2 == 0) ? 0 : I.ffd_vehicles() + 1;
        Solution a = brute_force_optimum(I, K);
        Solution b = bplex_setpart_optimum(I, K);
        CheckResult ca = check_solution(I, a), cb = check_solution(I, b);
        bool ok = a.cost >= 0 && b.cost >= 0 && ca.feasible && cb.feasible && ca.cost == a.cost &&
                  cb.cost == b.cost && a.cost == b.cost;
        {
            BcpParams P;
            P.log_level = 0;
            P.kmax = K;
            BcpResult R = solve_bcp(I, P);
            if (R.status != "optimal" || R.ub != a.cost) {
                ok = false;
                std::printf("     bcp: status %s lb %lld ub %lld nodes %d err=%s\n", R.status.c_str(), R.lb, R.ub,
                            R.nodes, R.error.c_str());
            }
        }
        done++;
        if (!ok) {
            fails++;
            std::printf("FAIL %s K=%d: brute %lld (%s) bplex %lld (%s)\n", I.name.c_str(), K, a.cost,
                        ca.feasible ? "feasible" : ca.error.c_str(), b.cost,
                        cb.feasible ? "feasible" : cb.error.c_str());
        } else {
            std::printf("ok   %-14s n=%d tw=%d K=%d opt=%lld routes=%zu\n", I.name.c_str(), n, (int)tw, K,
                        a.cost, a.routes.size());
        }
    }
    std::printf("selftest: %d/%d passed\n", done - fails, done);
    return fails ? 1 : 0;
}

static int cmd_info(const char *path)
{
    Instance I = read_instance(path);
    std::printf("%s\n", I.summary().c_str());
    std::printf("bplex %s\n", bplex_version());
    long long opt = comment_optimum(path);
    if (opt >= 0) std::printf("comment optimum: %lld\n", opt);
    return 0;
}

static int cmd_check(const char *ipath, const char *spath)
{
    Instance I = read_instance(ipath);
    Solution S = read_sol(spath, I.scale);
    CheckResult R = check_solution(I, S, false);
    long long opt = comment_optimum(ipath);
    if (!R.feasible) {
        std::printf("%s: INFEASIBLE: %s\n", I.name.c_str(), R.error.c_str());
        return 1;
    }
    bool match = (S.cost < 0 || S.cost == R.cost);
    std::printf("%s: feasible, %d routes, recomputed cost %lld, file cost %lld%s%s\n", I.name.c_str(),
                R.vehicles, R.cost, S.cost, match ? "" : " MISMATCH",
                (opt >= 0 && opt * I.scale != R.cost) ? " (differs from COMMENT optimum)" : "");
    return match ? 0 : 1;
}

static int cmd_root(int argc, char **argv)
{
    BcpParams P;
    parse_opts(argc, argv, 3, P, nullptr, nullptr);
    Instance I = read_instance(argv[2]);
    if (P.free_fleet) I.vehicles = -1;
    long long opt = comment_optimum(argv[2]);
    NodeResult R = solve_root(I, P, P.kmax);
    if (!R.ok) {
        std::printf("ROOT %s FAILED: %s\n", I.name.c_str(), R.error.c_str());
        return 1;
    }
    std::printf("ROOT %s lp %.4f lb %.4f lbint %lld opt %lld gap %.3f%% iters %d cutrounds %d cols %d rows %d "
                "integral %d fleetbind %d time %.2fs (pricing %.2f master %.2f cuts %.2f of which src %.2f)\n",
                I.name.c_str(), R.lp, R.lb, R.lb_int, opt, opt > 0 ? 100.0 * (opt - R.lb) / opt : 0.0,
                R.cg_iters, R.cut_rounds, R.ncols, R.ncuts, (int)R.integral, (int)R.fleet_binding, R.time,
                R.time_pricing, R.time_master, R.time_cuts, R.time_src);
    if (R.integral) std::printf("%s", format_sol(R.sol, I.scale).c_str());
    if (opt > 0 && R.lb_int > opt) { std::printf("BUG: lower bound exceeds the published optimum\n"); return 1; }
    return 0;
}

static int cmd_solve(int argc, char **argv)
{
    BcpParams P;
    std::string out, jsonl;
    parse_opts(argc, argv, 3, P, &out, &jsonl);
    Instance I = read_instance(argv[2]);
    if (P.free_fleet) I.vehicles = -1;
    long long opt = comment_optimum(argv[2]);
    // published optimum from a .sol next to the instance, if any
    std::string spath = argv[2];
    size_t dot = spath.rfind('.');
    if (dot != std::string::npos) spath = spath.substr(0, dot) + ".sol";
    long long published = opt > 0 ? opt * I.scale : -1;
    if (FILE *f = std::fopen(spath.c_str(), "r")) {
        std::fclose(f);
        try {
            Solution S = read_sol(spath, I.scale);
            CheckResult c = check_solution(I, S);
            if (c.feasible) published = c.cost;
        } catch (...) {
        }
    }
    BcpResult R = solve_bcp(I, P);
    const char *verdict = "";
    if (R.status == "optimal") {
        if (published >= 0) verdict = (R.ub == published) ? "MATCHES published optimum" : "DIFFERS from published optimum";
        else verdict = "no published optimum to compare";
    }
    std::printf("RESULT %s status %s lb %lld ub %lld published %lld nodes %d depth %d cols %d rows %d rootlp %.3f "
                "rootlb %.3f kmin %d kmax %d fleetcap %d time %.2fs root %.2fs %s\n",
                I.name.c_str(), R.status.c_str(), R.lb, R.ub, published, R.nodes, R.depth_max, R.ncols, R.nrows,
                R.root_lp, R.root_lb, R.kmin, R.kmax, (int)R.fleet_cap_binding, R.time, R.time_root, verdict);
    if (!R.error.empty()) std::printf("ERROR: %s\n", R.error.c_str());
    if (R.ub >= 0) std::printf("%s", format_sol(R.sol, I.scale).c_str());
    if (!out.empty() && R.ub >= 0) write_sol(out, R.sol, I.scale);
    if (!jsonl.empty()) {
        if (FILE *f = std::fopen(jsonl.c_str(), "a")) {
            std::fprintf(f,
                         "{\"instance\":\"%s\",\"n\":%d,\"status\":\"%s\",\"lb\":%lld,\"ub\":%lld,\"published\":%lld,"
                         "\"nodes\":%d,\"cols\":%d,\"rows\":%d,\"root_lp\":%.4f,\"root_lb\":%.4f,\"kmin\":%d,\"kmax\":%d,"
                         "\"fleet_cap_binding\":%s,\"time\":%.2f,\"time_root\":%.2f,\"scale\":%d,"
                         "\"bplex\":\"%s\"}\n",
                         I.name.c_str(), I.n, R.status.c_str(), R.lb, R.ub, published, R.nodes, R.ncols, R.nrows,
                         R.root_lp, R.root_lb, R.kmin, R.kmax, R.fleet_cap_binding ? "true" : "false", R.time,
                         R.time_root, I.scale, bplex_version());
            std::fclose(f);
        }
    }
    if (R.status == "error") return 3;
    if (R.status == "optimal" && published >= 0 && R.ub != published) return 1;
    if (R.status != "optimal" && published >= 0 && R.lb > published) { std::printf("BUG: lower bound exceeds the published optimum\n"); return 1; }
    return 0;
}

int main(int argc, char **argv)
{
    std::setvbuf(stdout, nullptr, _IOLBF, 0); // line-buffered even when piped: logs survive a kill
    if (argc < 2) return usage();
    try {
        std::string cmd = argv[1];
        if (cmd == "selftest") return cmd_selftest(argc > 2 ? std::atoi(argv[2]) : 12);
        if (cmd == "info" && argc == 3) return cmd_info(argv[2]);
        if (cmd == "check" && argc == 4) return cmd_check(argv[2], argv[3]);
        if (cmd == "root" && argc >= 3) return cmd_root(argc, argv);
        if (cmd == "solve" && argc >= 3) return cmd_solve(argc, argv);
        return usage();
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 3;
    }
}
