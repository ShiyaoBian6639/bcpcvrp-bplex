// ctest entry: brute force vs bplex set partitioning vs branch-cut-and-price
// on random tiny instances, plus reader/checker round trips on two CVRPLIB
// files (run from the repo root).
#include "instance/instance.h"
#include "tree/bcp.h"
#include "verify/brute.h"
#include "verify/verify.h"

#include <cstdio>
#include <string>

using namespace cvrpbcp;

static int fails = 0;
#define CHECK(c, msg)                                                                                                  \
    do {                                                                                                               \
        if (!(c)) { std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, msg); fails++; }                                \
    } while (0)

int main()
{
    for (int s = 1; s <= 10; s++) {
        bool tw = (s % 3 == 0);
        int n = 5 + (s % 4);
        Instance I = random_instance(n, 10, 4, 2000 + s, tw);
        int K = (s % 2 == 0) ? 0 : I.ffd_vehicles() + 1;
        Solution a = brute_force_optimum(I, K);
        Solution b = bplex_setpart_optimum(I, K);
        CHECK(a.cost >= 0 && a.cost == b.cost, "brute force and bplex set partitioning disagree");
        CHECK(check_solution(I, a).feasible, "brute force solution infeasible");
        {
            BcpParams P;
            P.log_level = 0;
            P.kmax = K;
            BcpResult R = solve_bcp(I, P);
            CHECK(R.status == "optimal" && R.ub == a.cost, "bcp optimum differs from brute force");
        }
    }
    const char *files[] = {"instances/cvrplib/A/A-n32-k5", "instances/cvrplib/E/E-n13-k4"};
    for (const char *f : files) {
        Instance I = read_instance(std::string(f) + ".vrp");
        Solution S = read_sol(std::string(f) + ".sol", I.scale);
        CheckResult c = check_solution(I, S);
        CHECK(c.feasible && c.cost == S.cost && c.cost == comment_optimum(std::string(f) + ".vrp"),
              "published solution does not verify");
        BcpParams P;
        P.log_level = 0;
        BcpResult R = solve_bcp(I, P);
        CHECK(R.status == "optimal" && R.ub == S.cost, "bcp does not reproduce the published optimum");
    }
    std::printf(fails ? "test_small: %d failures\n" : "test_small: all passed\n", fails);
    return fails ? 1 : 0;
}
