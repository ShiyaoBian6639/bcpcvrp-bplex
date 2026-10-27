/* Standalone check of bplex_chg_row_bounds on a ranged row whose slack is
 * basic (the cvrpbcp fleet row after a fleet branch):
 *   min -x1 - x2,  0 <= x1 <= 0.6, 0 <= x2 <= 0.5,  row: x1 + x2 in [lo, up]
 * Solve with [0, 2] (activity 1.1, slack basic), then tighten the upper bound
 * to 1.0 and re-solve: the optimum must move to activity 1.0, objective -1.0.
 * Then [1, 1] (both bounds), then back to [0, 2]. Exit 0 when all re-solves
 * respect the bounds, 1 otherwise. Build (WSL, after scripts/build_wsl.sh):
 *   gcc -O1 tests/bplex_rowbounds.c -Ibuild/bplex-install/include \
 *       build/bplex-install/lib/libbplex.a build/bplex-install/lib/libbasiclu.a -lm -lpthread -o build/bplex_rowbounds
 */
#include "bplex.h"

#include <math.h>
#include <stdio.h>

static int check(bplex_t *bp, const char *what, double lo, double up, double want_obj)
{
    int st = bplex_solve(bp);
    double x[2] = {0, 0}, act = 0;
    bplex_get_solution(bp, x, &act, NULL, NULL);
    double obj = bplex_objval(bp);
    int it = (int)bplex_get_info(bp, "iterations");
    int ok = st == BPLEX_OPTIMAL && act <= up + 1e-7 && act >= lo - 1e-7 && fabs(obj - want_obj) < 1e-6;
    printf("%-22s status %d iters %d x (%.3f %.3f) activity %.4f bounds [%.1f %.1f] obj %.4f (want %.4f) %s\n", what,
           st, it, x[0], x[1], act, lo, up, obj, want_obj, ok ? "ok" : "WRONG");
    return ok;
}

int main(void)
{
    printf("bplex %s\n", bplex_version());
    bplex_t *bp = bplex_create();
    bplex_set_param(bp, "log_level", 0);
    bplex_set_param(bp, "presolve", 0);
    double obj[2] = {-1, -1}, lo[2] = {0, 0}, up[2] = {0.6, 0.5}, rlo[1] = {0}, rup[1] = {2};
    int cs[3] = {0, 1, 2}, ri[2] = {0, 0};
    double v[2] = {1, 1};
    if (bplex_load(bp, 2, 1, obj, lo, up, NULL, rlo, rup, cs, ri, v, 0.0, 1) != 0) return 2;
    int ok = 1;
    ok &= check(bp, "initial [0,2]", 0, 2, -1.1);
    if (bplex_chg_row_bounds(bp, 0, 0, 1.0) != 0) return 2;
    ok &= check(bp, "upper -> 1.0", 0, 1.0, -1.0);
    if (bplex_chg_row_bounds(bp, 0, 0, 2.0) != 0) return 2;
    ok &= check(bp, "back to [0,2]", 0, 2, -1.1);
    if (bplex_chg_row_bounds(bp, 0, 1.0, 1.0) != 0) return 2;
    ok &= check(bp, "both -> [1,1]", 1.0, 1.0, -1.0);
    if (bplex_chg_row_bounds(bp, 0, 0.3, 0.3) != 0) return 2;
    ok &= check(bp, "both -> [0.3,0.3]", 0.3, 0.3, -0.3);
    /* the fleet-row pattern: slack basic in the interior, then only the upper bound moves below the activity */
    if (bplex_chg_row_bounds(bp, 0, 0, 2.0) != 0) return 2;
    ok &= check(bp, "again [0,2]", 0, 2, -1.1);
    if (bplex_chg_row_bounds(bp, 0, 0, 1.05) != 0) return 2;
    ok &= check(bp, "upper -> 1.05", 0, 1.05, -1.05);
    bplex_free(bp);
    printf(ok ? "rowbounds: all ok\n" : "rowbounds: FAILED\n");
    return ok ? 0 : 1;
}
