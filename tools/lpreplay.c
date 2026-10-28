/* lpreplay: replays an LP trace written by cvrpbcp (src/master/lptrace.h)
 * on an LP solver and reports the solve times, so that the master's LP work
 * can be compared between solvers on the identical operation sequence
 * (warm starts included).
 *
 *   lpreplay <trace> [max_solves]
 *
 * Backends (compile-time):
 *   default         bplex   (link bplex::bplex; built by the cvrpbcp CMake)
 *   -DLPREPLAY_CPLEX CPLEX  (functions loaded at run time from cplex2210.dll
 *                            or libcplex2210.so: no import library needed;
 *                            needs ilcplex/cplex.h for the constants)
 *     Windows, no MSVC:  python -m ziglang cc -O2 -DLPREPLAY_CPLEX
 *        -I"%CPLEX%\include" tools\lpreplay.c -o lpreplay_cplex.exe
 *     then run with cplex2210.dll on PATH or LPREPLAY_CPLEX_DLL=<path>.
 *
 * Output: one line per solve kind (LP, MIP, probe) with count, recorded
 * seconds (the solver that produced the trace), replayed seconds, iteration
 * totals, and the largest objective disagreement among solves both ended
 * optimal. Per-solve lines with -v.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void)
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---------------- trace reader ---------------- */
enum { LOAD = 1, ADDCOLS, DELCOLS, ADDROWS, DELROWS, CHGROW, CHGCOL, CHGCOLS, CHGTYPE, PARAM, SOLVE, RESULT, SLACKBASIS };
static FILE *in;
static int rep_ncols = 0, rep_nrows = 0; /* current sizes, for the slack basis */
static int rd_i(void)
{
    int x;
    if (fread(&x, sizeof x, 1, in) != 1) { fprintf(stderr, "truncated trace\n"); exit(2); }
    return x;
}
static double rd_d(void)
{
    double x;
    if (fread(&x, sizeof x, 1, in) != 1) { fprintf(stderr, "truncated trace\n"); exit(2); }
    return x;
}
static void *rd_arr(size_t sz, int *n)
{
    *n = rd_i();
    void *p = malloc(sz * (size_t)(*n > 0 ? *n : 1));
    if (*n > 0 && fread(p, sz, (size_t)*n, in) != (size_t)*n) { fprintf(stderr, "truncated trace\n"); exit(2); }
    return p;
}

/* ---------------- backend interface ---------------- */
typedef struct {
    void (*load)(int nc, int nr, const double *obj, const double *lo, const double *up, const double *rlo,
                 const double *rup, const int *cs, const int *ri, const double *v);
    void (*add_cols)(int n, const double *obj, const double *lo, const double *up, const int *cs, const int *ri,
                     const double *v);
    void (*del_cols)(int n, const int *idx);
    void (*add_rows)(int n, const double *lo, const double *up, const int *rs, const int *ci, const double *v);
    void (*del_rows)(int n, const int *idx);
    void (*chg_row)(int i, double lo, double up);
    void (*chg_col)(int j, double lo, double up);
    void (*chg_cols)(int n, const int *idx, const double *lo, const double *up);
    void (*chg_type)(int j, int is_int);
    void (*param)(const char *name, double v);
    void (*slack_basis)(void);
    /* returns 1 if optimal; fills obj and iterations */
    int (*solve)(int kind, double *obj, double *iters);
    const char *name;
} Backend;

#ifndef LPREPLAY_CPLEX
/* ======================= bplex backend ======================= */
#include "bplex.h"
static bplex_t *bp;
static void chk(int rc, const char *what)
{
    if (rc != 0) { fprintf(stderr, "bplex %s failed: %d\n", what, rc); exit(3); }
}
static void b_load(int nc, int nr, const double *obj, const double *lo, const double *up, const double *rlo,
                   const double *rup, const int *cs, const int *ri, const double *v)
{
    chk(bplex_load(bp, nc, nr, obj, lo, up, NULL, rlo, rup, cs, ri, v, 0.0, 1), "load");
}
static void b_add_cols(int n, const double *obj, const double *lo, const double *up, const int *cs, const int *ri,
                       const double *v)
{
    chk(bplex_add_cols(bp, n, obj, lo, up, cs, ri, v), "add_cols");
}
static void b_del_cols(int n, const int *idx) { chk(bplex_del_cols(bp, n, idx), "del_cols"); }
static void b_add_rows(int n, const double *lo, const double *up, const int *rs, const int *ci, const double *v)
{
    chk(bplex_add_rows(bp, n, lo, up, rs, ci, v), "add_rows");
}
static void b_del_rows(int n, const int *idx) { chk(bplex_del_rows(bp, n, idx), "del_rows"); }
static void b_chg_row(int i, double lo, double up) { chk(bplex_chg_row_bounds(bp, i, lo, up), "chg_row_bounds"); }
static void b_chg_col(int j, double lo, double up) { chk(bplex_chg_col_bounds(bp, j, lo, up), "chg_col_bounds"); }
static void b_chg_cols(int n, const int *idx, const double *lo, const double *up)
{
    chk(bplex_chg_col_bounds_many(bp, n, idx, lo, up), "chg_col_bounds_many");
}
static void b_chg_type(int j, int is_int) { chk(bplex_chg_col_type(bp, j, is_int), "chg_col_type"); }
static void b_param(const char *name, double v) { bplex_set_param(bp, name, v); }
static void b_slack_basis(void)
{
    int nc = rep_ncols, nr = rep_nrows;
    int *cst = malloc(sizeof(int) * (size_t)(nc > 0 ? nc : 1)), *rst = malloc(sizeof(int) * (size_t)(nr > 0 ? nr : 1));
    for (int j = 0; j < nc; j++) cst[j] = BPLEX_AT_LOWER;
    for (int i = 0; i < nr; i++) rst[i] = BPLEX_BASIC;
    chk(bplex_set_basis(bp, cst, rst), "set_basis");
    free(cst); free(rst);
}
static int b_solve(int kind, double *obj, double *iters)
{
    (void)kind;
    int st = bplex_solve(bp);
    *iters = bplex_get_info(bp, "iterations");
    *obj = (st == BPLEX_OPTIMAL) ? bplex_objval(bp) : 0.0;
    return st == BPLEX_OPTIMAL;
}
static Backend make_backend(void)
{
    bp = bplex_create();
    if (!bp) { fprintf(stderr, "bplex_create failed\n"); exit(3); }
    Backend B = {b_load, b_add_cols, b_del_cols, b_add_rows, b_del_rows, b_chg_row, b_chg_col,
                 b_chg_cols, b_chg_type, b_param, b_slack_basis, b_solve, "bplex"};
    return B;
}
#else
/* ======================= CPLEX backend (dynamic) ======================= */
#include <ilcplex/cplex.h>
#ifdef _WIN32
#include <windows.h>
static void *dl_open(const char *p) { return (void *)LoadLibraryA(p); }
static void *dl_sym(void *h, const char *s) { return (void *)GetProcAddress((HMODULE)h, s); }
#else
#include <dlfcn.h>
static void *dl_open(const char *p) { return dlopen(p, RTLD_NOW); }
static void *dl_sym(void *h, const char *s) { return dlsym(h, s); }
#endif
typedef CPXENVptr(CPXPUBLIC *f_open)(int *);
typedef CPXLPptr(CPXPUBLIC *f_createprob)(CPXCENVptr, int *, const char *);
typedef int(CPXPUBLIC *f_copylp)(CPXCENVptr, CPXLPptr, int, int, int, const double *, const double *, const char *,
                                 const int *, const int *, const int *, const double *, const double *,
                                 const double *, const double *);
typedef int(CPXPUBLIC *f_addcols)(CPXCENVptr, CPXLPptr, int, int, const double *, const int *, const int *,
                                  const double *, const double *, const double *, char **);
typedef int(CPXPUBLIC *f_addrows)(CPXCENVptr, CPXLPptr, int, int, int, const double *, const char *, const int *,
                                  const int *, const double *, char **, char **);
typedef int(CPXPUBLIC *f_delset)(CPXCENVptr, CPXLPptr, int *);
typedef int(CPXPUBLIC *f_chgbds)(CPXCENVptr, CPXLPptr, int, const int *, const char *, const double *);
typedef int(CPXPUBLIC *f_chgrhs)(CPXCENVptr, CPXLPptr, int, const int *, const double *);
typedef int(CPXPUBLIC *f_chgsense)(CPXCENVptr, CPXLPptr, int, const int *, const char *);
typedef int(CPXPUBLIC *f_chgctype)(CPXCENVptr, CPXLPptr, int, const int *, const char *);
typedef int(CPXPUBLIC *f_chgprobtype)(CPXCENVptr, CPXLPptr, int);
typedef int(CPXPUBLIC *f_opt)(CPXCENVptr, CPXLPptr);
typedef int(CPXPUBLIC *f_getstat)(CPXCENVptr, CPXCLPptr);
typedef int(CPXPUBLIC *f_getobjval)(CPXCENVptr, CPXCLPptr, double *);
typedef int(CPXPUBLIC *f_getitcnt)(CPXCENVptr, CPXCLPptr);
typedef int(CPXPUBLIC *f_getprobtype)(CPXCENVptr, CPXCLPptr);
typedef int(CPXPUBLIC *f_setint)(CPXENVptr, int, CPXINT);
typedef int(CPXPUBLIC *f_setdbl)(CPXENVptr, int, double);
typedef int(CPXPUBLIC *f_getnum)(CPXCENVptr, CPXCLPptr);
typedef int(CPXPUBLIC *f_copybase)(CPXCENVptr, CPXLPptr, const int *, const int *);
static struct {
    f_open open;
    f_createprob createprob;
    f_copylp copylp;
    f_addcols addcols;
    f_addrows addrows;
    f_delset delsetcols, delsetrows;
    f_chgbds chgbds;
    f_chgrhs chgrhs, chgrngval;
    f_chgsense chgsense;
    f_chgctype chgctype;
    f_chgprobtype chgprobtype;
    f_opt lpopt, mipopt;
    f_getstat getstat;
    f_getobjval getobjval;
    f_getitcnt getitcnt, getmipitcnt;
    f_getprobtype getprobtype;
    f_setint setint;
    f_setdbl setdbl;
    f_getnum getnumcols, getnumrows;
    f_copybase copybase;
} C;
static CPXENVptr env;
static CPXLPptr lp;
static void cchk(int rc, const char *what)
{
    if (rc != 0) { fprintf(stderr, "CPLEX %s failed: %d\n", what, rc); exit(3); }
}
static void *need(void *h, const char *s)
{
    void *p = dl_sym(h, s);
    if (!p) { fprintf(stderr, "symbol %s not found in the CPLEX library\n", s); exit(3); }
    return p;
}
static void row_sense(double lo, double up, char *sense, double *rhs, double *rng)
{
    *rng = 0;
    if (lo <= -1e29 && up >= 1e29) { *sense = 'R'; *rhs = -1e20; *rng = 2e20; }
    else if (fabs(lo - up) < 1e-12) { *sense = 'E'; *rhs = lo; }
    else if (lo <= -1e29) { *sense = 'L'; *rhs = up; }
    else if (up >= 1e29) { *sense = 'G'; *rhs = lo; }
    else { *sense = 'R'; *rhs = lo; *rng = up - lo; }
}
static void c_load(int nc, int nr, const double *obj, const double *lo, const double *up, const double *rlo,
                   const double *rup, const int *cs, const int *ri, const double *v)
{
    char *sense = malloc((size_t)nr);
    double *rhs = malloc(sizeof(double) * (size_t)nr), *rng = malloc(sizeof(double) * (size_t)nr);
    int *cnt = malloc(sizeof(int) * (size_t)nc);
    for (int i = 0; i < nr; i++) row_sense(rlo[i], rup[i], &sense[i], &rhs[i], &rng[i]);
    for (int j = 0; j < nc; j++) cnt[j] = cs[j + 1] - cs[j];
    cchk(C.copylp(env, lp, nc, nr, CPX_MIN, obj, rhs, sense, cs, cnt, ri, v, lo, up, rng), "copylp");
    free(sense); free(rhs); free(rng); free(cnt);
}
static void c_add_cols(int n, const double *obj, const double *lo, const double *up, const int *cs, const int *ri,
                       const double *v)
{
    cchk(C.addcols(env, lp, n, cs[n], obj, cs, ri, v, lo, up, NULL), "addcols");
}
static void delset(f_delset f, int total, int n, const int *idx, const char *what)
{
    int *del = calloc((size_t)(total > 0 ? total : 1), sizeof(int));
    for (int k = 0; k < n; k++) del[idx[k]] = 1;
    cchk(f(env, lp, del), what);
    free(del);
}
static void c_del_cols(int n, const int *idx) { delset(C.delsetcols, C.getnumcols(env, lp), n, idx, "delsetcols"); }
static void c_add_rows(int n, const double *lo, const double *up, const int *rs, const int *ci, const double *v)
{
    char *sense = malloc((size_t)n);
    double *rhs = malloc(sizeof(double) * (size_t)n), *rng = malloc(sizeof(double) * (size_t)n);
    int *ridx = malloc(sizeof(int) * (size_t)n);
    int first = C.getnumrows(env, lp), nrng = 0;
    for (int i = 0; i < n; i++) row_sense(lo[i], up[i], &sense[i], &rhs[i], &rng[i]);
    cchk(C.addrows(env, lp, 0, n, rs[n], rhs, sense, rs, ci, v, NULL, NULL), "addrows");
    for (int i = 0; i < n; i++)
        if (sense[i] == 'R') { ridx[nrng] = first + i; rng[nrng] = rng[i]; nrng++; }
    if (nrng) cchk(C.chgrngval(env, lp, nrng, ridx, rng), "chgrngval");
    free(sense); free(rhs); free(rng); free(ridx);
}
static void c_del_rows(int n, const int *idx) { delset(C.delsetrows, C.getnumrows(env, lp), n, idx, "delsetrows"); }
static void c_chg_row(int i, double lo, double up)
{
    char sense;
    double rhs, rng;
    row_sense(lo, up, &sense, &rhs, &rng);
    cchk(C.chgsense(env, lp, 1, &i, &sense), "chgsense");
    cchk(C.chgrhs(env, lp, 1, &i, &rhs), "chgrhs");
    if (sense == 'R') cchk(C.chgrngval(env, lp, 1, &i, &rng), "chgrngval");
}
static void c_chg_cols(int n, const int *idx, const double *lo, const double *up)
{
    int *ind = malloc(sizeof(int) * (size_t)(2 * n));
    char *lu = malloc((size_t)(2 * n));
    double *bd = malloc(sizeof(double) * (size_t)(2 * n));
    for (int k = 0; k < n; k++) {
        ind[2 * k] = idx[k]; lu[2 * k] = 'L'; bd[2 * k] = lo[k];
        ind[2 * k + 1] = idx[k]; lu[2 * k + 1] = 'U'; bd[2 * k + 1] = up[k] >= 1e29 ? CPX_INFBOUND : up[k];
    }
    cchk(C.chgbds(env, lp, 2 * n, ind, lu, bd), "chgbds");
    free(ind); free(lu); free(bd);
}
static void c_chg_col(int j, double lo, double up) { c_chg_cols(1, &j, &lo, &up); }
static void c_chg_type(int j, int is_int)
{
    char t = is_int ? 'I' : 'C';
    cchk(C.chgctype(env, lp, 1, &j, &t), "chgctype");
}
static void c_param(const char *name, double v)
{
    if (!strcmp(name, "time_limit")) C.setdbl(env, CPX_PARAM_TILIM, v >= 1e29 ? 1e75 : v);
    else if (!strcmp(name, "iter_limit")) C.setint(env, CPX_PARAM_ITLIM, v >= 2e9 ? 2100000000 : (CPXINT)v);
    else if (!strcmp(name, "threads")) C.setint(env, CPX_PARAM_THREADS, (CPXINT)(v < 1 ? 1 : v));
    else if (!strcmp(name, "method")) C.setint(env, CPX_PARAM_LPMETHOD, (int)v == 2 ? CPX_ALG_DUAL : CPX_ALG_AUTOMATIC);
    else if (!strcmp(name, "cutoff")) C.setdbl(env, CPX_PARAM_CUTUP, fabs(v) >= 1e20 ? 1e75 : v);
    else if (!strcmp(name, "mip_abs_gap")) C.setdbl(env, CPX_PARAM_EPAGAP, v);
    else if (!strcmp(name, "mip_rel_gap")) C.setdbl(env, CPX_PARAM_EPGAP, v);
    else if (!strcmp(name, "presolve")) C.setint(env, CPX_PARAM_PREIND, v != 0);
    /* log_level, relax: handled by the solve kind */
}
static void c_slack_basis(void)
{
    int nc = C.getnumcols(env, lp), nr = C.getnumrows(env, lp);
    int *cst = malloc(sizeof(int) * (size_t)(nc > 0 ? nc : 1)), *rst = malloc(sizeof(int) * (size_t)(nr > 0 ? nr : 1));
    for (int j = 0; j < nc; j++) cst[j] = CPX_AT_LOWER;
    for (int i = 0; i < nr; i++) rst[i] = CPX_BASIC;
    cchk(C.copybase(env, lp, cst, rst), "copybase");
    free(cst); free(rst);
}
static int c_solve(int kind, double *obj, double *iters)
{
    int st;
    if (kind == 1) {
        cchk(C.mipopt(env, lp), "mipopt");
        st = C.getstat(env, lp);
        *iters = C.getmipitcnt(env, lp);
        int ok = (st == CPXMIP_OPTIMAL || st == CPXMIP_OPTIMAL_TOL);
        *obj = 0;
        if (ok) C.getobjval(env, lp, obj);
        return ok;
    }
    if (C.getprobtype(env, lp) != CPXPROB_LP) cchk(C.chgprobtype(env, lp, CPXPROB_LP), "chgprobtype");
    cchk(C.lpopt(env, lp), "lpopt");
    st = C.getstat(env, lp);
    *iters = C.getitcnt(env, lp);
    *obj = 0;
    if (st == CPX_STAT_OPTIMAL) C.getobjval(env, lp, obj);
    return st == CPX_STAT_OPTIMAL;
}
static Backend make_backend(void)
{
    const char *path = getenv("LPREPLAY_CPLEX_DLL");
#ifdef _WIN32
    if (!path) path = "cplex2210.dll";
#else
    if (!path) path = "libcplex2210.so";
#endif
    void *h = dl_open(path);
    if (!h) { fprintf(stderr, "cannot load %s (set LPREPLAY_CPLEX_DLL)\n", path); exit(3); }
    C.open = (f_open)need(h, "CPXopenCPLEX");
    C.createprob = (f_createprob)need(h, "CPXcreateprob");
    C.copylp = (f_copylp)need(h, "CPXcopylp");
    C.addcols = (f_addcols)need(h, "CPXaddcols");
    C.addrows = (f_addrows)need(h, "CPXaddrows");
    C.delsetcols = (f_delset)need(h, "CPXdelsetcols");
    C.delsetrows = (f_delset)need(h, "CPXdelsetrows");
    C.chgbds = (f_chgbds)need(h, "CPXchgbds");
    C.chgrhs = (f_chgrhs)need(h, "CPXchgrhs");
    C.chgrngval = (f_chgrhs)need(h, "CPXchgrngval");
    C.chgsense = (f_chgsense)need(h, "CPXchgsense");
    C.chgctype = (f_chgctype)need(h, "CPXchgctype");
    C.chgprobtype = (f_chgprobtype)need(h, "CPXchgprobtype");
    C.lpopt = (f_opt)need(h, "CPXlpopt");
    C.mipopt = (f_opt)need(h, "CPXmipopt");
    C.getstat = (f_getstat)need(h, "CPXgetstat");
    C.getobjval = (f_getobjval)need(h, "CPXgetobjval");
    C.getitcnt = (f_getitcnt)need(h, "CPXgetitcnt");
    C.getmipitcnt = (f_getitcnt)need(h, "CPXgetmipitcnt");
    C.getprobtype = (f_getprobtype)need(h, "CPXgetprobtype");
    C.setint = (f_setint)need(h, "CPXsetintparam");
    C.setdbl = (f_setdbl)need(h, "CPXsetdblparam");
    C.getnumcols = (f_getnum)need(h, "CPXgetnumcols");
    C.getnumrows = (f_getnum)need(h, "CPXgetnumrows");
    C.copybase = (f_copybase)need(h, "CPXcopybase");
    int st = 0;
    env = C.open(&st);
    if (!env) { fprintf(stderr, "CPXopenCPLEX failed: %d\n", st); exit(3); }
    lp = C.createprob(env, &st, "replay");
    C.setint(env, CPX_PARAM_SCRIND, 0);
    C.setint(env, CPX_PARAM_THREADS, 1);
    Backend B = {c_load, c_add_cols, c_del_cols, c_add_rows, c_del_rows, c_chg_row, c_chg_col,
                 c_chg_cols, c_chg_type, c_param, c_slack_basis, c_solve, "cplex"};
    return B;
}
#endif

/* ---------------- replay ---------------- */
int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: lpreplay <trace> [max_solves] [-v]\n"); return 1; }
    long max_solves = argc > 2 && argv[2][0] != '-' ? atol(argv[2]) : 0;
    int verbose = 0;
    for (int a = 2; a < argc; a++) if (!strcmp(argv[a], "-v")) verbose = 1;
    in = fopen(argv[1], "rb");
    if (!in) { perror(argv[1]); return 1; }
    char magic[8];
    if (fread(magic, 1, 8, in) != 8 || memcmp(magic, "CVRPLPT1", 8)) { fprintf(stderr, "not a cvrpbcp LP trace\n"); return 1; }
    Backend B = make_backend();
    static const char *kname[3] = {"LP", "MIP", "probe"};
    long nsolve[3] = {0, 0, 0}, nboth[3] = {0, 0, 0}, nmismatch[3] = {0, 0, 0};
    double trec[3] = {0, 0, 0}, trep[3] = {0, 0, 0}, irec[3] = {0, 0, 0}, irep[3] = {0, 0, 0}, maxdev[3] = {0, 0, 0};
    int kind = 0, ok = 0;
    double obj = 0, iters = 0, tsolve = 0;
    long solves = 0;
    double t_total = now_s();
    for (;;) {
        int op = fgetc(in);
        if (op == EOF) break;
        int n, n2, n3;
        switch (op) {
        case LOAD: {
            int nc = rd_i(), nr = rd_i();
            double *obj_ = rd_arr(8, &n), *lo = rd_arr(8, &n), *up = rd_arr(8, &n), *rlo = rd_arr(8, &n), *rup = rd_arr(8, &n);
            int *cs = rd_arr(4, &n), *ri = rd_arr(4, &n2);
            double *v = rd_arr(8, &n3);
            B.load(nc, nr, obj_, lo, up, rlo, rup, cs, ri, v);
            rep_ncols = nc;
            rep_nrows = nr;
            free(obj_); free(lo); free(up); free(rlo); free(rup); free(cs); free(ri); free(v);
            break;
        }
        case ADDCOLS: {
            int nc = rd_i();
            double *obj_ = rd_arr(8, &n), *lo = rd_arr(8, &n), *up = rd_arr(8, &n);
            int *cs = rd_arr(4, &n), *ri = rd_arr(4, &n2);
            double *v = rd_arr(8, &n3);
            B.add_cols(nc, obj_, lo, up, cs, ri, v);
            rep_ncols += nc;
            free(obj_); free(lo); free(up); free(cs); free(ri); free(v);
            break;
        }
        case DELCOLS: { int nd = rd_i(); int *idx = rd_arr(4, &n); B.del_cols(nd, idx); rep_ncols -= nd; free(idx); break; }
        case ADDROWS: {
            int nr = rd_i();
            double *lo = rd_arr(8, &n), *up = rd_arr(8, &n);
            int *rs = rd_arr(4, &n), *ci = rd_arr(4, &n2);
            double *v = rd_arr(8, &n3);
            B.add_rows(nr, lo, up, rs, ci, v);
            rep_nrows += nr;
            free(lo); free(up); free(rs); free(ci); free(v);
            break;
        }
        case DELROWS: { int nd = rd_i(); int *idx = rd_arr(4, &n); B.del_rows(nd, idx); rep_nrows -= nd; free(idx); break; }
        case CHGROW: { int i = rd_i(); double lo = rd_d(), up = rd_d(); B.chg_row(i, lo, up); break; }
        case CHGCOL: { int j = rd_i(); double lo = rd_d(), up = rd_d(); B.chg_col(j, lo, up); break; }
        case CHGCOLS: {
            int nc = rd_i();
            int *idx = rd_arr(4, &n);
            double *lo = rd_arr(8, &n), *up = rd_arr(8, &n);
            B.chg_cols(nc, idx, lo, up);
            free(idx); free(lo); free(up);
            break;
        }
        case CHGTYPE: { int j = rd_i(), t = rd_i(); B.chg_type(j, t); break; }
        case PARAM: {
            int len = rd_i();
            char name[64];
            if (len >= (int)sizeof name) { fprintf(stderr, "bad param name\n"); return 2; }
            if (fread(name, 1, (size_t)len, in) != (size_t)len) return 2;
            name[len] = 0;
            double v = rd_d();
            B.param(name, v);
            break;
        }
        case SLACKBASIS: B.slack_basis(); break;
        case SOLVE: {
            kind = rd_i();
            if (kind < 0 || kind > 2) kind = 0;
            double t = now_s();
            ok = B.solve(kind, &obj, &iters);
            tsolve = now_s() - t;
            break;
        }
        case RESULT: {
            int st = rd_i();
            double robj = rd_d(), riters = rd_d(), rsec = rd_d();
            nsolve[kind]++;
            trec[kind] += rsec;
            trep[kind] += tsolve;
            irec[kind] += riters;
            irep[kind] += iters;
            if (st == 1 && ok) {
                nboth[kind]++;
                double dev = fabs(robj - obj) / (1.0 + fabs(robj));
                if (dev > maxdev[kind]) maxdev[kind] = dev;
                if (dev > 1e-6) nmismatch[kind]++;
            } else if ((st == 1) != ok) nmismatch[kind]++;
            if (verbose)
                printf("%s %ld: rec st %d obj %.6f it %.0f %.4fs | %s ok %d obj %.6f it %.0f %.4fs\n", kname[kind],
                       nsolve[kind], st, robj, riters, rsec, B.name, ok, obj, iters, tsolve);
            solves++;
            if (max_solves > 0 && solves >= max_solves) goto done;
            break;
        }
        default:
            fprintf(stderr, "unknown op %d\n", op);
            return 2;
        }
    }
done:
    t_total = now_s() - t_total;
    printf("lpreplay backend %s: %ld solves, %.1f s total replay (solve calls %.1f s)\n", B.name, solves, t_total,
           trep[0] + trep[1] + trep[2]);
    for (int k = 0; k < 3; k++)
        if (nsolve[k])
            printf("  %-5s n %6ld  recorded %8.2f s (%.0f iters)  replayed %8.2f s (%.0f iters)  both optimal %ld  max obj dev %.1e  mismatches %ld\n",
                   kname[k], nsolve[k], trec[k], irec[k], trep[k], irep[k], nboth[k], maxdev[k], nmismatch[k]);
    return 0;
}
