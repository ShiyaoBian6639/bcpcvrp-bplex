// LP trace: records every operation the master performs on its LP solver
// (load, column/row additions and deletions, bound changes, parameter
// changes, solves) in a compact binary file, so that the identical sequence
// can be replayed on another solver (tools/lpreplay.c: bplex or CPLEX).
// Enabled by the environment variable CVRPBCP_LPTRACE=<file>.
//
// Format (little endian): "CVRPLPT1" magic, then records
//   op (u8) followed by its payload; arrays are a count (i32) and the items.
//   LOAD      ncols nrows obj[ncols] lo[ncols] up[ncols] rlo[nrows] rup[nrows]
//             colstart[ncols+1] rowidx[nnz] val[nnz]
//   ADDCOLS   n obj[n] lo[n] up[n] colstart[n+1] rowidx[nnz] val[nnz]
//   DELCOLS   n idx[n]
//   ADDROWS   n lo[n] up[n] rowstart[n+1] colidx[nnz] val[nnz]
//   DELROWS   n idx[n]
//   CHGROW    i lo up
//   CHGCOL    j lo up
//   CHGCOLS   n idx[n] lo[n] up[n]
//   CHGTYPE   j is_integer
//   PARAM     name (i32 length + bytes) value (f64)
//   SOLVE     kind (0 LP, 1 MIP, 2 LP probe with iteration limit)
//   SLACKBASIS (no payload: restart from the all-slack basis)
//   RESULT    status objval iterations seconds   (what the recording solver did)
// Infinite bounds are written as +-1e30.
#pragma once

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

namespace cvrpbcp {

class LpTrace {
  public:
    enum Op : uint8_t { LOAD = 1, ADDCOLS, DELCOLS, ADDROWS, DELROWS, CHGROW, CHGCOL, CHGCOLS, CHGTYPE, PARAM, SOLVE, RESULT, SLACKBASIS };
    explicit LpTrace(const char *path) : f_(std::fopen(path, "wb"))
    {
        if (f_) std::fwrite("CVRPLPT1", 1, 8, f_);
    }
    ~LpTrace()
    {
        if (f_) std::fclose(f_);
    }
    LpTrace(const LpTrace &) = delete;
    LpTrace &operator=(const LpTrace &) = delete;
    bool on() const { return f_ != nullptr; }

    void load(int ncols, int nrows, const double *obj, const double *lo, const double *up, const double *rlo,
              const double *rup, const int *cs, const int *ri, const double *v)
    {
        op(LOAD);
        i(ncols);
        i(nrows);
        arr(obj, ncols);
        arr(lo, ncols);
        arr(up, ncols);
        arr(rlo, nrows);
        arr(rup, nrows);
        arr(cs, ncols + 1);
        arr(ri, cs[ncols]);
        arr(v, cs[ncols]);
    }
    // lo/up may be null (0 / +inf)
    void add_cols(int n, const double *obj, const double *lo, const double *up, const int *cs, const int *ri,
                  const double *v)
    {
        op(ADDCOLS);
        i(n);
        arr(obj, n);
        i(n);
        for (int j = 0; j < n; j++) d(lo ? lo[j] : 0.0);
        i(n);
        for (int j = 0; j < n; j++) d(up ? up[j] : 1e30);
        arr(cs, n + 1);
        arr(ri, cs[n]);
        arr(v, cs[n]);
    }
    void del_cols(int n, const int *idx)
    {
        op(DELCOLS);
        i(n);
        arr(idx, n);
    }
    void add_rows(int n, const double *lo, const double *up, const int *rs, const int *ci, const double *v)
    {
        op(ADDROWS);
        i(n);
        arr(lo, n);
        arr(up, n);
        arr(rs, n + 1);
        arr(ci, rs[n]);
        arr(v, rs[n]);
    }
    void del_rows(int n, const int *idx)
    {
        op(DELROWS);
        i(n);
        arr(idx, n);
    }
    void chg_row(int row, double lo, double up)
    {
        op(CHGROW);
        i(row);
        d(lo);
        d(up);
    }
    void chg_col(int col, double lo, double up)
    {
        op(CHGCOL);
        i(col);
        d(lo);
        d(up);
    }
    void chg_cols(int n, const int *idx, const double *lo, const double *up)
    {
        op(CHGCOLS);
        i(n);
        arr(idx, n);
        arr(lo, n);
        arr(up, n);
    }
    void chg_type(int col, int is_integer)
    {
        op(CHGTYPE);
        i(col);
        i(is_integer);
    }
    void param(const char *name, double value)
    {
        op(PARAM);
        int len = (int)std::strlen(name);
        i(len);
        std::fwrite(name, 1, len, f_);
        d(value);
    }
    void slack_basis() { op(SLACKBASIS); }
    void solve(int kind)
    {
        op(SOLVE);
        i(kind);
    }
    void result(int status, double objval, double iterations, double seconds)
    {
        op(RESULT);
        i(status);
        d(objval);
        d(iterations);
        d(seconds);
        std::fflush(f_);
    }

  private:
    void op(Op o) { std::fputc(o, f_); }
    void i(int x) { std::fwrite(&x, sizeof x, 1, f_); }
    void d(double x) { std::fwrite(&x, sizeof x, 1, f_); }
    void arr(const int *p, int n)
    {
        i(n);
        if (n > 0) std::fwrite(p, sizeof(int), (size_t)n, f_);
    }
    void arr(const double *p, int n)
    {
        i(n);
        if (n > 0) std::fwrite(p, sizeof(double), (size_t)n, f_);
    }
    std::FILE *f_;
};

} // namespace cvrpbcp
