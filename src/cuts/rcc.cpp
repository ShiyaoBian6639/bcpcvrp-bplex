#include "cuts/rcc.h"

#include "capsep.h"
#include "cnstrmgr.h"

#include <algorithm>
#include <cmath>

namespace cvrpbcp {

#define NEWC ((CnstrMgrPointer)new_cuts_)
#define OLDC ((CnstrMgrPointer)old_cuts_)

RccSeparator::RccSeparator(const Instance &I) : I_(I)
{
    CnstrMgrPointer a = nullptr, b = nullptr;
    CMGR_CreateCMgr(&a, 200);
    CMGR_CreateCMgr(&b, 2000);
    new_cuts_ = a;
    old_cuts_ = b;
    demand_.assign(I.n + 2, 0.0);
    for (int i = 1; i <= I.n; i++) demand_[i] = I.demand[i];
}

RccSeparator::~RccSeparator()
{
    CnstrMgrPointer a = NEWC, b = OLDC;
    if (a) CMGR_FreeMemCMgr(&a);
    if (b) CMGR_FreeMemCMgr(&b);
}

std::vector<RccCut> RccSeparator::separate(const std::vector<double> &x,
                                           const std::vector<const std::vector<char> *> &existing, int max_cuts,
                                           double min_violation, bool &integer_feasible)
{
    const int n = I_.n, nv = I_.nv();
    // CVRPSEP's "existing cuts" = the RCC rows the master holds now (rebuilt
    // each call so that deleted cuts can be found again)
    {
        CnstrMgrPointer b = OLDC;
        CMGR_FreeMemCMgr(&b);
        CMGR_CreateCMgr(&b, std::max<int>(100, (int)existing.size() + 100));
        old_cuts_ = b;
        list_.resize(n + 2);
        for (const std::vector<char> *s : existing) {
            int size = 0;
            long long dem = 0;
            for (int v = 1; v <= n; v++)
                if ((*s)[v]) { list_[++size] = v; dem += I_.demand[v]; }
            if (size == 0) continue;
            int k = (int)((dem + I_.capacity - 1) / I_.capacity);
            CMGR_AddCnstr(OLDC, CMGR_CT_CAP, 0, size, list_.data(), (double)(size - k));
        }
    }
    // edges with positive flow, 1-based arrays, depot = n+1
    std::vector<int> tail(1, 0), head(1, 0);
    std::vector<double> ex(1, 0.0);
    for (int i = 0; i < nv; i++)
        for (int j = i + 1; j < nv; j++) {
            double v = x[(size_t)i * nv + j];
            if (v < 1e-9) continue;
            tail.push_back(i == 0 ? n + 1 : i);
            head.push_back(j);
            ex.push_back(v);
        }
    const int nedges = (int)tail.size() - 1;
    char integer_and_feasible = 0;
    double max_violation = 0;
    CAPSEP_SeparateCapCuts(n, demand_.data(), (double)I_.capacity, nedges, tail.data(), head.data(), ex.data(),
                           OLDC, max_cuts, 1e-6, min_violation, &integer_and_feasible, &max_violation, NEWC);
    integer_feasible = integer_and_feasible != 0;
    std::vector<RccCut> out;
    for (int c = 0; c < NEWC->Size; c++) {
        CnstrPointer cp = NEWC->CPL[c];
        if (cp->CType != CMGR_CT_CAP) continue;
        RccCut cut;
        cut.in_s.assign(nv, 0);
        int size = 0;
        for (int j = 1; j <= cp->IntListSize; j++) {
            int v = cp->IntList[j];
            if (v >= 1 && v <= n) { cut.in_s[v] = 1; size++; }
        }
        cut.k = size - (int)std::llround(cp->RHS);
        if (cut.k <= 0 || size == 0) continue;
        double cross = 0;
        for (int i = 0; i < nv; i++)
            for (int j = i + 1; j < nv; j++)
                if (cut.in_s[i] != cut.in_s[j]) cross += x[(size_t)i * nv + j];
        cut.violation = 2.0 * cut.k - cross;
        out.push_back(std::move(cut));
    }
    // clear the new-cut manager for the next call (the master's rows are the memory)
    for (int c = 0; c < NEWC->Size; c++) CMGR_MoveCnstr(NEWC, OLDC, c, 0);
    NEWC->Size = 0;
    std::sort(out.begin(), out.end(), [](const RccCut &a, const RccCut &b) { return a.violation > b.violation; });
    return out;
}

} // namespace cvrpbcp
