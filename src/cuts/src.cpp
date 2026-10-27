#include "cuts/src.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace cvrpbcp {

std::vector<SrcCut> separate_src(const Instance &I, const Master &M, const std::vector<long long> &existing,
                                 int max_cuts, double min_violation)
{
    const int n = I.n, nv = I.nv();
    const auto &lam = M.lambda();
    std::unordered_set<long long> skip(existing.begin(), existing.end());

    // fractional routes and their distinct customer sets
    struct FR { int r; double l; std::vector<int> cust; };
    std::vector<FR> frs;
    for (size_t r = 0; r < lam.size(); r++) {
        if (lam[r] < 1e-6) continue;
        FR f{(int)r, lam[r], M.route((int)r).seq};
        std::sort(f.cust.begin(), f.cust.end());
        f.cust.erase(std::unique(f.cust.begin(), f.cust.end()), f.cust.end());
        frs.push_back(std::move(f));
    }
    // pair flow pf(i,j) = sum of lambda over routes containing both;
    // triple flow tf(i,j,k) likewise
    std::vector<double> pf((size_t)nv * nv, 0.0);
    std::unordered_map<long long, double> tf;
    for (const FR &f : frs) {
        const auto &c = f.cust;
        for (size_t a = 0; a < c.size(); a++)
            for (size_t b = a + 1; b < c.size(); b++) {
                pf[(size_t)c[a] * nv + c[b]] += f.l;
                for (size_t d = b + 1; d < c.size(); d++) tf[src_key(nv, c[a], c[b], c[d])] += f.l;
            }
    }
    auto PF = [&](int i, int j) { return i < j ? pf[(size_t)i * nv + j] : pf[(size_t)j * nv + i]; };
    // candidates: every triple containing a pair with positive pair flow
    std::vector<std::pair<int, int>> pairs;
    for (int i = 1; i <= n; i++)
        for (int j = i + 1; j <= n; j++)
            if (pf[(size_t)i * nv + j] > 1e-6) pairs.push_back({i, j});
    std::unordered_set<long long> seen;
    std::vector<SrcCut> out;
    for (auto [i, j] : pairs) {
        double pij = PF(i, j);
        for (int k = 1; k <= n; k++) {
            if (k == i || k == j) continue;
            double v = pij + PF(i, k) + PF(j, k);
            if (v <= 1.0 + min_violation) continue;
            long long key = src_key(nv, i, j, k);
            auto it = tf.find(key);
            if (it != tf.end()) v -= 2.0 * it->second;
            if (v <= 1.0 + min_violation) continue;
            if (skip.count(key) || !seen.insert(key).second) continue;
            SrcCut cut;
            cut.c[0] = std::min({i, j, k});
            cut.c[2] = std::max({i, j, k});
            cut.c[1] = i + j + k - cut.c[0] - cut.c[2];
            cut.in_c.assign(nv, 0);
            cut.in_m.assign(nv, 0);
            for (int t = 0; t < 3; t++) cut.in_c[cut.c[t]] = cut.in_m[cut.c[t]] = 1;
            cut.violation = v - 1.0;
            out.push_back(std::move(cut));
        }
    }
    std::sort(out.begin(), out.end(), [](const SrcCut &a, const SrcCut &b) { return a.violation > b.violation; });
    if ((int)out.size() > max_cuts) out.resize(max_cuts);
    // memory: on every fractional route with >= 2 visits in C, the vertices
    // between its first and last visit to C
    for (SrcCut &cut : out) {
        for (const FR &f : frs) {
            const auto &seq = M.route(f.r).seq;
            int first = -1, last = -1, cnt = 0;
            for (size_t p = 0; p < seq.size(); p++)
                if (cut.in_c[seq[p]]) { cnt++; if (first < 0) first = (int)p; last = (int)p; }
            if (cnt < 2) continue;
            for (int p = first; p <= last; p++) cut.in_m[seq[p]] = 1;
        }
    }
    return out;
}

} // namespace cvrpbcp
