#include "verify/verify.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cvrpbcp {

CheckResult check_solution(const Instance &I, const Solution &S, bool enforce_fleet)
{
    CheckResult R;
    std::vector<int> seen(I.nv(), 0);
    long long total = 0;
    int k = 0;
    for (const auto &r : S.routes) {
        if (r.empty()) continue;
        k++;
        long long load = 0, t = 0;
        int prev = 0;
        for (int c : r) {
            if (c < 1 || c > I.n) { R.error = "customer index out of range: " + std::to_string(c); return R; }
            if (seen[c]++) { R.error = "customer visited twice: " + std::to_string(c); return R; }
            load += I.demand[c];
            total += I.d(prev, c);
            if (I.has_tw) {
                t = std::max<long long>(t + I.d(prev, c), I.tw_lo[c]);
                if (t > I.tw_hi[c]) {
                    R.error = "time window violated at customer " + std::to_string(c) + " (arrival " +
                              std::to_string(t) + " > " + std::to_string(I.tw_hi[c]) + ")";
                    return R;
                }
                t += I.service[c];
            }
            prev = c;
        }
        if (load > I.capacity) {
            R.error = "capacity exceeded on route " + std::to_string(k) + " (" + std::to_string(load) + " > " +
                      std::to_string(I.capacity) + ")";
            return R;
        }
        total += I.d(prev, 0);
        if (I.has_tw) {
            t += I.d(prev, 0);
            if (t > I.tw_hi[0]) { R.error = "depot time window violated on route " + std::to_string(k); return R; }
        }
    }
    for (int c = 1; c <= I.n; c++)
        if (!seen[c]) { R.error = "customer not visited: " + std::to_string(c); return R; }
    if (enforce_fleet && I.vehicles > 0 && k > I.vehicles) {
        R.error = "too many vehicles: " + std::to_string(k) + " > " + std::to_string(I.vehicles);
        return R;
    }
    R.feasible = true;
    R.cost = total;
    R.vehicles = k;
    return R;
}

Solution read_sol(const std::string &path, int cost_scale)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    Solution S;
    std::string line;
    while (std::getline(f, line)) {
        if (line.compare(0, 5, "Route") == 0) {
            size_t c = line.find(':');
            if (c == std::string::npos) continue;
            std::istringstream is(line.substr(c + 1));
            std::vector<int> r;
            for (int v; is >> v;) r.push_back(v);
            S.routes.push_back(r);
        } else if (line.compare(0, 4, "Cost") == 0 || line.compare(0, 4, "cost") == 0) {
            double v = std::atof(line.c_str() + 4);
            S.cost = std::llround(v * cost_scale);
        }
    }
    return S;
}

std::string format_sol(const Solution &S, int cost_scale)
{
    std::ostringstream o;
    int k = 0;
    for (const auto &r : S.routes) {
        if (r.empty()) continue;
        o << "Route #" << ++k << ":";
        for (int c : r) o << ' ' << c;
        o << '\n';
    }
    if (S.cost >= 0) {
        if (cost_scale == 1) o << "Cost " << S.cost << '\n';
        else o << "Cost " << S.cost / cost_scale << '.' << std::abs(S.cost % cost_scale) << '\n';
    }
    return o.str();
}

void write_sol(const std::string &path, const Solution &S, int cost_scale)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << format_sol(S, cost_scale);
}

} // namespace cvrpbcp
