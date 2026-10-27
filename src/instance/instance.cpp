#include "instance/instance.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cvrpbcp {

namespace {

std::string trim(const std::string &s)
{
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string upper(std::string s)
{
    for (auto &ch : s) ch = (char)toupper((unsigned char)ch);
    return s;
}

// "KEY : value" -> key upper-cased, value trimmed; false if not of that shape
bool keyval(const std::string &line, std::string &key, std::string &val)
{
    size_t c = line.find(':');
    if (c == std::string::npos) return false;
    key = upper(trim(line.substr(0, c)));
    val = trim(line.substr(c + 1));
    return !key.empty() && key.find(' ') == std::string::npos;
}

bool is_section(const std::string &t)
{
    return t.size() >= 8 && t.compare(t.size() - 7, 7, "SECTION") == 0;
}

} // namespace

int Instance::min_vehicles() const
{
    long long tot = 0;
    for (int i = 1; i <= n; i++) tot += demand[i];
    return (int)((tot + capacity - 1) / capacity);
}

int Instance::ffd_vehicles() const
{
    std::vector<int> dem(demand.begin() + 1, demand.end());
    std::sort(dem.rbegin(), dem.rend());
    std::vector<int> bins;
    for (int q : dem) {
        bool placed = false;
        for (auto &b : bins)
            if (b + q <= capacity) { b += q; placed = true; break; }
        if (!placed) bins.push_back(q);
    }
    return (int)bins.size();
}

std::string Instance::summary() const
{
    std::ostringstream o;
    o << name << ": " << (has_tw ? "CVRPTW" : "CVRP") << " n=" << n << " Q=" << capacity
      << " K=" << vehicles << " scale=" << scale << " kmin=" << min_vehicles()
      << " ffd=" << ffd_vehicles();
    return o.str();
}

long long comment_optimum(const std::string &path)
{
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        std::string k, v;
        if (!keyval(line, k, v)) continue;
        if (k == "COMMENT") {
            size_t p = v.find("Optimal value");
            if (p == std::string::npos) return -1;
            p = v.find_first_of("0123456789", p);
            if (p == std::string::npos) return -1;
            return std::atoll(v.c_str() + p);
        }
        if (is_section(k)) break;
    }
    return -1;
}

Instance read_instance(const std::string &path, const Rounding *rounding_override)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    Instance I;
    std::string type, ewt = "EUC_2D", ewf;
    int dim = 0, depot_id = -1;
    double uniform_service = 0;
    bool have_uniform_service = false;
    std::vector<double> ex, ey, service_raw, twlo_raw, twhi_raw, ew;
    std::vector<int> edem;
    std::string line, section;
    std::vector<std::string> tokens;

    auto fail = [&](const std::string &m) { throw std::runtime_error(path + ": " + m); };

    while (std::getline(f, line)) {
        std::string t = trim(line);
        if (t.empty()) continue;
        if (upper(t) == "EOF") break;
        std::string k, v;
        if (keyval(t, k, v) && !is_section(k)) {
            section.clear();
            if (k == "NAME") I.name = v;
            else if (k == "COMMENT") {
                // CVRPLIB: "(Augerat et al, No of trucks: 5, Optimal value: 784)" and
                // "(Christophides and Eilon, Min no of trucks: 3, ...)" fix the fleet
                // (E-n30-k3's 534 needs exactly 3 trucks; 4 trucks give 503)
                size_t p = upper(v).find("NO OF TRUCKS");
                if (p != std::string::npos) {
                    p = v.find_first_of("0123456789", p);
                    if (p != std::string::npos) I.vehicles = std::atoi(v.c_str() + p);
                }
            }
            else if (k == "TYPE") type = upper(v);
            else if (k == "DIMENSION") dim = std::stoi(v);
            else if (k == "CAPACITY") I.capacity = (int)std::llround(std::stod(v));
            else if (k == "VEHICLES") I.vehicles = std::stoi(v);
            else if (k == "EDGE_WEIGHT_TYPE") ewt = upper(v);
            else if (k == "EDGE_WEIGHT_FORMAT") ewf = upper(v);
            else if (k == "SERVICE_TIME") { uniform_service = std::stod(v); have_uniform_service = true; }
            continue; // COMMENT, DISPLAY_DATA_TYPE, ... ignored
        }
        std::string up = upper(t);
        if (is_section(up)) {
            section = up;
            if (dim <= 0) fail("DIMENSION missing before " + section);
            if (section == "NODE_COORD_SECTION") { ex.assign(dim + 1, 0); ey.assign(dim + 1, 0); }
            if (section == "DEMAND_SECTION") edem.assign(dim + 1, 0);
            if (section == "SERVICE_TIME_SECTION") service_raw.assign(dim + 1, 0);
            if (section == "TIME_WINDOW_SECTION") { twlo_raw.assign(dim + 1, 0); twhi_raw.assign(dim + 1, 0); }
            continue;
        }
        std::istringstream is(t);
        tokens.clear();
        for (std::string w; is >> w;) tokens.push_back(w);
        auto id_of = [&](size_t need) {
            if (tokens.size() < need) fail("short line in " + section + ": " + t);
            int id = std::stoi(tokens[0]);
            if (id < 1 || id > dim) fail("id out of range in " + section + ": " + t);
            return id;
        };
        if (section == "NODE_COORD_SECTION") {
            int id = id_of(3);
            ex[id] = std::stod(tokens[1]);
            ey[id] = std::stod(tokens[2]);
        } else if (section == "DEMAND_SECTION") {
            int id = id_of(2);
            edem[id] = (int)std::llround(std::stod(tokens[1]));
        } else if (section == "SERVICE_TIME_SECTION") {
            int id = id_of(2);
            service_raw[id] = std::stod(tokens[1]);
        } else if (section == "TIME_WINDOW_SECTION") {
            int id = id_of(3);
            twlo_raw[id] = std::stod(tokens[1]);
            twhi_raw[id] = std::stod(tokens[2]);
        } else if (section == "DEPOT_SECTION") {
            int id = std::stoi(tokens[0]);
            if (id != -1 && depot_id == -1) depot_id = id;
        } else if (section == "EDGE_WEIGHT_SECTION") {
            for (auto &w : tokens) ew.push_back(std::stod(w));
        }
        // other sections ignored
    }
    if (dim <= 1) fail("no or bad DIMENSION");
    if (I.capacity <= 0) fail("no CAPACITY");
    if (depot_id == -1) depot_id = 1;
    I.has_tw = (type == "CVRPTW" || type == "VRPTW" || !twlo_raw.empty());
    I.n = dim - 1;
    I.rounding = I.has_tw ? Rounding::TruncOneDecimal : Rounding::Nearest;
    if (ewt == "EXPLICIT") I.rounding = Rounding::Explicit;
    if (rounding_override) I.rounding = *rounding_override;
    I.scale = (I.rounding == Rounding::TruncOneDecimal) ? 10 : 1;

    // map file ids (1..dim) to vertices: depot -> 0, others in file order
    std::vector<int> vert(dim + 1, -1);
    vert[depot_id] = 0;
    for (int id = 1, v = 1; id <= dim; id++)
        if (id != depot_id) vert[id] = v++;

    const int nv = dim;
    I.x.assign(nv, 0); I.y.assign(nv, 0); I.demand.assign(nv, 0);
    I.service.assign(nv, 0); I.tw_lo.assign(nv, 0); I.tw_hi.assign(nv, 0);
    for (int id = 1; id <= dim; id++) {
        int v = vert[id];
        if (!ex.empty()) { I.x[v] = ex[id]; I.y[v] = ey[id]; }
        if (!edem.empty()) I.demand[v] = edem[id];
        if (I.has_tw) {
            double s = 0;
            if (!service_raw.empty()) s = service_raw[id];
            else if (have_uniform_service && v != 0) s = uniform_service;
            I.service[v] = (int)std::llround(s * I.scale);
            I.tw_lo[v] = (int)std::llround(twlo_raw[id] * I.scale);
            I.tw_hi[v] = (int)std::llround(twhi_raw[id] * I.scale);
        }
    }
    I.demand[0] = 0;
    I.service[0] = 0;

    I.dist.assign((size_t)nv * nv, 0);
    if (ewt == "EXPLICIT") {
        std::vector<std::vector<double>> M(dim, std::vector<double>(dim, 0));
        size_t p = 0;
        auto get = [&]() {
            if (p >= ew.size()) fail("EDGE_WEIGHT_SECTION too short");
            return ew[p++];
        };
        if (ewf == "FULL_MATRIX") {
            for (int i = 0; i < dim; i++) for (int j = 0; j < dim; j++) M[i][j] = get();
        } else if (ewf == "LOWER_ROW") {
            for (int i = 1; i < dim; i++) for (int j = 0; j < i; j++) M[i][j] = M[j][i] = get();
        } else if (ewf == "LOWER_DIAG_ROW") {
            for (int i = 0; i < dim; i++) for (int j = 0; j <= i; j++) M[i][j] = M[j][i] = get();
        } else if (ewf == "UPPER_ROW") {
            for (int i = 0; i < dim; i++) for (int j = i + 1; j < dim; j++) M[i][j] = M[j][i] = get();
        } else if (ewf == "UPPER_DIAG_ROW") {
            for (int i = 0; i < dim; i++) for (int j = i; j < dim; j++) M[i][j] = M[j][i] = get();
        } else fail("unsupported EDGE_WEIGHT_FORMAT " + ewf);
        for (int a = 1; a <= dim; a++)
            for (int b = 1; b <= dim; b++)
                I.dist[(size_t)vert[a] * nv + vert[b]] = (int)std::llround(M[a - 1][b - 1]);
    } else if (ewt == "EUC_2D") {
        if (ex.empty()) fail("EUC_2D without NODE_COORD_SECTION");
        for (int i = 0; i < nv; i++)
            for (int j = 0; j < nv; j++) {
                double dx = I.x[i] - I.x[j], dy = I.y[i] - I.y[j];
                double e = std::sqrt(dx * dx + dy * dy);
                int c = (I.rounding == Rounding::TruncOneDecimal) ? (int)std::floor(e * 10.0 + 1e-9)
                                                                   : (int)std::llround(e);
                I.dist[(size_t)i * nv + j] = c;
            }
    } else fail("unsupported EDGE_WEIGHT_TYPE " + ewt);
    for (int i = 0; i < nv; i++) I.dist[(size_t)i * nv + i] = 0;

    if (I.name.empty()) {
        size_t s = path.find_last_of("/\\");
        I.name = path.substr(s == std::string::npos ? 0 : s + 1);
        size_t e = I.name.rfind('.');
        if (e != std::string::npos) I.name = I.name.substr(0, e);
    }
    return I;
}

} // namespace cvrpbcp
