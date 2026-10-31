// The C++ interface against the certified grid (conformance/eckart_certified.csv) and its errors.
//   cargo build --release
//   c++ -std=c++17 -Wall -Wextra -Wpedantic -Iinclude tests/test_cpp.cpp -Ltarget/release -leliah -o test_cpp
//   LD_LIBRARY_PATH=target/release ./test_cpp conformance/eckart_certified.csv
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "eliah.hpp"

int main(int argc, char** argv) {
    const char* path = argc > 1 ? argv[1] : "conformance/eckart_certified.csv";
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path);
        return 2;
    }
    std::string line;
    auto next = [&in, &line]() {
        if (!std::getline(in, line)) return false;
        if (!line.empty() && line.back() == '\r') line.pop_back();  // CRLF files
        return true;
    };
    next();
    std::map<std::string, std::size_t> col;
    {
        std::stringstream ss(line);
        std::string name;
        for (std::size_t i = 0; std::getline(ss, name, ','); ++i) col[name] = i;
    }
    for (const char* need : {"v1", "v2", "ustar", "kappa_mid", "cos_branch"}) {
        if (!col.count(need)) {
            std::fprintf(stderr, "column %s missing\n", need);
            return 2;
        }
    }
    const double tol = 1e-12;
    std::size_t n = 0, bad = 0, cos_rows = 0, cos_flag_wrong = 0;
    double worst = 0.0;
    while (next()) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) f.push_back(cell);
        const double v1 = std::strtod(f[col["v1"]].c_str(), nullptr);
        const double v2 = std::strtod(f[col["v2"]].c_str(), nullptr);
        const double u = std::strtod(f[col["ustar"]].c_str(), nullptr);
        const double ref = std::strtod(f[col["kappa_mid"]].c_str(), nullptr);
        const bool cos_ref = f[col["cos_branch"]] == "1";
        ++n;
        try {
            const eliah::Kappa k = eliah::kappa(v1, v2, u);
            const double rel = std::isfinite(k.value) ? std::fabs(k.value / ref - 1.0) : INFINITY;
            if (rel > worst) worst = rel;
            if (rel > tol) ++bad;
            if (cos_ref) ++cos_rows;
            if (k.cos_branch != cos_ref) ++cos_flag_wrong;
        } catch (const eliah::error& e) {
            ++bad;
            std::printf("exception at v1=%g v2=%g u=%g: %s\n", v1, v2, u, e.what());
        }
    }
    std::printf("%zu certified points: %zu within %.0e, worst relative error %.2e; cos branch %zu rows, "
                "flag mismatches %zu\n", n, n - bad, tol, worst, cos_rows, cos_flag_wrong);

    int err_ok = 0;
    try { eliah::kappa(0.0, 1.0, 1.0); } catch (const eliah::error& e) { err_ok += e.code() == 1; }
    try { eliah::kappa_si(NAN, 1e5, 1000.0, 300.0); } catch (const eliah::error& e) { err_ok += e.code() == 1; }
    try { eliah::log_transmission(1.0, 1.0, -1.0, 0.5); } catch (const eliah::error& e) { err_ok += e.code() == 1; }
    const double p = eliah::transmission_si(100e3, 160e3, 160e3, 1200.0);
    const bool p_ok = p > 0.0 && p < 1.0;
    std::printf("errors raised for invalid input: %d of 3; P(E) in (0, 1): %s\n", err_ok, p_ok ? "yes" : "no");

    const bool ok = n == 1197 && bad == 0 && cos_flag_wrong == 0 && err_ok == 3 && p_ok;
    std::printf("%s\n", ok ? "C++ interface: PASS" : "C++ interface: FAIL");
    return ok ? 0 : 1;
}
