// ELIAH from C++.
//   cargo build --release
//   c++ -std=c++17 -Iinclude examples/cpp_example.cpp -Ltarget/release -leliah -o cpp_example
//   LD_LIBRARY_PATH=target/release ./cpp_example        (macOS: DYLD_LIBRARY_PATH)
#include <cstdio>

#include "eliah.hpp"

int main() {
    // 160 kJ/mol from both sides, imaginary frequency 1200 cm^-1 (sign ignored), 200 K
    const eliah::Kappa k = eliah::kappa_si(160e3, 160e3, -1200.0, 200.0);
    std::printf("kappa = %.10g  (ln %.10g, quad err %.1e, %s branch, %llu evaluations)\n", k.value,
                k.ln_value, k.quad_err, k.cos_branch ? "cos" : "cosh",
                static_cast<unsigned long long>(k.evals));

    try {
        eliah::kappa(0.0, 1.0, 1.0);  // a barrier height must be positive
    } catch (const eliah::error& e) {
        std::printf("invalid input: code %d, %s\n", e.code(), e.what());
    }

    std::printf("P(E = 100 kJ/mol) = %.6e\n", eliah::transmission_si(100e3, 160e3, 160e3, 1200.0));
    return 0;
}
