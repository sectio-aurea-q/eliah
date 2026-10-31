// ELIAH C++ interface: header-only, C++17, a thin layer over the C interface in eliah.h.
// Link with -leliah (target/release/libeliah.{so,dylib}). Errors become eliah::error.
#ifndef ELIAH_HPP
#define ELIAH_HPP

#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "eliah.h"

namespace eliah {

/// Result of a kappa evaluation; the fields are those of EliahResult.
struct Kappa {
    double value;         ///< kappa; +inf when it exceeds double range, use ln_value then
    double ln_value;      ///< ln kappa, finite even when value overflows
    double quad_err;      ///< quadrature error estimate (absolute); an estimate, not a bound
    double tail_bound;    ///< rigorous bound on the truncated tail (absolute)
    bool cos_branch;      ///< true if 4 sqrt(V1 V2) < h nu, where D = cos(2 pi |d|)
    std::uint64_t evals;  ///< number of integrand evaluations
};

/// Thrown on invalid input (code 1: non-finite or non-positive) and on failure
/// (code 2: no convergence or overflow).
class error : public std::runtime_error {
public:
    error(int code, const char* what) : std::runtime_error(what), code_(code) {}
    int code() const noexcept { return code_; }

private:
    int code_;
};

namespace detail {
inline Kappa checked(std::int32_t rc, const EliahResult& r) {
    if (rc == 1) throw error(1, "eliah: invalid input (non-finite or non-positive)");
    if (rc != 0) throw error(static_cast<int>(rc), "eliah: no convergence or overflow");
    return Kappa{r.value, r.ln_value, r.quad_err, r.tail_bound, r.cos_branch != 0, r.evals};
}
}  // namespace detail

/// Reduced units: barrier heights v1, v2 in k_B T (reactant and product side), u = h nu / k_B T.
inline Kappa kappa(double v1, double v2, double u) {
    EliahResult r;
    const std::int32_t rc = eliah_kappa(v1, v2, u, &r);
    return detail::checked(rc, r);
}

/// SI units: barriers in J/mol, |imaginary frequency| in cm^-1 (the sign is ignored), T in K.
/// Barriers are used as given; pass zero-point-corrected values for the ZPE-corrected model.
inline Kappa kappa_si(double v1_j_mol, double v2_j_mol, double nu_cm, double t_k) {
    EliahResult r;
    const std::int32_t rc = eliah_kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k, &r);
    return detail::checked(rc, r);
}

/// ln P(eps), eps = energy above the higher asymptote in k_B T.
inline double log_transmission(double v1, double v2, double u, double eps) {
    const double v = eliah_log_transmission(v1, v2, u, eps);
    if (std::isnan(v)) throw error(1, "eliah: invalid input");
    return v;
}

/// P(E), E in J/mol above the higher asymptote, barriers in J/mol, nu in cm^-1.
inline double transmission_si(double e_j_mol, double v1_j_mol, double v2_j_mol, double nu_cm) {
    const double p = eliah_transmission_si(e_j_mol, v1_j_mol, v2_j_mol, nu_cm);
    if (std::isnan(p)) throw error(1, "eliah: invalid input");
    return p;
}

}  // namespace eliah

#endif  // ELIAH_HPP
