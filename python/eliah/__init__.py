"""ELIAH in pure Python: accurate Eckart tunnelling coefficients, no compiler needed.

    pip install "git+https://github.com/sectio-aurea-q/eliah"

    import eliah
    eliah.kappa(v1, v2, u)                    # reduced units: barriers in k_B T, u = h nu / k_B T
    eliah.kappa_si(V1, V2, nu_cm, T)          # J/mol, J/mol, cm^-1 (sign ignored), K
    eliah.kappa_full(v1, v2, u)               # dict: value, ln_value, quad_err, tail_bound, cos_branch, evals
    eliah.log_transmission(v1, v2, u, eps)    # ln P(eps), eps above the higher asymptote in k_B T
    eliah.transmission_si(E, V1, V2, nu_cm)   # P(E), E from the higher asymptote, all in J/mol

This is a line-by-line port of the Rust library (src/lib.rs) and gives the same results:
it passes the certified conformance grid (1197 points) within 1e-12, like the Rust version.
Conventions: Johnston & Heicklen (1962), Brown (1981); see the README.

Cite: S. Hackfort, ELIAH (formerly NETI), doi:10.5281/zenodo.23485066.
"""
import math

__version__ = "0.3.0"
__all__ = ["kappa", "kappa_si", "kappa_full", "log_transmission", "transmission_si",
           "Barrier", "EckartError", "REL_TOL", "TAIL", "MAX_EVALS"]

PI = math.pi
LN_2 = math.log(2.0)
EPS = 2.220446049250313e-16          # f64::EPSILON
MIN_POSITIVE = 2.2250738585072014e-308

# exact SI constants (CODATA 2018 = 2022)
H = 6.62607015e-34
C = 299792458.0
KB = 1.380649e-23
NA = 6.02214076e23
R = KB * NA

REL_TOL = 1e-13      # relative tolerance requested from the quadrature
TAIL = 80.0          # truncation distance above v_min (in k_B T)
MAX_EVALS = 200_000  # maximum integrand evaluations per kappa


class EckartError(ValueError):
    pass


def _check(x, name):
    if not math.isfinite(x):
        raise EckartError(f"{name} must be finite")
    if not x > 0.0:
        raise EckartError(f"{name} must be > 0")


def _log_sinh(x):
    return x - LN_2 + math.log1p(-math.exp(-2.0 * x)) if x > 1.0 else math.log(math.sinh(x))


def _log_cosh(x):
    return x - LN_2 + math.log1p(math.exp(-2.0 * x)) if x > 1.0 else math.log(math.cosh(x))


def _log_add_exp(x, y):
    if x > y:
        return x + math.log1p(math.exp(y - x))
    return y + math.log1p(math.exp(x - y))


def _log_sinh_or_neginf(x):
    return _log_sinh(x) if x > 0.0 else -math.inf


class Barrier:
    """Precomputed barrier geometry (reduced units)."""

    def __init__(self, v1, v2, u):
        _check(v1, "v1"); _check(v2, "v2"); _check(u, "u")
        lo, hi = (v1, v2) if v1 <= v2 else (v2, v1)
        k = 2.0 * PI / u
        a_lo, a_hi = k * lo, k * hi
        q = a_lo * a_hi - PI * PI / 4.0
        if q >= 0.0:
            self.log_d, self.log_cos2 = _log_cosh(2.0 * math.sqrt(q)), math.nan
        else:
            # cos(pi|d|) = sin(pi/2 - pi|d|), pi/2 - sqrt(pi^2/4 - p) = p / (pi/2 + sqrt(pi^2/4 - p))
            p = a_lo * a_hi
            w = p / (PI / 2.0 + math.sqrt(PI * PI / 4.0 - p))
            self.log_d, self.log_cos2 = math.nan, 2.0 * math.log(math.sin(w))
        self.v_min = lo
        self.k = k
        self.d_alpha = k * (hi - lo)
        self.inv_s = 1.0 / (1.0 / math.sqrt(a_lo) + 1.0 / math.sqrt(a_hi))
        self.q = q

    @property
    def cos_branch(self):
        return self.q < 0.0

    def _log_cosh_plus_d(self, x):
        if self.q >= 0.0:
            return _log_add_exp(_log_cosh(x), self.log_d)
        # cosh(x) + cos(2 pi |d|) = 2 sinh^2(x/2) + 2 cos^2(pi |d|)
        return LN_2 + _log_add_exp(2.0 * _log_sinh_or_neginf(0.5 * x), self.log_cos2)

    def log_p_and_q(self, eps):
        """(ln P, ln(1 - P)) at eps (energy above the higher asymptote, in k_B T)."""
        if math.isnan(eps):
            return math.nan, math.nan
        if eps <= 0.0:
            return -math.inf, 0.0
        a = 2.0 * math.sqrt(self.k * eps) * self.inv_s
        b = 2.0 * math.sqrt(self.d_alpha + self.k * eps) * self.inv_s
        log_den = self._log_cosh_plus_d(a + b)
        log_p = LN_2 + _log_sinh(a) + _log_sinh(b) - log_den
        if log_p < -LN_2:
            return log_p, math.log1p(-math.exp(log_p))
        log_q = self._log_cosh_plus_d(abs(a - b)) - log_den
        if log_q > 0.0:
            log_q = 0.0
        return math.log(-math.expm1(log_q)), log_q

    def log_transmission(self, eps):
        return self.log_p_and_q(eps)[0]

    def transmission(self, eps):
        return math.exp(self.log_transmission(eps))

    def _log_integrand_t(self, t):
        if t <= 0.0:
            return -math.inf
        e = t * t
        return math.log(2.0 * t) + self.log_transmission(e) + self.v_min - e


# Gauss-Kronrod 7/15 nodes and weights (QUADPACK qk15)
_XGK = (0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
        0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
        0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
        0.207784955007898467600689403773245, 0.0)
_WGK = (0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
        0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
        0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
        0.204432940075298892414161999234649, 0.209482141084727828012999174891714)
_WG = (0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
       0.381830050505118944950369775488975, 0.417959183673469387755102040816327)


def _gk15(f, a, b):
    c = 0.5 * (a + b)
    h = 0.5 * (b - a)
    fv1 = [0.0] * 7
    fv2 = [0.0] * 7
    fc = f(c)
    resk = fc * _WGK[7]
    resg = fc * _WG[3]
    resabs = abs(fc) * _WGK[7]
    for j in range(7):
        x = h * _XGK[j]
        f1, f2 = f(c - x), f(c + x)
        fv1[j], fv2[j] = f1, f2
        resk += _WGK[j] * (f1 + f2)
        resabs += _WGK[j] * (abs(f1) + abs(f2))
        if j % 2 == 1:
            resg += _WG[j // 2] * (f1 + f2)
    reskh = resk * 0.5
    resasc = _WGK[7] * abs(fc - reskh)
    for j in range(7):
        resasc += _WGK[j] * (abs(fv1[j] - reskh) + abs(fv2[j] - reskh))
    result, resabs, resasc = resk * h, resabs * abs(h), resasc * abs(h)
    err = abs((resk - resg) * h)
    if resasc != 0.0 and err != 0.0:
        err = resasc * min(1.0, (200.0 * err / resasc) ** 1.5)
    floor = 50.0 * EPS * resabs
    if resabs > MIN_POSITIVE / (50.0 * EPS):
        err = max(err, floor)
    return result, err, resabs


def _qag(f, breaks, rel_tol):
    panels = []
    evals = 0
    for a, b in zip(breaks[:-1], breaks[1:]):
        v, e, _ = _gk15(f, a, b)
        evals += 15
        panels.append([a, b, v, e])
    while True:
        total = 0.0
        err = 0.0
        for p in panels:
            total += p[2]
        for p in panels:
            err += p[3]
        if not (math.isfinite(total) and math.isfinite(err)):
            raise EckartError("quadrature did not converge within the evaluation limit")
        floor = 64.0 * EPS * abs(total)
        if err <= max(rel_tol * abs(total), floor):
            return total, err, evals
        if evals >= MAX_EVALS:
            if err <= 1e3 * floor:
                return total, err, evals
            raise EckartError("quadrature did not converge within the evaluation limit")
        i, best = 0, -1.0
        for j, p in enumerate(panels):
            if p[3] > best:
                i, best = j, p[3]
        # swap_remove, as in the Rust version (keeps the summation order identical)
        p = panels[i]
        panels[i] = panels[-1]
        panels.pop()
        a, b = p[0], p[1]
        m = 0.5 * (a + b)
        if not (a < m < b):
            panels.append([a, b, p[2], 0.0])
            continue
        v1, e1, _ = _gk15(f, a, m)
        v2, e2, _ = _gk15(f, m, b)
        evals += 30
        panels.append([a, m, v1, e1])
        panels.append([m, b, v2, e2])


def kappa_full(v1, v2, u):
    """kappa in reduced units, with diagnostics (dict)."""
    bar = Barrier(v1, v2, u)
    vm = bar.v_min
    eps_max = vm + TAIL
    cand = [0.25 * vm, 0.5 * vm, 0.75 * vm, 0.9 * vm, vm, 1.1 * vm,
            vm + 2.0, vm + 8.0, vm + 25.0, 0.5, 2.0, 5.0, 12.0, 30.0]
    w = u / (2.0 * PI)
    for n in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
        cand.append(vm - n * w)
        cand.append(vm + n * w)
    cand = sorted(x for x in cand if 0.0 < x < eps_max)
    br = [0.0]
    for x in cand:
        if x > br[-1] * (1.0 + 1e-9):
            br.append(x)
    br.append(eps_max)
    ts = [math.sqrt(e) for e in br]
    c = -math.inf
    for t0, t1 in zip(ts[:-1], ts[1:]):
        for j in range(1, 8):
            c = max(c, bar._log_integrand_t(t0 + (t1 - t0) * j / 8.0))
    if not math.isfinite(c):
        raise EckartError("quadrature did not converge within the evaluation limit")
    val, err, evals = _qag(lambda t: math.exp(bar._log_integrand_t(t) - c), ts, REL_TOL)
    ln_value = c + math.log(val)
    try:
        value = math.exp(ln_value)
    except OverflowError:
        value = math.inf
    try:
        scale = math.exp(c)
        quad_err = err * scale
    except OverflowError:
        quad_err = math.inf
    return {"value": value, "ln_value": ln_value, "quad_err": quad_err,
            "tail_bound": math.exp(vm - eps_max), "cos_branch": bar.cos_branch, "evals": evals}


def kappa(v1, v2, u):
    """Eckart kappa; barriers v1, v2 in k_B T (reactant and product side), u = h nu / k_B T."""
    return kappa_full(v1, v2, u)["value"]


def kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k, full=False):
    """Eckart kappa from barriers in J/mol, imaginary frequency in cm^-1 (sign ignored), T in K.
    Barriers are used as given (pass ZPE-corrected barriers for the usual ZPE-corrected model)."""
    _check(t_k, "T"); _check(v1_j_mol, "V1"); _check(v2_j_mol, "V2")
    if not math.isfinite(nu_cm):
        raise EckartError("nu must be finite")
    nu = abs(nu_cm)
    if not nu > 0.0:
        raise EckartError("|nu| must be > 0")
    rt = R * t_k
    u = H * C * 100.0 * nu / (KB * t_k)
    r = kappa_full(v1_j_mol / rt, v2_j_mol / rt, u)
    return r if full else r["value"]


def log_transmission(v1, v2, u, eps):
    """ln P(eps) in reduced units (eps above the higher asymptote, in k_B T)."""
    return Barrier(v1, v2, u).log_transmission(eps)


def transmission_si(e_j_mol, v1_j_mol, v2_j_mol, nu_cm):
    """Microcanonical P(E): E from the higher asymptote, barriers in J/mol, nu in cm^-1."""
    _check(v1_j_mol, "V1"); _check(v2_j_mol, "V2")
    nu = abs(nu_cm)
    _check(nu, "|nu|")
    hnu = H * C * 100.0 * nu * NA
    return Barrier(v1_j_mol / hnu, v2_j_mol / hnu, 1.0).transmission(e_j_mol / hnu)
