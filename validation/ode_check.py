"""Direct check of the closed-form P(E) used by independent_mpmath_ref.py.

The stationary Schroedinger equation (hbar = m = 1) is integrated through the same potential
V(x) = A y/(1+y) + B y/(1+y)^2, y = exp(x/L), from the transmitted side (pure outgoing wave e^{i k2 x})
to the incident side, where psi = c+ e^{i k1 x} + c- e^{-i k1 x} gives P = (k2/k1) / |c+|^2.
k1 = sqrt(2E), k2 = sqrt(2(E - A)), E from the left asymptote, as in the reference. Both branches of D
(cosh and cos) are covered. Needs numpy and scipy.

    python3 validation/ode_check.py
"""
import math
import os
import sys

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from independent_mpmath_ref import P_of_E, setup  # noqa: E402


def P_ode(E, A, B, L, X=45.0):
    """Transmission probability by integrating psi'' = 2 (V - E) psi from x = +X L to x = -X L."""
    k1, k2 = math.sqrt(2 * E), math.sqrt(2 * (E - A))

    def rhs(x, s):
        y = math.exp(x / L)
        f = 2 * (A * y / (1 + y) + B * y / (1 + y) ** 2 - E)
        return [s[2], s[3], f * s[0], f * s[1]]

    x0 = X * L
    psi0 = complex(math.cos(k2 * x0), math.sin(k2 * x0))
    dpsi0 = 1j * k2 * psi0
    sol = solve_ivp(rhs, [x0, -x0], [psi0.real, psi0.imag, dpsi0.real, dpsi0.imag],
                    method="DOP853", rtol=1e-13, atol=1e-15)
    psi = sol.y[0, -1] + 1j * sol.y[1, -1]
    dpsi = sol.y[2, -1] + 1j * sol.y[3, -1]
    cplus = 0.5 * (psi + dpsi / (1j * k1))
    return k2 / (k1 * abs(cplus) ** 2)


if __name__ == "__main__":
    mp.mp.dps = 30
    worst = 0.0
    # barriers v1, v2 and u = h nu / k_B T in k_B T units, as in the reference; t^2 = energy above the higher asymptote
    for v1, v2, u in ((4, 4, 2), (1, 3, 1.5), (10, 2, 4), (0.05, 0.05, 2), (0.3, 0.1, 5)):
        A, B, L, D = setup(v1, v2, u)
        branch = "cos" if 8 * B * L ** 2 < 1 else "cosh"
        E0 = max(mp.mpf(0), A)
        print(f"v1 = {v1}, v2 = {v2}, u = {u}: L = {mp.nstr(L, 6)}, {branch} branch")
        for t2 in (0.05, 0.3, 1.0, 3.0, 8.0):
            E = E0 + mp.mpf(t2) * max(v1, v2)
            p_ref = P_of_E(E, A, L, D)
            p_ode = P_ode(float(E), float(A), float(B), float(L))
            rel = abs(p_ode / float(p_ref) - 1)
            worst = max(worst, rel)
            print(f"   E - E0 = {float(E - E0):9.4f}   P closed form {float(p_ref):.13e}   ODE {p_ode:.13e}   rel. {rel:.1e}")
    print(f"largest relative difference: {worst:.1e}")
    assert worst < 1e-9, worst
    print("closed-form P(E) agrees with the direct integration")
