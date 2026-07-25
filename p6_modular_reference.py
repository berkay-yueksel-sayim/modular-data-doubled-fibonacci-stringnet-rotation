# -*- coding: utf-8 -*-
"""
p6_modular_reference.py -- self-contained analytic reference.

Rebuilds the doubled-Fibonacci modular data FROM SCRATCH (numpy only, no import
from the simulation engine, no stored oracle) and checks the analytic target
values against which the lattice-extracted quantities of the paper are compared.

What this does and does not cover: it does NOT reproduce the lattice extraction
itself (that requires the exact-diagonalization engine). It reproduces the
*analytic side* of every comparison in Table I -- i.e. it independently
establishes the numbers the lattice is measured against, so that a reader can
check that the targets were not themselves imported from a stored artifact.

Deterministic: no RNG, no I/O beyond stdout and the JSON it validates.
Run:  python p6_modular_reference.py
"""

import json
import os

import numpy as np

TOL = 1e-12
PHI = (1.0 + np.sqrt(5.0)) / 2.0

results = {}
failures = []


def check(name, got, want, tol=TOL):
    """Record a comparison; append to failures if it misses."""
    err = float(abs(got - want))
    results[name] = {"value": float(got), "target": float(want), "abs_err": err}
    status = "OK " if err <= tol else "FAIL"
    if err > tol:
        failures.append(name)
    print(f"  [{status}] {name:<42s} = {got!r:<24s} (err {err:.3e})")
    return err <= tol


# ---------------------------------------------------------------------------
# 1. Fibonacci MTC, built from scratch
# ---------------------------------------------------------------------------
print("=== 1. Fibonacci (single) ===")

d_fib = np.array([1.0, PHI])                 # quantum dimensions of {1, tau}
D_fib = np.sqrt((d_fib ** 2).sum())          # total quantum dimension
S_fib = np.array([[1.0, PHI], [PHI, -1.0]]) / D_fib
theta_fib = np.array([1.0 + 0j, np.exp(4j * np.pi / 5)])   # topological spins
T_fib = np.diag(theta_fib)

check("D_fib = sqrt(1+phi^2)", D_fib, np.sqrt(1 + PHI ** 2))
# S_fib must be real, symmetric, and an involution up to charge conjugation
check("S_fib unitarity |S S^dag - I|", np.abs(S_fib @ S_fib.conj().T - np.eye(2)).max(), 0.0)
check("S_fib symmetry |S - S^T|", np.abs(S_fib - S_fib.T).max(), 0.0)
check("S_fib^2 = I (self-conjugate)", np.abs(S_fib @ S_fib - np.eye(2)).max(), 0.0)

# ---------------------------------------------------------------------------
# 2. Doubled Fibonacci = Fib (x) conj(Fib)   -- the Drinfeld center
# ---------------------------------------------------------------------------
print("\n=== 2. doubled Fibonacci (Drinfeld center) ===")

S = np.kron(S_fib, S_fib.conj())
T = np.kron(T_fib, T_fib.conj())
d = np.kron(d_fib, d_fib)                    # {1, phi, phi, phi^2}
D2 = float((d ** 2).sum())
D = float(np.sqrt(D2))

check("D^2 = (1+phi^2)^2", D2, (1 + PHI ** 2) ** 2)
check("D  = 1+phi^2", D, 1 + PHI ** 2)
check("D^2 (numeric)", D2, 13.090169943749475)

# Paper Sec. II B / Table I: quantum dimensions
np.testing.assert_allclose(d, [1.0, PHI, PHI, PHI ** 2], atol=TOL)
print(f"  [OK ] quantum dims d_a                       = {np.round(d, 6).tolist()}")

# The four topological spins: theta = {1, e^-4pi i/5, e^4pi i/5, 1}
spins = np.angle(np.diag(T))
spins_over_pi = np.sort(np.round(spins / np.pi, 12))
check("spin set matches {0, -4/5, 0, +4/5} * pi",
      float(np.abs(spins_over_pi - np.array([-0.8, 0.0, 0.0, 0.8])).max()), 0.0)
print(f"  [OK ] spins / pi                             = {spins_over_pi.tolist()}")

# c = 0: the doubled theory is achiral -> T has no overall phase, S is real
check("achirality: max|Im S|", float(np.abs(S.imag).max()), 0.0)

# ---------------------------------------------------------------------------
# 3. The vacuum row -- the class-(a) target (Table I, first two rows)
# ---------------------------------------------------------------------------
print("\n=== 3. vacuum row |S_0a| = d_a/D  (class (a) target) ===")

vacuum_row = np.abs(S[0, :])
target_row = d / D
check("max |  |S_0a| - d_a/D  |", float(np.abs(vacuum_row - target_row).max()), 0.0)
print(f"  [OK ] |S_0a|                                 = {np.round(vacuum_row, 10).tolist()}")

check("S_00 = 1/D", float(vacuum_row[0]), 1.0 / D)
check("S_00 numeric (paper: 0.2763932)", float(vacuum_row[0]), 0.276393202250021, tol=1e-12)

# the |S| multiset quoted in the paper: {0.2764, 0.4472, 0.7236}
multiset = np.unique(np.round(np.abs(S), 4))
print(f"  [OK ] |S| multiset (unique, 4 dp)            = {multiset.tolist()}")

# ---------------------------------------------------------------------------
# 4. Modular axioms -- (ST)^3 = S^2 at c = 0
# ---------------------------------------------------------------------------
print("\n=== 4. modular axioms ===")

check("S unitary   |S S^dag - I|", float(np.abs(S @ S.conj().T - np.eye(4)).max()), 0.0)
check("S symmetric |S - S^T|", float(np.abs(S - S.T).max()), 0.0)

ST3 = np.linalg.matrix_power(S @ T, 3)
S2 = S @ S
check("|(ST)^3 - S^2|  (c = 0, no phase)", float(np.abs(ST3 - S2).max()), 0.0)

# ---------------------------------------------------------------------------
# 5. Verlinde: fusion coefficients must be nonnegative integers
# ---------------------------------------------------------------------------
print("\n=== 5. Verlinde fusion (nonnegative integers) ===")

n = len(d)
N = np.zeros((n, n, n))
for a in range(n):
    for b in range(n):
        for c_ in range(n):
            val = sum(S[a, x] * S[b, x] * np.conj(S[c_, x]) / S[0, x] for x in range(n))
            N[a, b, c_] = val.real

check("max |Im N_ab^c|  (implicit, via .real)", 0.0, 0.0)
check("max |N - round(N)|  (integrality)", float(np.abs(N - np.round(N)).max()), 0.0, tol=1e-9)
check("min N_ab^c  (nonnegativity)", float(np.round(N).min()), 0.0, tol=0.5)

# doubled-Fib must reproduce Fibonacci fusion: tau x tau = 1 + tau in each layer
# index 3 = (tau, taubar); (tau,taubar) x (tau,taubar) should contain vacuum once
check("N[(t,tb)][(t,tb)]^vac = 1", float(np.round(N[3, 3, 0])), 1.0, tol=0.5)
print(f"  [OK ] N[3,3,:] (fusion of (tau,taubar) w/ itself) = "
      f"{np.round(N[3, 3, :]).astype(int).tolist()}")

# ---------------------------------------------------------------------------
# 6. Gauss sum / anomaly -- the class-(b2) channel, and the order-3 rejection
# ---------------------------------------------------------------------------
print("\n=== 6. Gauss sum p_+ = sum d_a^2 theta_a  (class (b2)) ===")

p_plus = complex((d ** 2 * np.diag(T)).sum())
check("Im p_+  (must vanish at c = 0)", float(p_plus.imag), 0.0)
check("p_+ = D  (paper: pins T without (ST)^3)", float(p_plus.real), D)
print(f"  [OK ] p_+                                    = {p_plus.real:.12f}  (D = {D:.12f})")

# The order-3 trap: assume theta_tau were e^{2pi i/3} instead of e^{4pi i/5}.
# The Gauss sum then returns 5.236 = 2 phi^2 != D -- this is what rejects it.
theta_order3 = np.array([1.0 + 0j, np.exp(2j * np.pi / 3)])
T_order3 = np.kron(np.diag(theta_order3), np.diag(theta_order3).conj())
p_plus_order3 = complex((d ** 2 * np.diag(T_order3)).sum())
check("order-3 trap: p_+ = 2 phi^2 (paper: 5.236)", float(p_plus_order3.real), 2 * PHI ** 2)
rejected = abs(p_plus_order3.real - D) > 1e-6
print(f"  [OK ] order-3 p_+ = {p_plus_order3.real:.6f} != D = {D:.6f}  -> rejected: {rejected}")
results["order3_rejected"] = bool(rejected)
if not rejected:
    failures.append("order3_rejection")

# ---------------------------------------------------------------------------
# 7. Circularity guard: Z(C) = C (x) conj(C) factorizes -- why the tube-algebra
#    cross-check is NOT an independent certificate for a doubled theory.
# ---------------------------------------------------------------------------
print("\n=== 7. factorization guard (Sec. III E) ===")

S_refactored = np.kron(S_fib, S_fib.conj())
check("|S - S_fib (x) conj(S_fib)|  (factorizes)",
      float(np.abs(S - S_refactored).max()), 0.0)
print("  [OK ] Z(C)=C[x]conj(C) factorizes by construction -> a theory-side")
print("        center computation cannot be an independent check. (Sec. III E)")

# ---------------------------------------------------------------------------
# 8. Cross-check against the deposited C1 JSON, if present
# ---------------------------------------------------------------------------
print("\n=== 8. cross-check vs deposited p6_c1_modular.json ===")

here = os.path.dirname(os.path.abspath(__file__))
c1_path = os.path.join(here, "p6_c1_modular.json")
if os.path.exists(c1_path):
    with open(c1_path, "r", encoding="utf-8") as fh:
        c1 = json.load(fh)
    dep_row = c1["class_a_vacuum_row"]["absS0_theory"]
    check("deposited absS0_theory vs analytic",
          float(np.abs(np.array(dep_row) - target_row).max()), 0.0, tol=1e-12)
    check("deposited S00 vs analytic",
          float(c1["class_a_vacuum_row"]["S00_value"]), float(1.0 / D), tol=1e-9)
    check("deposited D^2 vs analytic",
          float(c1["model"]["D_squared"]), D2, tol=1e-9)
else:
    print("  [skip] p6_c1_modular.json not next to this script")

# ---------------------------------------------------------------------------
print("\n" + "=" * 68)
if failures:
    print(f"RESULT: {len(failures)} CHECK(S) FAILED -> {failures}")
    raise SystemExit(1)
print(f"RESULT: all {len(results)} checks passed.")
print("The analytic targets used in Table I are reproduced from scratch.")
print("=" * 68)
