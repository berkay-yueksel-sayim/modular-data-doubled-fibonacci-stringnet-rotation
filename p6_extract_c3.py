"""
p6_extract_c3.py — rotation-channel extractor (single-model test, 2026-07-16).
UNVERIFIED until p6_verify_patch.py reports 5/5 + P5.

WHAT THIS FIX CHANGES OVER THE PREDECESSOR:
  The old constraint set ("S symmetric") was underdetermined — many symmetric
  unitary S solve it. Two GENERIC MTC bindings were missing:

  (a) FUSION-CHARACTER BINDING (theta-free!):
      The measured Wilson loop has eigenvalues lambda_a = c * S_{x,a}/S_{v,a}
      (v = vacuum position, x = anyon of the loop, c = normalization).
      Since S_{x,a}/S_{v,a} = A_{x,a}/A_{v,a} (the spins theta cancel), this becomes
          lambda_a * A_{v,a} = c * A_{x,a}          (A = rotation matrix in the MES basis)
      — a condition DIRECTLY on A, tying the measured lambda to the rows of S.

  (b) VACUUM-ROW POSITIVITY: S_{v,a} = d_a * S_{v,v} > 0 real (quantum dimensions).
      Kills the row-phase freedom that inflated the solution set.

  Both are axioms of EVERY unitary MTC — no doubled-Fib input, no theory S/T.

CIRCULARITY ARCHITECTURE (unchanged): extract(M_wilson, M_rot, ...) receives
ONLY measured matrices. Theory lives solely in the separate theory_*() helpers,
which extract() never calls (harness gate P4 verifies this by poisoning them).

Algorithm (Li-Mong 2203.04329 V.1, Eqs. 60/62/63):
  A := B^dag (E^dag M_rot E) B  should equal  S*conj(T)  (Eq. 63) =>  S = A*T.
  S symmetric  <=>  A_ab t_b = A_ba t_a.
  Parameters: B = block unitary on the eigenspaces of the Wilson loop (Li-Mong's
  "phase degree of freedom") + 3 spins (t_v = 1 fixed, vacuum-spin axiom).
  Outer loop over the vacuum position v (round-robin over the starts) — v is
  NOT guessed from theory; every hypothesis must solve the constraints.
"""
import itertools
import numpy as np
from scipy.optimize import minimize, least_squares
from scipy.linalg import expm

# ---------------------------------------------------------------- building blocks
def _blocks(evals, tol=1e-6):
    b, i = [], 0
    while i < len(evals):
        j = i
        while j + 1 < len(evals) and abs(evals[j + 1] - evals[i]) < tol: j += 1
        b.append(list(range(i, j + 1))); i = j + 1
    return b

def _B_of(p, blocks, n=4):
    B = np.zeros((n, n), dtype=complex); k = 0
    for blk in blocks:
        d = len(blk)
        h = np.zeros((d, d), dtype=complex); idx = 0
        for a in range(d):
            h[a, a] = p[k + idx]; idx += 1
        for a in range(d):
            for b in range(a + 1, d):
                h[a, b] = p[k + idx] + 1j * p[k + idx + 1]
                h[b, a] = np.conj(h[a, b]); idx += 2
        k += d * d
        Bl = expm(1j * h)
        for x, a in enumerate(blk):
            for y, b in enumerate(blk): B[a, b] = Bl[x, y]
    return B

def _char_cost(A, lam, v):
    """(a) Fusion-character binding, theta-free: lambda_a A_{v,a} = c A_{x,a}.
    c analytic via least squares, min over the candidate row x."""
    w = lam * A[v, :]
    best = np.inf
    for x in range(A.shape[0]):
        u = A[x, :]
        den = np.vdot(u, u).real
        if den < 1e-300: continue
        c = np.vdot(u, w) / den
        r = np.linalg.norm(w - c * u) ** 2 / max(np.linalg.norm(w) ** 2, 1e-300)
        best = min(best, float(r))
    return best

def _assemble(p, R0, blocks, lam, v, nB):
    B = _B_of(p[:nB], blocks)
    A = B.conj().T @ R0 @ B
    t = np.empty(4, dtype=complex)
    phases = iter(p[nB:])
    for a in range(4):
        t[a] = 1.0 if a == v else np.exp(1j * next(phases))   # t_v = 1 (vacuum-spin axiom)
    return A, t, A * t[None, :]

def _resid(p, R0, blocks, lam, v, nB):
    """Residual VECTOR (for Levenberg-Marquardt): sym + vac + char, all stacked real."""
    A, t, S = _assemble(p, R0, blocks, lam, v, nB)
    M = S - S.T
    row = S[v, :]
    # character binding: lambda_a A_{v,a} = c A_{x,a}; x + c analytic (best pair)
    w = lam * A[v, :]
    best, bx, bc = np.inf, 0, 0.0
    for x in range(4):
        u = A[x, :]
        den = np.vdot(u, u).real
        if den < 1e-300: continue
        c = np.vdot(u, w) / den
        rr = float(np.linalg.norm(w - c * u) ** 2) / max(float(np.linalg.norm(w) ** 2), 1e-300)
        if rr < best: best, bx, bc = rr, x, c
    ch = (w - bc * A[bx, :]) / max(np.linalg.norm(w), 1e-150)
    # S^2 = C (charge conjugation = real permutation) — encoded SMOOTHLY (16 Jul, fix 3):
    # without these terms, {sym+unitary+vac+char} has an entire branch of EXACT non-MTC
    # solutions (an LM candidate with fun 4e-30 but S2-perm 0.59, (ST)^3 1.58, Verlinde 2.6)
    # — the optimizer lands there instead of on the modular branch. Generic: this is the
    # DEFINITION of the object sought, applied identically to every input (gate P3!).
    U2 = S @ S
    m2 = np.abs(U2) ** 2
    perm_mag = (m2 * (1.0 - m2)).ravel()                  # |S^2| entries -> {0,1}
    perm_im  = U2.imag.ravel()                            # S^2 real
    # (ST)^3 = phi * S^2 with |phi| = 1 (modularity; phi = central charge, free)
    ST = S * t[None, :]                                   # == S @ diag(t)
    ST3 = ST @ ST @ ST
    den = np.vdot(U2, U2).real
    phi = np.vdot(U2, ST3) / den if den > 1e-300 else 1.0
    st3 = (ST3 - phi * U2)
    return np.concatenate([M.real.ravel(), M.imag.ravel(),
                           row.imag, np.minimum(row.real, 0.0),
                           ch.real, ch.imag,
                           perm_mag, perm_im,
                           st3.real.ravel(), st3.imag.ravel(),
                           [abs(phi) - 1.0]])

def _cost_weak(p, R0, blocks, lam, v, nB):
    """Stage A (global search): ONLY sym + vac + char — smooth landscape, ~27% basin
    hit rate. Its solution branch CONTAINS the true point; the modular terms
    (perm/st3) are too rough for BFGS (1/25 instead of 8/30 hits, measured 16 Jul)."""
    A, t, S = _assemble(p, R0, blocks, lam, v, nB)
    M = S - S.T
    row = S[v, :]
    return float(np.sum(np.abs(M) ** 2) + np.sum(row.imag ** 2)
                 + np.sum(np.minimum(row.real, 0.0) ** 2) + _char_cost(A, lam, v))

def _cost(p, R0, blocks, lam, v, nB):
    r = _resid(p, R0, blocks, lam, v, nB)
    return float(np.dot(r, r))

# ---------------------------------------------------------------- EXTRACTOR
def extract(M_wilson, M_rot, n_starts=400, seed=2026, input_scale=None):
    """ONLY measured matrices in. No theory argument — by design.

    DATA-ADAPTIVE tolerances: the noise level eps is estimated from the
    unitarity deviation of the MEASURED rotation (a property of the data, not
    of the desired result). Exact input => eps~1e-15 => strict thresholds;
    noisy input => thresholds scale along (required for the pre-registered
    scaling gate P5: the residual must move with the input noise)."""
    wv, E = np.linalg.eigh(M_wilson)
    R0 = E.conj().T @ M_rot @ E
    lam = wv.astype(complex)
    eps = max(1e-13, float(np.linalg.norm(R0.conj().T @ R0 - np.eye(4))))
    # input_scale: optional DEPOSITED, data-derived input error (e.g. the loop
    # ideality 3.2668e-3 from p6_loop_ideality.json for the lattice run, where the
    # Wilson side is invisible in the R0 unitarity). NEVER a desired value.
    if input_scale is not None:
        eps = max(eps, float(input_scale))
    # EPS-AWARE block detection (16 Jul, fix 5): noise splits the degenerate pair
    # by ~noise; a fixed 1e-6 tolerance broke up the block -> the U(2) freedom of the
    # MES basis was lost -> 0 solutions under ANY noise (P5 finding). 3*eps stays
    # well below the true gaps (dFib: 1.38/2.24).
    blocks = _blocks(wv, tol=max(1e-6, 3 * eps))
    tol_accept = max(1e-14, (10 * eps) ** 2)   # headroom: residual vector has ~100 components of O(eps) each
    tol_polish = max(1e-6, (30 * eps) ** 2)   # BFGS-endpoint trigger for LM (junk ~0.4 stays out)
    tol_S      = max(1e-7, 5 * eps)
    tol_perm   = max(1e-6, 10 * eps)
    tol_st3    = max(1e-6, 10 * eps)
    tol_verl   = max(1e-4, 50 * eps)
    dec        = int(max(2, min(7, -np.log10(max(eps, 1e-8)))))
    nB = sum(len(b) ** 2 for b in blocks)
    npar = nB + 3
    rng = np.random.default_rng(seed)
    sols = []
    for s in range(n_starts):
        v = s % 4                                              # round-robin vacuum hypothesis
        p0 = rng.uniform(-np.pi, np.pi, npar)
        # STAGE A (global): BFGS on the WEAK cost function (smooth, high hit rate).
        r = minimize(_cost_weak, p0, args=(R0, blocks, lam, v, nB), method="BFGS",
                     options={"maxiter": 2000, "gtol": 1e-12})
        if r.fun > tol_polish: continue
        # STAGE B (local): LEVENBERG-MARQUARDT on the AUGMENTED residual vector
        # (sym+vac+char+S2perm+(ST)3). LM slides along the weak solution branch
        # to the modular point; BFGS strands there via precision loss (16 Jul, fixes 2+3).
        ls = least_squares(_resid, r.x, args=(R0, blocks, lam, v, nB), method="lm",
                           xtol=2.3e-16, ftol=2.3e-16, gtol=2.3e-16, max_nfev=20000)
        fun = float(np.dot(ls.fun, ls.fun))
        # LM RESTART for the last decade: LM occasionally strands just above the
        # acceptance (observed: true dFib solution at 2.8e-14 vs. acceptance 1e-14).
        # Restart with a fresh trust-region radius, NO threshold change.
        for _ in range(2):
            if fun < tol_accept or fun > 1e-8: break
            ls = least_squares(_resid, ls.x, args=(R0, blocks, lam, v, nB), method="lm",
                               xtol=2.3e-16, ftol=2.3e-16, gtol=2.3e-16, max_nfev=20000)
            fun = float(np.dot(ls.fun, ls.fun))
        if fun >= tol_accept: continue
        A, t, S = _assemble(ls.x, R0, blocks, lam, v, nB)
        # ---- post-filter: ONLY generic MTC axioms (no theory S/T, no F/R) ----
        if np.linalg.norm(S - S.T) > tol_S: continue                       # S symmetric
        if np.linalg.norm(S.conj().T @ S - np.eye(4)) > tol_S: continue    # S unitary
        if np.min(S[v, :].real) < 1e-9 or np.max(np.abs(S[v, :].imag)) > tol_S: continue  # vacuum row > 0
        S2 = np.abs(S @ S)                                                 # S^2 = charge conjugation
        if np.linalg.norm(S2 - (S2 > 0.5).astype(float)) > tol_perm: continue
        if not np.allclose((S2 > 0.5).sum(axis=0), 1) or not np.allclose((S2 > 0.5).sum(axis=1), 1):
            continue
        # (ST)^3 = e^{i phi} S^2 — modularity relation (generic; phase = central
        # charge, left free). NEW 16 Jul: kills the P1 junk class (violation there 0.6-1.4).
        ST = S @ np.diag(t)
        ST3 = ST @ ST @ ST; S2m = S @ S
        idx = np.unravel_index(np.argmax(np.abs(S2m)), S2m.shape)
        phz = ST3[idx] / S2m[idx] if abs(S2m[idx]) > 1e-9 else 1.0
        if abs(abs(phz) - 1) > tol_S or np.linalg.norm(ST3 - phz * S2m) > tol_st3: continue
        # Verlinde integrality: N_ab^c = sum_x S_ax S_bx S*_cx / S_vx must be non-negative
        # integer (generic). NEW 16 Jul: kills the P1 junk class (there 0.48..2e8).
        # Verlinde: N_ab^c = sum_x S_ax S_bx S*_cx / S_{v,x} — division INSIDE the
        # x sum! (16 Jul, fix 4: the old version divided by S_{v,c} AFTER the sum;
        # identical for a constant vacuum row (toric, all 1/2) -> undetected, but for a
        # non-constant one (dFib) it filtered out the TRUE solution with spurious dev 0.48.)
        N = np.einsum('ax,bx,cx,x->abc', S, S, np.conj(S), 1.0 / S[v, :])
        Nr = np.round(N.real)
        if np.abs(N - Nr).max() > tol_verl or Nr.min() < -0.5: continue
        # dedup (tolerance-aware)
        key = (np.round(S, dec).tobytes(), np.round(t, dec).tobytes())
        if any(k == key for k, *_ in sols): continue
        sols.append((key, {"S": S, "T": np.diag(t), "vac": v, "cost": fun,
                           "eps_input": eps}))
    return [d for _, d in sols], blocks, R0

def canonical(S, T, tol=1e-6):
    th = np.diag(T)
    return (tuple(np.round(np.sort_complex(th), 6)),
            tuple(np.round(np.sort(np.abs(S).ravel()), 6)),
            tuple(np.round(np.sort_complex(np.linalg.eigvals(S)), 6)))

# ------------------------------------------------- theory/comparison (SEPARATE)
def theory_doubled_fib():
    phi = (1 + np.sqrt(5)) / 2; D = np.sqrt(1 + phi ** 2)
    Sf = np.array([[1, phi], [phi, -1]]) / D; Tf = np.diag([1, np.exp(4j * np.pi / 5)])
    return np.kron(Sf, np.conj(Sf)), np.kron(Tf, np.conj(Tf))

def theory_toric():
    S = 0.5 * np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]], dtype=complex)
    return S, np.diag([1, 1, 1, -1]).astype(complex)

def compare_to_theory(S_ex, T_ex, S_th, T_th):
    best = (1e9, None, None)
    for perm in itertools.permutations(range(4)):
        P = np.eye(4)[list(perm)]
        Sp, Tp = P @ S_ex @ P.T, P @ T_ex @ P.T
        for conj in (False, True):
            Sc, Tc = (np.conj(Sp), np.conj(Tp)) if conj else (Sp, Tp)
            d = np.linalg.norm(Sc - S_th) + np.linalg.norm(Tc - T_th)
            if d < best[0]: best = (float(d), perm, conj)
    return best

def synth(S, T, x=1):
    lam = S[x, :] / S[0, :]
    M_w = np.diag(lam); R = S @ np.conj(T)
    W = np.linalg.qr(np.random.default_rng(7).normal(size=(4, 4)) +
                     1j * np.random.default_rng(8).normal(size=(4, 4)))[0]
    return W @ M_w @ W.conj().T, W @ R @ W.conj().T
