"""
p6_verify_patch.py — REFEREE HARNESS for any extractor patch.

RULE (16 Jul): a patch counts as UNVERIFIED until it passes this
harness — no matter how convincing the diff looks.

TRUST BOUNDARY: this harness imports from the patch module ONLY
extract(M_wilson, M_rot, ...). Theory matrices, synthetic-data generation,
comparison, canonicalization, and ALL thresholds live HERE — a patch cannot
soften them. (synth()/compare_to_theory() from the module are deliberately NOT
used either: otherwise a patch could supply its own test data or its own yardstick.)

GATES (pre-registered acceptance gates; P1' re-registered 16 Jul):
  P1' toric-code synthetic  -> (a) toric among the solutions (<1e-6)
                               (b) EVERY solution in the two SOLUTION CLASSES {toric, double semion}
                               (c) every solution Verlinde-/(ST)^3-valid
                               (d) both classes populated (toric AND double semion
                                   must appear — a toric-only extractor would carry a
                                   silent selection bias).
                                   NOTE — CORRECTION 16 Jul (internal-review finding, confirmed):
                                   This is NOT the fusion-phase orbit of toric. The true
                                   IV.3 orbit is {toric, 3-FERMION}; double semion does NOT
                                   lie in it (separable by loops, theta_s^2 = -1 vs +1). {toric, DS}
                                   are the two SOLUTION CLASSES of OUR inverse problem
                                   (one loop spectrum + M_R): DS' (M_w, S T*) is unitarily
                                   equivalent to toric's, while the true orbit partner 3-fermion is
                                   excluded by the rotation data (tr(M_w R^2)=-1 vs +1).
                                   "DS self-conjugate => 2 classes" was the WRONG justification
                                   (two distinct ambiguities happen to both give 2 classes).
  P2  dFib synthetic        -> UNCHANGED strict: match < 1e-6, ALL solutions equal the true one
                               (fusion-phase freedom trivial, Abelian subcategory A = {1})
  P3  random orthogonal     -> must FAIL (0 solutions)
  P3b random unitary        -> must FAIL (0 solutions)
  P4  circularity canary: module theory functions patched to junk;
      extract() must return the IDENTICAL result (no theory inside the extractor)
  P5  scaling gate (pre-registered): input noise 1e-4/1e-3/1e-2 -> the residual moves along
      (run only if P1+P2 pass; mandatory BEFORE any lattice run)

API CONTRACT with the patch: extract(M_wilson, M_rot, n_starts=..., seed=...) ->
  (sols, blocks, R0) with sols = list of dicts with keys 'S' (4x4 complex)
  and 'T' (4x4 diagonal complex). Everything else may change.

Usage:  python p6_verify_patch.py [modulename]     (default: p6_extract_c3)
"""
import sys, os, json, itertools, importlib
import numpy as np

MOD = sys.argv[1] if len(sys.argv) > 1 else "p6_extract_c3"
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
patch = importlib.import_module(MOD)

rng_global = np.random.default_rng(20260716)

# ---------------- OWN references (not from the patch module!) ----------------
def th_toric():
    S = 0.5 * np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]], dtype=complex)
    return S, np.diag([1, 1, 1, -1]).astype(complex)

def th_dfib():
    phi = (1 + np.sqrt(5)) / 2; D = np.sqrt(1 + phi ** 2)
    Sf = np.array([[1, phi], [phi, -1]]) / D; Tf = np.diag([1, np.exp(4j * np.pi / 5)])
    return np.kron(Sf, np.conj(Sf)), np.kron(Tf, np.conj(Tf))

def my_synth(S, T, x, seed_w=7, noise=0.0):
    """Own synthetic data: Wilson loop = fusion character, rotation = S T* (Li-Mong Eq.63),
    both in a randomly twisted basis. Optional relative noise on BOTH."""
    lam = S[x, :] / S[0, :]
    Mw = np.diag(lam); R = S @ np.conj(T)
    rw = np.random.default_rng(seed_w)
    W = np.linalg.qr(rw.normal(size=(4, 4)) + 1j * rw.normal(size=(4, 4)))[0]
    Mw, R = W @ Mw @ W.conj().T, W @ R @ W.conj().T
    if noise > 0:
        rn = np.random.default_rng(seed_w + 1000)
        Nw = rn.normal(size=(4, 4)) + 1j * rn.normal(size=(4, 4))
        Nr = rn.normal(size=(4, 4)) + 1j * rn.normal(size=(4, 4))
        Mw = Mw + noise * np.linalg.norm(Mw) / np.linalg.norm(Nw) * (Nw + Nw.conj().T) / 2
        R = R + noise * np.linalg.norm(R) / np.linalg.norm(Nr) * Nr
    return Mw, R

def my_match(S_ex, T_ex, S_th, T_th):
    best = 1e9
    for perm in itertools.permutations(range(4)):
        P = np.eye(4)[list(perm)]
        Sp, Tp = P @ S_ex @ P.T, P @ T_ex @ P.T
        for conj in (False, True):
            Sc, Tc = (np.conj(Sp), np.conj(Tp)) if conj else (Sp, Tp)
            d = np.linalg.norm(Sc - S_th) + np.linalg.norm(Tc - T_th)
            best = min(best, float(d))
    return best

results = {}
def gate(name, ok, detail):
    results[name] = {"pass": bool(ok), "detail": detail}
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

print("=" * 76)
print(f"REFEREE HARNESS against module '{MOD}'  (thresholds live in the harness, not in the patch)")
print("=" * 76)

N_STARTS = 150

# ---------------- references for P1' (each self-validated BEFORE use — review lesson:
# ---------------- a "wrongly guessed S convention" is caught by self-validation)
def th_dsemion():
    Ss = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    Ts = np.diag([1, 1j])
    return np.kron(Ss, np.conj(Ss)), np.kron(Ts, np.conj(Ts))

def mtc_valid(S, T, tol=1e-9):
    """Generic self-validation: unitary, symmetric, S^2=perm, (ST)^3=phi S^2, Verlinde-integer."""
    if np.linalg.norm(S.conj().T @ S - np.eye(4)) > tol: return False, "unitary"
    if np.linalg.norm(S - S.T) > tol: return False, "symmetric"
    P2 = np.abs(S @ S)
    if np.linalg.norm(P2 - (P2 > 0.5)) > tol: return False, "S2perm"
    ST = S @ T; ST3 = ST @ ST @ ST; S2m = S @ S
    idx = np.unravel_index(np.argmax(np.abs(S2m)), S2m.shape)
    ph = ST3[idx] / S2m[idx]
    if np.linalg.norm(ST3 - ph * S2m) > tol or abs(abs(ph) - 1) > tol: return False, "ST3"
    v = None
    for a in range(4):
        if np.all(S[a, :].real > 1e-9) and np.max(np.abs(S[a, :].imag)) < tol: v = a; break
    if v is None: return False, "no vacuum"
    N = np.einsum('ax,bx,cx,x->abc', S, S, np.conj(S), 1.0 / S[v, :])
    if np.abs(N - np.round(N.real)).max() > 1e-6 or np.round(N.real).min() < -0.5: return False, "Verlinde"
    return True, "ok"

S_ds, T_ds = th_dsemion()
ok_ds, why_ds = mtc_valid(S_ds, T_ds)
S_tcv, T_tcv = th_toric()
ok_tc, why_tc = mtc_valid(S_tcv, T_tcv)
print(f"  [ref] toric self-validated: {ok_tc} ({why_tc}) | double semion: {ok_ds} ({why_ds})")
ds_selfconj = my_match(np.conj(S_ds), np.conj(T_ds), S_ds, T_ds) < 1e-9
print(f"  [ref] conj(DS) ~ DS per Permutation: {ds_selfconj} "
      f"-> naive class count would be {'2' if ds_selfconj else '3'} classes (diagnostic; see correction above, NOT the reason for 2 classes)")
assert ok_ds and ok_tc, "reference invalid — harness abort instead of a wrong gate"

# ---------------- P1' (re-registered 16.7.) ----------------
S_th, T_th = th_toric()
Mw, Mr = my_synth(S_th, T_th, x=1)
sols, _, _ = patch.extract(Mw, Mr, n_starts=N_STARTS)
if not sols:
    gate("P1_toric_classes", False, "0 solutions")
else:
    d_tc = [my_match(s["S"], s["T"], S_th, T_th) for s in sols]
    d_ds = [my_match(s["S"], s["T"], S_ds, T_ds) for s in sols]
    valid = [mtc_valid(s["S"], s["T"], tol=1e-6)[0] for s in sols]
    in_classes = [min(a, b) < 1e-6 for a, b in zip(d_tc, d_ds)]
    a_ok = min(d_tc) < 1e-6                       # toric found
    b_ok = all(in_classes)                          # every solution in one of the two classes
    c_ok = all(valid)                             # every solution MTC-valid
    d_ok = a_ok and (min(d_ds) < 1e-6)            # both classes populated: toric AND DS
    gate("P1_toric_classes", a_ok and b_ok and c_ok and d_ok,
         f"{len(sols)} solutions | (a) toric best {min(d_tc):.2e} | (b) all in the 2 classes: {b_ok} "
         f"| (c) all MTC-valid: {c_ok} | (d) DS best {min(d_ds):.2e} -> both classes: {d_ok}")

# ---------------- P2 (unchanged, strict) ----------------
S_th, T_th = th_dfib()
Mw, Mr = my_synth(S_th, T_th, x=3)
sols, _, _ = patch.extract(Mw, Mr, n_starts=N_STARTS)
if not sols:
    gate("P2_dfib", False, "0 solutions")
else:
    ds = [my_match(s["S"], s["T"], S_th, T_th) for s in sols]
    best, worst = min(ds), max(ds)
    gate("P2_dfib", best < 1e-6 and worst < 1e-6,
         f"{len(sols)} solutions, match best {best:.2e} / worst {worst:.2e} (both < 1e-6 required; "
         f"worst probes uniqueness: EVERY solution must be the true one)")

# ---------------- P3 / P3b ----------------
S_th, T_th = th_dfib(); Mw, _ = my_synth(S_th, T_th, x=3)
Q = np.linalg.qr(rng_global.normal(size=(4, 4)))[0]
sb, _, _ = patch.extract(Mw, Q, n_starts=N_STARTS)
gate("P3_antifudge_orthogonal", len(sb) == 0, f"{len(sb)} solutions from junk (0 required)")
Qu = np.linalg.qr(rng_global.normal(size=(4, 4)) + 1j * rng_global.normal(size=(4, 4)))[0]
sbu, _, _ = patch.extract(Mw, Qu, n_starts=N_STARTS)
gate("P3b_antifudge_unitary", len(sbu) == 0, f"{len(sbu)} solutions from junk (0 required)")

# ---------------- P4: circularity canary ----------------
import inspect
sig = list(inspect.signature(patch.extract).parameters)
sig_ok = not any(("theo" in p.lower()) or p.startswith("S_") or p.startswith("T_") for p in sig)
Mw2, Mr2 = my_synth(*th_dfib(), x=3)
ref_sols, _, _ = patch.extract(Mw2, Mr2, n_starts=40)
saved = {}
for fn in ("theory_doubled_fib", "theory_toric"):
    if hasattr(patch, fn):
        saved[fn] = getattr(patch, fn)
        setattr(patch, fn, lambda *a, **k: (_ for _ in ()).throw(RuntimeError("theory poisoned")))
try:
    poisoned_sols, _, _ = patch.extract(Mw2, Mr2, n_starts=40)
    same = len(ref_sols) == len(poisoned_sols) and all(
        np.allclose(a["S"], b["S"]) and np.allclose(a["T"], b["T"])
        for a, b in zip(ref_sols, poisoned_sols))
    gate("P4_circularity_canary", sig_ok and same,
         f"signature free of theory args: {sig_ok}; result identical under poisoned theory: {same}")
except RuntimeError as e:
    gate("P4_circularity_canary", False, f"extract() CALLED the theory function: {e}")
finally:
    for fn, f in saved.items():
        setattr(patch, fn, f)

# ---------------- P5: scaling gate (pre-registered) — only if P1+P2 pass ----------------
if results.get("P1_toric_classes", {}).get("pass") and results.get("P2_dfib", {}).get("pass"):
    print("\n  [P5] scaling check: the residual must move with the input noise")
    S_th, T_th = th_dfib()
    scal = {}
    for noise in (1e-4, 1e-3, 1e-2):
        Mw, Mr = my_synth(S_th, T_th, x=3, noise=noise)
        sols, _, _ = patch.extract(Mw, Mr, n_starts=N_STARTS)
        d = min((my_match(s["S"], s["T"], S_th, T_th) for s in sols), default=np.nan)
        scal[noise] = float(d)
        print(f"        noise {noise:.0e} -> residual {d:.3e}  ({len(sols)} solutions)")
    vals = [scal[n] for n in (1e-4, 1e-3, 1e-2)]
    monotone = (not any(np.isnan(v) for v in vals)) and vals[0] < vals[1] < vals[2]
    within = all(not np.isnan(v) and v < 30 * n for n, v in scal.items())
    gate("P5_scaling", monotone and within,
         f"monotone {monotone}, residual ~ O(noise) {within}: {scal}")
else:
    gate("P5_scaling", False, "skipped — P1/P2 not green")

# ---------------- verdict ----------------
required = ["P1_toric_classes", "P2_dfib", "P3_antifudge_orthogonal", "P3b_antifudge_unitary", "P4_circularity_canary"]
ok = all(results[g]["pass"] for g in required)
lattice_ready = ok and results["P5_scaling"]["pass"]
print("\n" + "=" * 76)
print(f"VERDICT: patch {'ACCEPTED' if ok else 'REJECTED'} "
      f"({sum(results[g]['pass'] for g in required)}/{len(required)} mandatory gates)")
print(f"LATTICE CLEARANCE (incl. P5): {'YES' if lattice_ready else 'NO'}")
with open(os.path.join(H, f"verify_{MOD}.json"), "w") as f:
    json.dump({"module": MOD, "gates": results, "accepted": ok,
               "lattice_ready": lattice_ready}, f, indent=2)
print(f"EXPORT OK -> verify_{MOD}.json")
