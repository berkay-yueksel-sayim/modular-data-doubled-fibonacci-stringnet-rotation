"""
p6_reconciliation_test.py -- reconciliation, second attempt (16 Jul). Replaces a
predecessor diagnostic (not part of this record) whose construction was
physically inconsistent: there lam was SORTED while the rotation was left in the
anyon ordering -> a basis vector carried the Wilson eigenvalue of anyon a but the
rotation matrix elements of anyon b. No MTC has that pairing; the 0 solutions
there were the extractor CORRECTLY refusing impossible data (not a physics
finding). Here: everything consistent in the ANYON FRAME.

PRE-REGISTERED PREDICTION (fixed before the run, mechanism argument):
  lam enters ONLY through the character-binding term. At the true point its
  residual is irreducibly ~strength*delta and NORM-INVARIANT under the whole
  B freedom, BECAUSE the real delta form is equal on the degenerate pair
  (measured: 5.165338968e-3 vs 5.165339017e-3, equal to 5e-11). No gradient
  -> no trade-off -> S,T come from the (exact) rotation alone.
  EXPECTATION part 1 (eigenvalue channel, real form): residual FLAT ~LM floor,
    independent of the strength.
  EXPECTATION part 2 (eigenvector channel):           residual SCALES ~ strength.
    That is the control that CAN FAIL -- without it part 1 would be decoration.
  EXPECTATION part 3 (invariance demo):               character residual constant
    under random B, sym residual varies O(1). Mechanism shown directly.

Deterministic. Read-only on the deposits.
"""
import os, json, itertools
import numpy as np
from scipy.linalg import expm
from p6_extract_c3 import extract, _blocks, _B_of

H = os.path.dirname(os.path.abspath(__file__))
dep = json.load(open(os.path.join(H, "p6_loop_ideality.json")))
rA = dep["load_bearing"]["route_A_eigh_Mj"]
INPUT_DEV = rA["rel_residual"]                                  # 3.2668e-3

phi = (1 + np.sqrt(5)) / 2; D = np.sqrt(1 + phi ** 2)
Sf = np.array([[1, phi], [phi, -1]]) / D; Tf = np.diag([1, np.exp(4j * np.pi / 5)])
S_th, T_th = np.kron(Sf, np.conj(Sf)), np.kron(Tf, np.conj(Tf))
chi = (S_th[3, :] / S_th[0, :]).real                            # ANYON ordering [2.618,-1,-1,0.382]
R_true = S_th @ np.conj(T_th)

# Real error form (deposited, SORTED ordering) -> map into the anyon frame
delta_sorted = np.array(rA["eigenvalues_scaled"]) - np.array(rA["ideal_character"])
idx = np.argsort(chi)                                           # anyon positions, ascending
delta_anyon = np.empty(4); delta_anyon[idx] = delta_sorted
print("chi (anyon frame)   :", np.round(chi, 4))
print("delta (anyon frame) :", np.round(delta_anyon, 6),
      " | pair equality:", abs(delta_anyon[1] - delta_anyon[2]))

def synth_anyon(lam_anyon, seed=7, vec_err=0.0, vec_seed=101):
    """CONSISTENT: basis vector a carries the Wilson eigenvalue lam_anyon[a] AND
    the rotation matrix elements (S T*)_{a,:}. The eigenvector error rotates
    ONLY Mw (Q Mw Q^dag); the rotation itself stays untouched."""
    rw = np.random.default_rng(seed)
    W = np.linalg.qr(rw.normal(size=(4, 4)) + 1j * rw.normal(size=(4, 4)))[0]
    Mw = W @ np.diag(lam_anyon.astype(complex)) @ W.conj().T
    Mr = W @ R_true @ W.conj().T
    if vec_err > 0:
        rv = np.random.default_rng(vec_seed)
        K = rv.normal(size=(4, 4)) + 1j * rv.normal(size=(4, 4))
        K = (K - K.conj().T) / 2; K = K / np.linalg.norm(K)
        Q = expm(vec_err * K)
        Mw = Q @ Mw @ Q.conj().T
    return Mw, Mr

def evaluate(S_ex, T_ex):
    best = 1e9
    for perm in itertools.permutations(range(4)):
        P = np.eye(4)[list(perm)]
        Sp, Tp = P @ S_ex @ P.T, P @ T_ex @ P.T
        for conj in (False, True):
            Sc, Tc = (np.conj(Sp), np.conj(Tp)) if conj else (Sp, Tp)
            best = min(best, float(np.linalg.norm(Sc - S_th) + np.linalg.norm(Tc - T_th)))
    return best

def run(lam, seed, iscale, vec_err=0.0, n_starts=160):
    Mw, Mr = synth_anyon(lam, seed=seed, vec_err=vec_err)
    sols, blocks, _ = extract(Mw, Mr, n_starts=n_starts, input_scale=iscale)
    if not sols: return None, 0, blocks
    return min(evaluate(s["S"], s["T"]) for s in sols), len(sols), blocks

out = {"input_dev": INPUT_DEV, "delta_anyon": list(delta_anyon)}

# ---------------- GATE: baseline (strength 0) MUST find solutions ----------------
print("\n[GATE] baseline strength=0 (consistent construction):")
ok_gate = True
for seed in (7, 23):
    d, n, blk = run(chi.copy(), seed, 1e-13)
    print(f"  seed={seed}: {n} solutions, best={'%.2e' % d if d else 'NONE'}, blocks={blk}")
    ok_gate &= (d is not None and d < 1e-9)
out["gate_baseline"] = ok_gate
if not ok_gate:
    print(">>> GATE FAIL -- construction faulty again. ABORT, no physics finding.")
    json.dump(out, open(os.path.join(H, "p6_reconciliation_test.json"), "w"), indent=2)
    raise SystemExit(1)

# ---------------- PART 1: eigenvalue channel, REAL form ----------------
print("\n[PART 1] real error form on the EIGENVALUES, rotation exact:")
t1 = {}
for s in (0.3, 1.0, 3.0, 10.0, 30.0):
    vals = []
    for seed in (7, 23):
        d, n, _ = run(chi + s * delta_anyon, seed, max(1e-13, s * INPUT_DEV))
        vals.append(d)
        print(f"  {s:5.1f}x  seed={seed}: {n:2d} solutions, best={'%.3e' % d if d else 'NONE'}")
    t1[s] = [v for v in vals if v is not None]
out["part1_eigenvalue_channel"] = {str(k): v for k, v in t1.items()}

# ---------------- PART 2: eigenvector channel (the control that can fail) ----------------
print("\n[PART 2] eigenvector frame error (rotation exact, eigenvalues ideal):")
t2 = {}
for v in (1e-10, 1e-6, 1e-4, 3.3e-3):
    d, n, _ = run(chi.copy(), 7, max(1e-13, v), vec_err=v)
    t2[v] = d
    print(f"  vec_err={v:.1e}: {n:2d} solutions, best={'%.3e' % d if d else 'NONE'}")
out["part2_eigenvector_channel"] = {f"{k:.1e}": v for k, v in t2.items()}

# ---------------- PART 3: invariance demo (mechanism direct) ----------------
print("\n[PART 3] character residual under random B freedom (strength 1x):")
lam1 = chi + delta_anyon
Mw, Mr = synth_anyon(lam1, seed=7)
wv, E = np.linalg.eigh(Mw); R0 = E.conj().T @ Mr @ E
blocks = _blocks(wv, tol=max(1e-6, 3 * INPUT_DEV))
nB = sum(len(b) ** 2 for b in blocks)
lam_e = wv.astype(complex)
v_pos = int(np.argmax(wv))                                       # vacuum = largest eigenvalue
rng = np.random.default_rng(5)
chars, syms = [], []
for _ in range(300):
    B = _B_of(rng.uniform(-np.pi, np.pi, nB), blocks)
    A = B.conj().T @ R0 @ B
    w = lam_e * A[v_pos, :]
    best = min(
        float(np.linalg.norm(w - (np.vdot(u, w) / max(np.vdot(u, u).real, 1e-300)) * u) ** 2
              / max(float(np.linalg.norm(w) ** 2), 1e-300))
        for u in (A[xx, :] for xx in range(4)))
    chars.append(best)
    syms.append(float(np.linalg.norm(A - A.T)))                  # raw asymmetry proxy
chars = np.array(chars); syms = np.array(syms)
print(f"  character residual: min {chars.min():.6e}  max {chars.max():.6e}  "
      f"rel. scatter {(chars.max()-chars.min())/chars.mean():.2e}")
print(f"  asymmetry proxy   : min {syms.min():.3f}  max {syms.max():.3f}  (varies O(1))")
out["part3_invariance"] = {"char_min": float(chars.min()), "char_max": float(chars.max()),
                          "char_rel_spread": float((chars.max() - chars.min()) / chars.mean()),
                          "sym_min": float(syms.min()), "sym_max": float(syms.max())}

# ---------------- VERDICT ----------------
print("\n" + "=" * 78)
flat1 = all(len(v) > 0 and max(v) < 1e-9 for v in t1.values())
tracks2 = all(t2[k] is not None for k in t2) and t2[3.3e-3] > 100 * max(t2[1e-10], 1e-15)
inv3 = out["part3_invariance"]["char_rel_spread"] < 1e-6
verdict = ("A_STRUCTURALLY_PROJECTED" if (flat1 and tracks2 and inv3) else
           "B_SCALED" if not flat1 else "MIXED_SEE_RAWDATA")
print(f"part1 flat: {flat1} | part2 scales: {tracks2} | part3 invariant: {inv3}")
print(f">>> VERDICT: {verdict}")
out["verdict"] = verdict
json.dump(out, open(os.path.join(H, "p6_reconciliation_test.json"), "w"), indent=2)
print("EXPORT OK -> p6_reconciliation_test.json")
