"""
p6_center_independence.py -- the projected C3 rotation does not depend on the
rotation center, and the extra centers are therefore not extra evidence.

Why this exists. The 3x3 extraction was run for three rotation centers, and an
earlier version of this manuscript counted that as three repetitions. It is not
three repetitions. On an L x L torus the nine possible centers generate only
THREE distinct edge permutations: the permutation depends on the center only
through (j - i) mod 3. Centers that share a permutation perform literally the
same computation, so their projected rotations are identical bit for bit -- two
of the three centers used in the run fall in the same class, which is why two of
the twelve recorded solution sets coincide entry for entry and only eight
solutions are distinct. Across all nine centers the projected rotation agrees to
about 1e-15, the floating-point floor.

Stated plainly: running further centers is a consistency check on the optimizer,
not an independent repetition of the measurement. The independent repetition in
this work is the second lattice size.

The three deposited matrices are one representative per permutation class:
  p6_M_R_3x3.npy            center (0,0)   class (j-i) mod 3 = 0
  p6_M_R_3x3_center20.npy   center (2,0)   class 1
  p6_M_R_3x3_center10.npy   center (1,0)   class 2
The first two are the matrices behind the recorded runs, so the residuals this
script prints for them can be compared against p6_lattice_extraction.json
directly. Class 2 was not part of the recorded run; it is deposited so that all
three classes can be checked. The construction of these matrices from the
string-net basis is not part of this record (that is the exact-diagonalization
engine, a separate and heavier computation).

Load the deposited matrices, do not recompute them. A rotation matrix rebuilt
from the basis can differ in the last bit from one environment to the next, and
the multi-start search may then resolve a different number of minima; we have
seen four and three on inputs agreeing to 1e-15. This is not a tolerance effect:
with the deposited input scale the extractor's noise estimate eps is 3.2668e-3
for all three matrices, identical to the last digit, because the unitarity
deviation of the input (about 3e-15) never reaches it. It is the optimizer
landscape itself. Residual values and block structure are unaffected; only the
number of distinct minima the restart search separates is. We record this as an
honest limit of reproducibility rather than tune it away.

Deterministic within one run. Read-only: this script writes nothing.
Run:  python p6_center_independence.py     (a few minutes)
"""
import json
import os

import numpy as np

from p6_extract_c3 import extract, theory_doubled_fib, compare_to_theory

H = os.path.dirname(os.path.abspath(__file__))

CLASSES = [
    ("center (0,0), class 0", "p6_M_R_3x3.npy", "(0, 0)"),
    ("center (2,0), class 1", "p6_M_R_3x3_center20.npy", "(2, 0)"),
    ("center (1,0), class 2", "p6_M_R_3x3_center10.npy", None),
]
INPUT_SCALE = 3.2668e-3          # deposited loop ideality of the 3x3 run
N_STARTS = 240                   # same as the README extraction snippet

mats = [(lab, np.load(os.path.join(H, fn)), key) for lab, fn, key in CLASSES]
Mj = np.load(os.path.join(H, "p6_Mgs_j_3x3.npy"))

print("=" * 74)
print("1. the three deposited rotation matrices, pairwise")
print("=" * 74)
worst = 0.0
for a in range(len(mats)):
    for b in range(a + 1, len(mats)):
        d = float(np.abs(mats[a][1] - mats[b][1]).max())
        worst = max(worst, d)
        print(f"  max|{mats[a][0]} - {mats[b][0]}| = {d:.4e}")
print(f"\n  worst pairwise deviation: {worst:.4e}")
print("  -> the projected rotation is center independent at the numerical floor")

print()
print("=" * 74)
print("2. the extractor on each class, scored against the analytic target")
print("=" * 74)
S_th, T_th = theory_doubled_fib()
dep = json.load(open(os.path.join(H, "p6_lattice_extraction.json"), encoding="utf-8"))


def recorded_residuals(block_key):
    """residuals recorded for one center in p6_lattice_extraction.json"""
    found = []

    def walk(o, inside):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == "residual" and inside:
                    found.append(float(v))
                else:
                    walk(v, inside or k == block_key)
        elif isinstance(o, list):
            for v in o:
                walk(v, inside)

    walk(dep, False)
    return sorted(found)


ok_struct = []
for lab, M, key in mats:
    sols, blocks, _ = extract(Mj, M, n_starts=N_STARTS, input_scale=INPUT_SCALE)
    res = sorted(compare_to_theory(s["S"], s["T"], S_th, T_th)[0] for s in sols)
    tdiag = np.round(np.diag(sols[0]["T"]), 6)
    ok_struct.append((len(sols), str(blocks)))
    print(f"\n  {lab}")
    print(f"    {len(sols)} solutions, blocks {blocks}")
    print(f"    T diagonal: {tdiag}")
    print(f"    residuals:  {['%.6e' % r for r in res]}")
    if key is not None:
        rec = recorded_residuals(key)
        if len(rec) == len(res):
            dmax = max(abs(a - b) for a, b in zip(res, rec))
            print(f"    recorded in p6_lattice_extraction.json for {key}:")
            print(f"                {['%.6e' % r for r in rec]}")
            print(f"    max deviation from the record: {dmax:.3e}")
        else:
            print(f"    (record for {key} has {len(rec)} entries, run has {len(res)})")
    else:
        print("    this class is not part of the recorded run")

same = len(set(ok_struct)) == 1
print()
print("=" * 74)
print(f"  identical solution count and block structure on all three: {same}")
if same and worst < 1e-12:
    print("\nRESULT: center independence confirmed from the deposited files.")
    print("Further rotation centers are an optimizer consistency check, not an")
    print("independent repetition; the independent repetition is the 2x2 lattice.")
else:
    print("\nRESULT: deviation - inspect the numbers above.")
print("=" * 74)
