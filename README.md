# Resolving the Full Modular Data of a Doubled-Fibonacci String-Net on a Finite Torus by Point-Group Rotation, and an Obstruction to Quasi-One-Dimensional Anyon Exchange

**Author:** Berkay Yüksel Sayim
**ORCID:** [0009-0004-4993-7352](https://orcid.org/0009-0004-4993-7352)
**DOI (this version):** [10.5281/zenodo.21362246](https://doi.org/10.5281/zenodo.21362246)

## Abstract

We report a real-space resolution of the full modular data of a
doubled-Fibonacci Levin–Wen string-net on a finite torus, together with a
bounded negative result for the dynamical half of the same program. (i) The
vacuum row of the modular S-matrix, |S_0a| = d_a/D, is confirmed by two
structurally distinct lattice routes, agreeing at <1e-9 (Verlinde–Casimir vs.
theory at 1.1e-16). (ii) GSD = 4 is a separate check, and |S| is 24×O(4)
gauge-invariant and size-robust to ≤2.5e-11. (iii) The Wilson-loop algebra
*alone* does not fix the chiral sign — measured: the loop algebra is
ℂ⊕M₃(ℂ) (dim 10 = 1²+3²), its antisymmetric part vanishes on the real-symmetric
holonomy, and the chiral splitting capacity is zero — a concrete instance of a
theorem of Li and Mong. Their named remedy, a point-group rotation of the
ground-state manifold (real but *not* symmetric, so it carries the phase the
loops cannot), resolves the full signed S and T from the lattice on **two**
sizes (3×3 and 2×2), up to the Z₂ relabeling intrinsic to an achiral theory;
a validated two-channel error budget with no free parameters accounts for the
residual on both. The modular data themselves are not new (Francuz–Dziarmaga,
2020); new here are the finite-lattice demonstration of the rotation channel,
the quantified error budget, and the measured loop-algebra limit as an instance
of the Li–Mong bound. Complementing this, the Hamiltonian transport of a single
mobile Ising anyon on the one-dimensionally mapped ladder cannot realize the
non-Abelian braid matrix: [U(C1),U(C2)] = 0, established by two independent
obstructions, three null escape routes, and a machinery-alive control.

## The evidence classes (read this first)

The classes are **not interchangeable**, and every entry in `p6_c1_modular.json`
and Table I of the paper is labeled with its class:

| Class | What it means | Where it applies |
|---|---|---|
| **(a)** | two structurally distinct **lattice** routes | the vacuum row `|S_0a| = d_a/D` **only** |
| **(b)** | the topological spins: loop overlaps give only algebraic consistency; the direct lattice measurement is class (d) | the loop route and its limit |
| **(c)** | a separate lock | GSD = 4 |
| **(d)** | the **point-group rotation channel** that resolves the full signed S and T | 3×3 and 2×2, `p6_c3_results.json` / `p6_lattice_extraction.json` / `p6_extraction_2x2.json` |

**The honest cut moved, it did not vanish.** An earlier framing said "the lattice does not
carry the chiral sign." That was true of the **Wilson-loop route** and false of
the **lattice**: the rotation channel — real but not symmetric — carries it.
The measured loop-algebra limit (dim 10, split 0) is now framed as a concrete
instance of the Li–Mong bound, and the rotation channel is that paper's named
remedy, not our discovery.

## Contents

| File | What it is |
|---|---|
| `main_v1.0.tex` / `main_v1.0.pdf` | the manuscript (revtex4-2) |
| `p6_c1_modular.json` | the class-(a)–(c) quantities and the measured loop-route limit, each with evidence class and per-entry provenance |
| `p6_c3_results.json` | the C₃ symmetry on 3×3 (orthogonality dev. 2.96e-15, order-3 dev. 4.35e-15) and the algebra jump 10→16; also the `structurally_forced` block (spectra that are forced, not measured). The 2×2 orthogonality dev. 4.41e-15 that the paper quotes alongside it lives in `p6_extraction_2x2.json`. |
| `p6_lattice_extraction.json` | the signed S,T extraction on 3×3: four solutions per computed rotation center, twelve entries over three centers, of which **eight are distinct** — two of the three centers fall in the same permutation class and give identical results (see `p6_center_independence.py`) |
| `p6_M_R_3x3_center20.npy` / `p6_M_R_3x3_center10.npy` | the projected C₃ rotation for centers (2,0) and (1,0), one per permutation class; `p6_M_R_3x3.npy` is center (0,0). The permutation depends on the center only through (j−i) mod 3 |
| `p6_center_independence.py` | checks from this record that the projected rotation does not depend on the center (agreement ~1e-15) and that the extractor returns the same block structure and the recorded residuals on each |
| `p6_extraction_2x2.json` | the signed S,T extraction on 2×2 (9 solutions / 3 centers), plus the two-channel budget and the Wfat_s0 vacuum control |
| `p6_reconciliation_test.json` | the two-channel error budget (transfer 2.081, quadratic k, first-order projection) |
| `p6_reconciliation_test.py` | the deterministic generator of the row above (all RNGs seeded); takes only `p6_loop_ideality.json` and `p6_extract_c3.py` as input |
| `p6_loop_ideality.json` | the input scale (loop ideality 2.2e-2 → 3.3e-3) |
| `p6_extract_c3.py` | the rotation-channel extractor (takes only measured matrices) |
| `p6_verify_patch.py` | the pre-registered referee harness for the extractor (five acceptance gates + scaling gate; synthetic test data, no engine import) |
| `p6_V_gs_3x3.npy` / `p6_Mgs_i_3x3.npy` / `p6_Mgs_j_3x3.npy` / `p6_M_R_3x3.npy` | the measured 3×3 input matrices (ground-space embedding, the two Wilson-loop matrices) and the C₃ rotation matrix of center (0,0), so the 3×3 extraction runs from this record |
| `p6_V_gs_2x2.npy` / `p6_Mgs_i_2x2.npy` / `p6_Mgs_j_2x2.npy` / `p6_M_R_2x2.npy` | the measured 2×2 input matrices (ground-space embedding, the two Wilson-loop matrices) and the C₃ rotation matrix of center (0,0) |
| `p6_manifest_2x2.json` | provenance of the 2×2 set: source file per export, checksums, shapes, and the reproduction recipe |
| `p6_nogo.json` | the C0 no-go: structural route, positive controls, machinery-alive control |
| `p6_manifold.json` | the N=3 dense-ED spectrum (the oracle's energetics) |
| `p6_oracle.json` | the algebra oracle (`source: algebra_oracle`) — explicitly **not** a measurement |
| `p6_modular_reference.py` | self-contained analytic reference (numpy only), 27 checks |
| `LICENSE` | CC BY 4.0 (paper, figures, data) |
| `LICENSE-CODE` | MIT (source code, `*.py`) |

## License
- Paper, figures, and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — see `LICENSE`
- Source code (`*.py`): [MIT License](https://opensource.org/licenses/MIT) — see `LICENSE-CODE`

## Reproducing

These scripts are meant to be run, not imported: none of them has a
`if __name__ == "__main__":` guard, so an `import` executes them. Two of them
write, and both write into their own directory: `p6_verify_patch.py` creates a
new file, `verify_p6_extract_c3.json`; `p6_reconciliation_test.py` overwrites
the deposited `p6_reconciliation_test.json`. If you want to keep the deposited
files untouched, work on a copy of the archive. `p6_extract_c3.py` and
`p6_center_independence.py` write nothing, so the import the snippets below
use, `from p6_extract_c3 import extract`, is safe; `p6_modular_reference.py`
also writes nothing, but importing it still runs its 27 checks.

The analytic reference runs standalone:

```
python p6_modular_reference.py
```

`numpy` only, deterministic; expected output `all 27 checks passed`. It rebuilds
the analytic doubled-Fibonacci modular data from scratch (no engine import, no
stored oracle) and checks every analytic target, including the Gauss sum
p₊ = D and the order-3 rejection at 5.236. The lattice extractions
(`p6_lattice_extraction.json`, `p6_extraction_2x2.json`) are the deposited outputs of
runs that require the exact-diagonalization engine; the extractor
`p6_extract_c3.py` takes only the measured Wilson-loop and rotation matrices as
input — no analytic S/T, no F/R, no solvability selector.

The two-channel error budget runs from this record as well:

```
python p6_reconciliation_test.py
```

All RNGs are seeded (7, 23, 101, 5), so the run is deterministic. It reads
`p6_loop_ideality.json`, imports `p6_extract_c3.py`, and reproduces every numeric
field of `p6_reconciliation_test.json`: the transfer factor 2.081
(2.08057/2.08052/2.08053/2.08088 across four and a half decades, scatter 0.018%)
and the quadratic eigenvalue channel. ⚠️ The script rewrites
`p6_reconciliation_test.json` in place, and the `interpretation` fields of the
deposited file are added by hand and are not regenerated — run it on a copy if you
want to keep them. It takes a few minutes.

On both lattice sizes the extraction step itself runs from this record.

### ⚠️ The extraction input is not the same file on both sizes

**3×3 uses `p6_Mgs_j_3x3.npy`, 2×2 uses `p6_Mgs_i_2x2.npy`.** The convention does not
carry over, and neither call fails loudly if you pick the wrong one — you get
plausible solutions of the same doubled-Fibonacci class with residuals that do
not match the deposit. The reason is historical: the deposited files are named
after their **source operators**, while the original run scripts named their
local variables independently of that, so on the 2×2 size the matrix deposited
as `p6_Mgs_i_2x2.npy` is the one the run script calls `Mj` internally. Check your
residuals against the values quoted below before drawing any conclusion.
`n_starts=240` is likewise part of both recipes: the `extract()` default of 400
returns a larger solution set.

**3×3:**

```python
import numpy as np
from p6_extract_c3 import extract
sols, blocks, _ = extract(np.load("p6_Mgs_j_3x3.npy"), np.load("p6_M_R_3x3.npy"),
                          n_starts=240, input_scale=3.2668e-3)
print(len(sols), blocks)     # -> 4 [[0, 1], [2], [3]]
```

These are the four solutions of rotation center `(0, 0)` recorded in
`p6_lattice_extraction.json`; their residuals, sign structure, and T deviations
reproduce the deposited entries (residuals 2.9426946078281404e-10,
3.262291310487419e-10, 3.628959392477765e-10, and 4.801729437638451e-10) to all
printed digits. The T diagonal carries e^(∓4πi/5), 1, 1. The input scale
3.2668e-3 is the deposited loop ideality from `p6_loop_ideality.json`, not a fitted
quantity.

The two further rotation centers `(1, 1)` and `(2, 0)` in that file are **not
independent repetitions.** On an L×L torus the nine possible centers generate
only three distinct edge permutations, and centers sharing a permutation run the
same computation: `(0, 0)` and `(1, 1)` fall in the same class and their four
solutions are identical entry for entry, so the twelve entries are eight
distinct solutions. Across all nine centers the projected rotation agrees to
about 1e-15. You can check this from the record without the lattice basis:

```
python p6_center_independence.py
```

It loads `p6_M_R_3x3.npy`, `p6_M_R_3x3_center20.npy` and
`p6_M_R_3x3_center10.npy` — one matrix per permutation class — and reports the
pairwise agreement, the solutions the extractor finds on each, and, for the two
centers that were part of the recorded run, the comparison against the residuals
stored in `p6_lattice_extraction.json` (takes a few minutes; it writes nothing).
The construction of these matrices from the string-net basis is not part of this
record.

**2×2:**

```python
import numpy as np
from p6_extract_c3 import extract
sols, _, _ = extract(np.load("p6_Mgs_i_2x2.npy"), np.load("p6_M_R_2x2.npy"),
                     n_starts=240, input_scale=0.022118297120066792)
print(len(sols))          # -> 3
```

These are the three solutions of rotation center `(0, 0)` recorded in
`p6_extraction_2x2.json`. Scored against the doubled-Fibonacci target they
reproduce the deposited residuals 7.24444564537034e-10, 8.575583316540904e-10,
and 1.7860418785978925e-09 to all printed digits. Feeding `p6_Mgs_j_2x2.npy`
instead returns three solutions of the same class with residuals near
1.07e-9, 2.34e-9, and 5.17e-9 — plausible, but not the deposited ones. The
remaining six of the nine entries in that file belong to the two further
rotation centers, whose rotation matrices the run script builds from the lattice
basis rather than reading them from a deposited file.

On both sizes the residual is not an output of `p6_extract_c3.py`. It is formed in
a separate scoring step against the analytic doubled-Fibonacci target, which
`p6_modular_reference.py` rebuilds from scratch; the driver script that
performed the original scoring is not part of this record.

## What the extractor does, precisely (it is released with the paper)

It solves a nonlinear least-squares problem whose residual encodes **generic
MTC axioms** — S symmetric, real nonnegative vacuum row, |θ_a|=1,
fusion-character consistency, S² = charge conjugation, and (ST)³ = ζS² with
|ζ|=1 (the central-charge phase, left free). These axioms **define the object
sought** and are applied identically to every input (anti-fudge). No analytic
S/T, no F/R symbol, no solvability selector is fed in (a poisoning test confirms
the theory helpers are never reached). Verlinde integrality and S-unitarity are
**output gates after the solve**; an independent referee harness rechecks the
output. The chiral sign comes from the measured rotation M_R; the modular axioms
only sieve out nonmodular branches.

## Two objects that must not be confused

- **M_ττ = −1/φ² = −0.382** is the **target value** of the dynamic half. It was
  a no-go and was **never measured**. Do not confuse it with the static Wilson
  eigenvalue **+1/φ²** (opposite sign).
- **E_0 = −4** is the **2×2** ground-state energy; on 3×3 it is **−9**. E_0 is
  extensive (one unit per plaquette); only the GSD and the gap are invariant.
- The **Z₂ relabeling** (τ↔τ̄, from complex conjugation in an achiral theory) is
  a different thing from the **fusion-phase freedom** (trivial for doubled
  Fibonacci because |𝒜|=1). The rotation channel *separates* the chiral pair but
  does not *label* it — for this method or any other.
