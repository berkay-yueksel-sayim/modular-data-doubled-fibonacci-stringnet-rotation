#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The scoring driver: forms the reference error from this record, on both sizes.

WHAT THIS CLOSES
The paper states that the residual is "not an output of the solver" but is formed
in a separate scoring step against the analytic doubled-Fibonacci target. The
extractor (p6_extract_c3.py) and the analytic reference (p6_modular_reference.py)
are both deposited; the driver that joins them was not. This is that driver.

WHAT THE REFERENCE ERROR IS
For one extracted solution (S_ex, T_ex) and the analytic target (S_th, T_th):

    d = min over the 24 label permutations P and complex conjugation of
        || P S_ex P^T - S_th ||_F  +  || P T_ex P^T - T_th ||_F

The minimisation over permutation and conjugation is not a fit: it removes the
label freedom and the Z2 relabeling (tau <-> taubar) intrinsic to an achiral
theory, neither of which the extraction can fix. No analytic S or T enters the
extraction itself -- the theory values are used here, after the solve, and only
to measure a distance. compare_to_theory() in p6_extract_c3.py implements it.

RUNNING IT

    python p6_scoring_driver.py            both sizes, 240 starts   (~2.5 min)
    python p6_scoring_driver.py --quick     both sizes,  40 starts   (smoke test)

--quick finds FEWER solutions than the deposit records, and says so; it checks
that the ones it finds match, not that all of them were found. Use the full run
to reproduce the record.

This script writes nothing -- including no bytecode cache: it disables .pyc
writing before importing its sibling module, so running it leaves the record
byte-identical. It reads only files from its own directory.

THE INPUT CONVENTION DIFFERS BETWEEN THE TWO SIZES
3x3 scores from p6_Mgs_j_3x3.npy, 2x2 from p6_Mgs_i_2x2.npy. This is recorded in
p6_extraction_2x2.json as wilson_choice = "ops_wfat_i[1]" and is explained in the
README. Picking the other matrix does not fail loudly: it yields plausible
solutions of the same doubled-Fibonacci class whose residuals do not match the
deposit. The sizes are therefore kept in an explicit table below rather than
derived from a naming rule.
"""
import json
import os
import sys

# Importing a sibling module writes a bytecode cache (__pycache__/) next to it.
# In a deposited record that is an unlisted build artefact, so it is switched off
# before the first import -- which is why this line sits above them rather than
# with the other module-level settings.
sys.dont_write_bytecode = True

import numpy as np

from p6_extract_c3 import extract, theory_doubled_fib, compare_to_theory

H = os.path.dirname(os.path.abspath(__file__))

# One row per lattice size. The Wilson matrix differs by size on purpose (see the
# module docstring); the input scale is the deposited loop ideality of that run.
SIZES = [
    {"tag": "2x2",
     "wilson": "p6_Mgs_i_2x2.npy",
     "rotation": "p6_M_R_2x2.npy",
     "scale": 0.022118297120066792,     # p6_extraction_2x2.json: ideality_2x2
     "record": "p6_extraction_2x2.json",
     "center": "(0, 0)"},
    {"tag": "3x3",
     "wilson": "p6_Mgs_j_3x3.npy",
     "rotation": "p6_M_R_3x3.npy",
     "scale": 0.0032668,                # p6_loop_ideality.json: the 3x3 run
     "record": "p6_lattice_extraction.json",
     "center": "(0, 0)"},
]

FULL_STARTS = 240      # same as the README extraction snippet
QUICK_STARTS = 40


def recorded_residuals(path, center):
    """Every 'residual' recorded under one center key, in file order.

    The two records nest differently, so this walks the parsed structure rather
    than assuming a path. Reading them with a regex would also match residuals
    belonging to a neighbouring center.
    """
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    found = []

    def walk(node, inside):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "residual" and inside:
                    found.append(float(value))
                walk(value, inside or key == center)
        elif isinstance(node, list):
            for value in node:
                walk(value, inside)

    walk(doc, False)
    return sorted(found)


def score_one(size, n_starts, S_th, T_th):
    """Extract, then score against the analytic target. Returns sorted residuals."""
    wilson = np.load(os.path.join(H, size["wilson"]))
    rotation = np.load(os.path.join(H, size["rotation"]))
    solutions, blocks, _ = extract(wilson, rotation,
                                   n_starts=n_starts, input_scale=size["scale"])
    residuals = sorted(compare_to_theory(s["S"], s["T"], S_th, T_th)[0]
                       for s in solutions)
    return residuals, blocks


def main():
    quick = "--quick" in sys.argv
    n_starts = QUICK_STARTS if quick else FULL_STARTS

    print("=" * 78)
    print("scoring driver -- the reference error, formed from this record")
    print("=" * 78)
    print("mode: %s (%d starts)%s"
          % ("quick" if quick else "full", n_starts,
             "   NOTE: finds fewer solutions than the record" if quick else ""))

    S_th, T_th = theory_doubled_fib()
    print("analytic target rebuilt by theory_doubled_fib(); "
          "no theory value enters the extraction itself\n")

    problems = []
    for size in SIZES:
        residuals, blocks = score_one(size, n_starts, S_th, T_th)
        deposited = recorded_residuals(os.path.join(H, size["record"]),
                                       size["center"])
        print("%s   center %s   blocks %s" % (size["tag"], size["center"], blocks))
        print("   wilson input : %s" % size["wilson"])
        print("   scored       : %s" % ["%.9e" % r for r in residuals])
        print("   deposited    : %s" % ["%.9e" % r for r in deposited])

        if len(residuals) == len(deposited):
            worst = max(abs(a - b) for a, b in zip(residuals, deposited))
            ok = worst < 1e-18
            print("   max deviation: %.3e   %s" % (worst, "MATCH" if ok else "MISMATCH"))
            if not ok:
                problems.append("%s: deviation %.3e" % (size["tag"], worst))
        elif quick and len(residuals) < len(deposited):
            # Expected in quick mode: fewer starts reach fewer basins. Check that
            # the ones we did find are among the deposited values.
            unmatched = [r for r in residuals
                         if not any(abs(r - d) < 1e-18 for d in deposited)]
            print("   found %d of %d solutions (quick mode); unmatched: %d"
                  % (len(residuals), len(deposited), len(unmatched)))
            if unmatched:
                problems.append("%s: %d residuals not in the record"
                                % (size["tag"], len(unmatched)))
        else:
            print("   solution count %d does not match the record's %d"
                  % (len(residuals), len(deposited)))
            problems.append("%s: %d solutions vs %d recorded"
                            % (size["tag"], len(residuals), len(deposited)))
        print()

    print("=" * 78)
    if problems:
        print("%d problem(s):" % len(problems))
        for p in problems:
            print("   %s" % p)
    elif quick:
        # A summary must not claim more than the run checked. The quick mode
        # reaches fewer basins, so it can only report that what it found agrees.
        print("quick mode: every residual found matches the record, but not all")
        print("recorded solutions were reached. Run without --quick to reproduce")
        print("the deposited residuals in full.")
    else:
        print("the reference error reproduces the deposited residuals on both sizes")
    print("=" * 78)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
