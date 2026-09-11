#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The minimal exact-diagonalization path: lattice geometry and ground space, from scratch.

WHAT THIS CLOSES
The deposited matrices (p6_V_gs_2x2.npy and the Wilson/rotation matrices projected onto it) were
produced by a run whose generator is not part of this record -- p6_manifest_2x2.json lists
lw_torus.py as "original run artifact (path not part of this record)". Two things therefore could
not be checked from the record alone: which lattice configuration each ROW of V_gs corresponds to,
and whether the ground space itself is reproducible. This script supplies both. It is not the
production engine and does not try to be: it builds the honeycomb torus, the allowed sector and
the ground space, and nothing else.

WHY THE ROW LABELS MATTER
Three of the four deposited ground-state columns vanish on the same five rows. Those five rows are
the vacuum configuration and the four elementary plaquette boundaries -- but that identification
needs the edge numbering, which lived only in the generator. Without it a reader sees a pattern and
cannot say what it is about.

WHAT IT WRITES
p6_geometry_2x2.json -- for every row of V_gs, the edge configuration it stands for, plus the
plaquette rings and the vertex incidence that define the numbering. Written only under --write.

RUNNING IT

    python p6_ed_path.py              geometry + sector + ground space, verify, write nothing
    python p6_ed_path.py --write      the same, and write p6_geometry_2x2.json
    python p6_ed_path.py --fast       geometry + sector only, skip the diagonalization

Default is read-only on purpose: a diagnostic run that writes is not a diagnostic run.

This script imports no sibling module, so running it creates no __pycache__ in the record. It reads
p6_V_gs_2x2.npy only to verify against it.

CONVENTIONS (they are the point of this file, so they are stated, not implied)
Edges are numbered 3*(i*M + j) + d with d = 0,1,2 for the x, y and z bond of unit cell (i,j).
A configuration is an integer whose bit for edge q is (cfg >> (E-1-q)) & 1 -- the FIRST edge is the
MOST significant bit. Getting this backwards yields a valid-looking sector of the same size whose
rows mean something else, and nothing downstream fails loudly. This is why the check against the
deposited array below is not optional.
"""
import json
import os
import sys

import numpy as np

H = os.path.dirname(os.path.abspath(__file__))

PHI = (1.0 + 5.0 ** 0.5) / 2.0
INV_PHI = 1.0 / PHI
D2 = PHI + 2.0                      # the total quantum dimension squared
SQ_PHI = PHI ** 0.5
N_CELLS = M_CELLS = 2               # the 2x2 torus of the deposited set


def edge(i, j, d, N, M):
    """Edge id of bond d in {0:x, 1:y, 2:z} of cell (i, j)."""
    return 3 * ((i % N) * M + (j % M)) + d


def vertex_edges(N, M):
    """The three edges meeting at each vertex. A vertices first, then B."""
    NM = N * M
    inc = {}
    for v in range(2 * NM):
        if v < NM:
            i, j = divmod(v, M)
            inc[v] = [edge(i, j, 0, N, M), edge(i, j, 1, N, M), edge(i, j, 2, N, M)]
        else:
            i, j = divmod(v - NM, M)
            inc[v] = [edge(i, j, 0, N, M),
                      edge(i + 1, j, 1, N, M),
                      edge(i, j + 1, 2, N, M)]
    return inc


def plaquette(i, j, N, M):
    """The six edges and six vertices of hexagon (i, j), in cyclic order A-B-A-B-A-B.

    The ring is DERIVED, not tabulated: two consecutive vertices of the cycle share exactly one
    edge, and the incidence lists already fix which. Writing the six edge ids out by hand is the
    kind of formula that is wrong the second time it is copied -- and a wrong ring still produces
    six plausible edge ids, so it does not announce itself. Deriving them makes the construction
    self-checking: "exactly one shared edge" fails the moment the geometry is off.
    """
    NM = N * M
    ip, jm = (i + 1) % N, (j - 1) % M
    verts = [(i % N) * M + (j % M),                 # A(i, j)
             NM + (i % N) * M + (j % M),            # B(i, j)
             ip * M + (j % M),                      # A(i+1, j)
             NM + ip * M + jm,                      # B(i+1, j-1)
             ip * M + jm,                           # A(i+1, j-1)
             NM + (i % N) * M + jm]                 # B(i, j-1)
    inc = vertex_edges(N, M)
    edges = []
    for k in range(6):
        shared = set(inc[verts[k]]) & set(inc[verts[(k + 1) % 6]])
        if len(shared) != 1:
            raise RuntimeError("hexagon (%d,%d): vertices %d and %d share %d edges, expected 1"
                               % (i, j, verts[k], verts[(k + 1) % 6], len(shared)))
        edges.append(shared.pop())
    if len(set(edges)) != 6:
        raise RuntimeError("hexagon (%d,%d): ring is not six distinct edges" % (i, j))
    return edges, verts


def fusion_ok(a, b, c):
    """Fibonacci fusion: tau x tau = 1 + tau, so a triple is allowed unless exactly one leg
    carries tau. Same condition the vertex constraint enforces -- stated once, used twice."""
    return (a + b + c) != 1


def allowed_sector(N, M):
    """Every edge configuration satisfying the vertex constraint, as sorted integers."""
    E = 3 * N * M
    idx = np.arange(2 ** E, dtype=np.int64)
    bits = ((idx[:, None] >> (E - 1 - np.arange(E))[None, :]) & 1).astype(np.int8)
    keep = np.ones(len(idx), dtype=bool)
    for qs in vertex_edges(N, M).values():
        keep &= bits[:, qs].sum(axis=1) != 1
    return idx[keep], bits[keep]


def f_symbol(leg, old_prev, old_cur, new_cur, new_prev, F2):
    """The Levin-Wen F-symbol for one hexagon corner, doubled Fibonacci, string type tau.

    Four fusion nodes touch this corner and all four must be allowed -- the two that exist before
    the move and the two that exist after it. They are spelled out rather than packed into
    positional arguments: a six-argument symbol invites exactly the mix-up that produced
    E0 = -3.056 and GSD 1 here on the first attempt, and a wrong assignment still returns
    plausible numbers for every single corner. The ground-state energy is what catches it.
    """
    before_outer = fusion_ok(leg, old_prev, old_cur)      # leg, incoming old, current old
    before_inner = fusion_ok(1, new_cur, old_cur)         # the tau being fused in
    after_outer = fusion_ok(leg, new_prev, new_cur)       # leg, incoming new, current new
    after_inner = fusion_ok(old_prev, 1, new_prev)        # the tau, on the incoming side
    if not (before_outer and before_inner and after_outer and after_inner):
        return 0.0
    # Only the all-tau corner carries a nontrivial F-symbol; the index pair is
    # (current OLD, incoming NEW), not (current old, current new).
    if leg == 1 and old_prev == 1 and new_cur == 1:
        return F2[old_cur, new_prev]
    return 1.0


def plaquette_operator(sector, pos, ring, legs, E):
    """B_p^tau restricted to the allowed sector, built directly there.

    The full-Hilbert-space form would be a 2^E x 2^E array -- 134 MB per plaquette at E=12, of
    which all but 175x175 entries are discarded. Building inside the sector keeps the same operator
    and stays under 250 kB. `pos` maps a configuration integer to its row index.
    """
    F2 = np.array([[INV_PHI, 1.0 / SQ_PHI], [1.0 / SQ_PHI, -INV_PHI]])
    n = len(sector)
    Bt = np.zeros((n, n))
    shifts = [E - 1 - q for q in range(E)]
    for col, cfg_int in enumerate(sector):
        cfg = [(int(cfg_int) >> s) & 1 for s in shifts]
        old = [cfg[q] for q in ring]
        legl = [cfg[q] for q in legs]
        for new in range(64):
            nb = [(new >> k) & 1 for k in range(6)]
            amp = 1.0
            for k in range(6):
                amp *= f_symbol(legl[k], old[k - 1], old[k], nb[k], nb[(k - 1) % 6], F2)
                if amp == 0.0:
                    break
            if amp == 0.0:
                continue
            out = int(cfg_int)
            for k in range(6):
                bit = 1 << (E - 1 - ring[k])
                out = (out | bit) if nb[k] else (out & ~bit)
            row = pos.get(out)
            if row is None:          # would mean B_p leaves the sector -- it must not
                raise RuntimeError("B_p left the allowed sector at column %d" % col)
            Bt[row, col] += amp
    return Bt


def ground_space(N, M, verbose=True):
    """The doubled-Fibonacci ground space of the honeycomb torus, by exact diagonalization."""
    E = 3 * N * M
    sector, _ = allowed_sector(N, M)
    pos = {int(c): k for k, c in enumerate(sector)}
    n = len(sector)
    H_op = np.zeros((n, n))
    for i in range(N):
        for j in range(M):
            ring, verts = plaquette(i, j, N, M)
            ring_set = set(ring)
            legs = []
            inc = vertex_edges(N, M)
            for v in verts:
                rest = [q for q in inc[v] if q not in ring_set]
                if len(rest) != 1:
                    raise RuntimeError("plaquette (%d,%d) is not clean" % (i, j))
                legs.append(rest[0])
            Bt = plaquette_operator(sector, pos, ring, legs, E)
            H_op -= (np.eye(n) + PHI * Bt) / D2
    H_op = (H_op + H_op.T) / 2.0
    evals, evecs = np.linalg.eigh(H_op)
    E0 = float(evals[0])
    gsd = int(np.sum(evals <= E0 + 1e-8))
    gap = float(evals[gsd] - evals[gsd - 1])
    if verbose:
        print("  sector %d   E0 %.6f   GSD %d   gap %.4f" % (n, E0, gsd, gap))
    return sector, evecs[:, :gsd], dict(sector_dim=n, E0=E0, GSD=gsd, gap=gap)


def main():
    write = "--write" in sys.argv
    fast = "--fast" in sys.argv
    N, M = N_CELLS, M_CELLS
    E = 3 * N * M

    print("=" * 78)
    print("minimal ED path -- honeycomb torus %dx%d, doubled Fibonacci" % (N, M))
    print("=" * 78)

    sector, bits = allowed_sector(N, M)
    print("\n[1] geometry and allowed sector")
    print("    %d edges, %d vertices, %d plaquettes" % (E, 2 * N * M, N * M))
    print("    allowed configurations: %d of %d" % (len(sector), 2 ** E))

    rings = []
    for i in range(N):
        for j in range(M):
            ring, _ = plaquette(i, j, N, M)
            rings.append(sorted(ring))
    print("    plaquette rings: %s" % rings)

    # The five rows the deposited columns vanish on, named rather than merely located.
    pos = {int(c): k for k, c in enumerate(sector)}
    vac = pos.get(0)
    ring_rows = []
    for ring in rings:
        cfg = 0
        for q in ring:
            cfg |= 1 << (E - 1 - q)
        ring_rows.append(pos.get(cfg))
    print("\n[2] the distinguished rows, by construction")
    print("    vacuum (empty configuration) -> row %s" % vac)
    print("    plaquette boundaries         -> rows %s" % ring_rows)

    ok = True
    V = None
    gs = None
    info = {}
    dep = os.path.join(H, "p6_V_gs_2x2.npy")
    if os.path.isfile(dep):
        V = np.load(dep)
        rows = sorted([vac] + ring_rows)
        null_rows = sorted(np.where((np.abs(V) <= 1e-12).sum(axis=1) >= 3)[0].tolist())
        same = rows == null_rows
        print("\n[3] against the deposited array")
        print("    rows where three columns vanish: %s" % null_rows)
        print("    rows named above               : %s" % rows)
        print("    %s" % ("MATCH -- the numbering agrees with the deposited object"
                          if same else "MISMATCH -- the edge numbering does NOT agree"))
        ok &= same
    else:
        print("\n[3] p6_V_gs_2x2.npy not found -- verification skipped, stated not hidden")

    if not fast:
        print("\n[4] ground space by exact diagonalization")
        sector2, gs, info = ground_space(N, M)
        if V is not None:
            P1, P2 = gs @ gs.T, V @ V.T
            dev = float(np.abs(P1 - P2).max())
            print("    max|P_computed - P_deposited| = %.3e" % dev)
            print("    %s" % ("MATCH -- same ground space (the basis within it may differ)"
                              if dev < 1e-10 else "MISMATCH -- different subspace"))
            ok &= dev < 1e-10
    else:
        print("\n[4] skipped (--fast)")

    if write:
        rec = {
            "what": "row -> edge configuration for the rows of p6_V_gs_2x2.npy",
            "lattice": {"N": N, "M": M, "edges": E, "vertices": 2 * N * M,
                        "edge_numbering": "3*(i*M + j) + d, d = 0,1,2 for the x, y, z bond",
                        "bit_order": "bit of edge q is (cfg >> (E-1-q)) & 1; first edge is MSB"},
            "plaquette_rings": rings,
            "vertex_edges": {str(k): v for k, v in vertex_edges(N, M).items()},
            "row_to_configuration": [int(c) for c in sector],
            "distinguished_rows": {"vacuum": vac, "plaquette_boundaries": ring_rows},
            "ground_space": info,
            "generated_by": "p6_ed_path.py --write",
        }
        out = os.path.join(H, "p6_geometry_2x2.json")
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=2)
        print("\nwritten: %s (%d bytes)" % (os.path.basename(out), os.path.getsize(out)))
    else:
        print("\n(nothing written; pass --write to produce p6_geometry_2x2.json)")

    print("\n" + "=" * 78)
    print("all checks passed" if ok else "CHECKS FAILED -- see above")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
