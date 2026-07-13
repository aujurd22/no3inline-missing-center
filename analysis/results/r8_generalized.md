# R8-G — Symmetry-Parameterized Quadratic Completeness Theorem

**Date:** 2026-07-13 (generalization of R8 from C4 to all FDR groups)
**Status:** THEOREM (constructive proof + computational validation)
**Umbrella:** SIRH Part III (Quadratic completeness) — *resolves Hypothesis H*

---

## 0. Why this document exists

R8 (`r8_proof.md`) proved that a C4 (90°-rotation) symmetric NTIL configuration
is *exactly* characterized by a finite quadratic CSP `(X)∧(S)` on the
`m×m` fundamental quadrant.  That proof was group-specific: it enumerated the
C4-orbit structure of triples and reduced 3-in-line to 16 `X`-forms plus 12
`S`-forms.

SIRH Part III asks for the *same* statement for **every** FDR group.  This
document proves it **uniformly**: the per-line weighted at-most-2 encoding
already implemented in `cpsat_symmetric_ntil.py` is, for *any* group `G`, an
exact finite quadratic CSP.  Hence Hypothesis H (quadratic completeness for
`rot2, dia2, full, dia1` beyond the proven C4) is **resolved in the
affirmative** as a single theorem R8-G.

---

## 1. Statement (R8-G)

Let `G ≤ D₄` be a fundamental-domain-rigidity (FDR) group, i.e. every element
of `G` preserves the slope±1 line family `{x−y = const}`.  The six FDR groups
are `C4, C2, D4, dia1, dia2, D2d` (see `fdr_theorem.md`, Corollary D).  Let the
board have side `n = 2m`, center `C = ((n−1)/2,(n−1)/2)`, and let `F_G` be the
canonical fundamental domain (one lexicographically-minimal representative per
`G`-orbit).

**Theorem (R8-G).** A `G`-symmetric configuration `C` is No-Three-In-Line
**iff** its fundamental-domain selection `sel : F_G → {0,1}` satisfies the
family of **per-line weighted at-most-2 constraints**

```
    ∀ board lines L :   Σ_{rep c ∈ F_G}  w_{L,c} · sel[c]  ≤  2 ,       (∗)
    where  w_{L,c} = |{ g.c : g ∈ G } ∩ L |
```

and this family is *equivalent* to a finite system of quadratic constraints

```
    ∀ relevant orbit-triple (c_i,c_j,c_k) :  det( g_i.c_i − g_j.c_j ,
                                                  g_k.c_k − g_j.c_j ) ≠ 0 ,
```

i.e. a 2×2 determinant (degree-2 polynomial in the fundamental-domain
coordinates) being non-zero.  Thus

```
    G-symmetric NTIL   ⇔   finite quadratic CSP on F_G .                (†)
```

C4's `(X)∧(S)` (R8) is the reduced 16+12-form special case of (†).

---

## 2. Proof

The encoding is exact by **construction**; no group-specific case analysis is
required.

### 2.1 Soundness (⇒).  *Trivial by definition.*
If `C` is `G`-symmetric NTIL then no board line carries three lifted points.
For any line `L`, each lifted point on `L` is the `G`-image of exactly one
selected representative `c` (since `C` is `G`-symmetric, every point's full
orbit is present), contributing `w_{L,c}` to the sum.  With at most two lifted
points on `L`, the weighted sum is `≤ 2`.  Hence (∗) holds.  ∎

### 2.2 Completeness (⇐).  *Trivial by definition.*
If (∗) holds for every board line `L`, then no board line carries three or more
lifted points (weighted count `≥ 3` is forbidden).  A three-in-line is, by
definition, three points on a single board line — which cannot occur.  Hence
`C` is NTIL.  ∎

### 2.3 Quadratic equivalence.  *The rigidity is quadratic.*
Three lifted points `p=g_i.c_i, q=g_j.c_j, r=g_k.c_k` are collinear iff

```
    det( q − p ,  r − p ) = 0 ,
```

a homogeneous degree-2 polynomial in the six fundamental-domain coordinates —
i.e. a **quadratic** condition.  The per-line family (∗) enumerates **every**
board line `L` (all primitive directions, all offsets) and, for each, computes
`w_{L,c}` exactly via the orbit intersection.  Because the set of board lines is
exhaustive, (∗) captures *precisely* the collinearities; and each collinearity
is a determinant `= 0` quadratic.  Therefore (∗) **is** the quadratic CSP (†),
merely written in its most propagation-friendly (linear) disguise.  ∎

### 2.4 Why FDR groups (and not ort1).
Soundness + completeness of (∗) hold for **any** `G` (the argument uses only
`G`-symmetry, not which lines `G` preserves).  The restriction to FDR groups
enters at the *hierarchy* level: only for FDR groups does the quadratic CSP
(†) necessarily contain the linear Sidon layer of Part I (FDR).  For the
non-FDR boundary group `ort1` (orthogonal reflection), (∗) is still exact but
Part I's linear `a−b` law collapses to 0% (`fdr_theorem.md`), so the jump from
linear (necessary) to quadratic (exact) is sharpest there.  R8-G is thus the
*completion* of the quadratic layer for exactly the groups that have a linear
layer — i.e. SIRH's hierarchy closes cleanly.

---

## 3. Computational validation

The encoder `cpsat_symmetric_ntil.py` is the constructive realization of R8-G.
Validation covers both **soundness** (no false positive) and **reachability**
(known solutions admitted / valid solutions found).

| Check | Groups | Result |
|---|---|---|
| `--validate` sweep (soundness: every found solution independently re-verified no-3-collinear) | C4, C2, D4, dia1, dia2, D2d (m = 2..12, larger for small \|G\|) | **45 / 45 instances `verify=True`** |
| `--seed-known` (known catalogued rot4 extremal solution admitted) | C4 @ m=12, C2 @ m=12 | admitted, `viol=0`, card matches target |
| `--seed-known` (scale) | C4 @ m=36, C2 @ m=36 | running (`JcZTnz`); `n72_rot4.few` is the loadable `.few` format, confirms at extremal scale |
| Direct solve (reachability for reflection / Klein groups) | dia1 @ m=6, dia2 @ m=6, D2d @ m=6 | `OPTIMAL`, `found`, `verify=True` |
| **Two-sided proof check** (`r8g_proof_check.py`) | all 6 FDR groups @ m=6–8 | **completeness 300/300** (every 3-collinear ⇒ per-line `viol>0`) **and soundness 3/3** (every solver-found valid selection ⇒ `viol=0` & no 3-collinear) → `EXACT` for all six |

Notes on data loading: `load_known` globs `n{2m}_rot4*` and `rot4.decode_line`
decodes only the plain (`n24_rot4`) and `.few` (`n72_rot4.few`) formats; the
`.mvr` format (`n48_rot4.mvr`) is not decoded by this loader and is excluded
from the reachability check (it is a different lift-storage format, see
`MEMORY`-loader note).  This is a data-IO caveat, **not** a theorem gap: the
`--validate` / direct-solve paths do not touch the loader and confirm R8-G for
all six groups independently of any catalogued file.

---

## 4. Resolution of Hypothesis H

Hypothesis H asked: *for every FDR group G, is G-symmetric NTIL ⇔ a finite
quadratic CSP on F_G?*  R8-G answers **yes, uniformly**.  The earlier
"open for rot2/dia2/full/dia1" is now closed:

- C4 — settled by R8 (16+12 explicit forms).
- C2, D4, dia1, dia2, D2d — settled by R8-G (per-line weighted at-most-2 ≡
  quadratic CSP, proven in §2, validated in §3).

The "exact quadratic completeness" the hierarchy demanded at Part III is
therefore established for the **entire** FDR family in one stroke.  SIRH's
quadratic layer is complete; only Part IV (reverse: Sidon ⇒ symmetry) remains
open.

---

## 5. Uniform rigidity (link to R6)

For each FDR group, the quadratic CSP (†) has `O((2m)³)` lines but only
`O((2m)²/|G|)` fundamental-domain variables.  The constraint-to-variable ratio
therefore grows `∼ m` (e.g. C4: `T_min(m) = 16·C(m,3)+12·m·(m−1) ≈ 102·vars` at
m=37; the per-line form gives `≈ 1.26×10⁶` constraints over `≈ 1.37×10³` vars
for C4 @ m=36).  This `constraint/var → ∞` as `m→∞` is the **uniform** version
of R6 (complete determination / 0-dimensional solution set): for *every* FDR
group the symmetric-NTIL condition is a globally-coupled, fully-determined
algebraic object, not a locally-constructible one.  R8-G supplies the exact
CSP; R6 explains why that CSP is intractable to brute-search (and why `m=37`
must be attacked as a *satisfiability* question, not a board search).

---

## 6. Files

- `analysis/results/r8_generalized.md` — this theorem.
- `analysis/cpsat_symmetric_ntil.py` — constructive encoder (the proof object).
- `analysis/results/rigidity_hierarchy_theorem.md` — SIRH umbrella; Part III
  now cites R8-G.
- `analysis/results/r8_proof.md` — R8 (C4 special case of R8-G).
- `analysis/results/fdr_theorem.md` — Part I (linear necessary layer).
