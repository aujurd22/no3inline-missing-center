# Research Status — no3inline-rigidity

> **What this repo is.** A computational notebook / research diary for an
> AI-assisted exploration of *No-Three-In-Line* that has grown into a theory of
> **symmetry-induced rigidity in extremal lattice configurations**.  NTIL is the
> testbed; the real subject is the two-layer rigidity structure (linear FDR →
> quadratic R8).  This file is a **navigation map only** — it marks what is
> proven, what is computational evidence, what is conjectured, and what was
> revised/abandoned.  It does **not** reorganise the content.
>
> See `two_layer_rigidity.md` for the unifying framework, and
> `README.md` for the full main-axis spine.

---

## ✅ Proven (Theorem)

| Result | Statement | Doc |
|--------|-----------|-----|
| **R8** (Quadratic-Sidon Completeness) | rot4 NTIL ⇔ `(X)∧(S)` quadratic CSP on an m-subset of the m×m fundamental quadrant | `analysis/results/r8_proof.md`, `quadratic_sidon_completeness.md` |
| **R7** (Quadratic Gap) | cross-quadrant collinearity is an irreducible *quadratic* condition; linear Sidon forms are insufficient | `analysis/results/quadratic_gap_theorem.md` |
| **FDR** (Fundamental Domain Rigidity) | slope±1-line-preserving symmetry `G≤D₄` ⇒ a−b Sidon law on `F_G` (ort1/iden are the boundary) | `analysis/results/fdr_theorem.md` |
| Th-56 | rot4 C4-orbit intercept 4-fold ⇒ `count(d)+count(−d) ≤ 2` (base case of FDR) | `fdr_theorem.md` |
| Th-39 | dia1 ⇔ 2-factor pseudograph reduction | (in README §2) |
| Th-42 | three-body resonance `E = 0` | (in README §2) |
| Th-48 | diagonal points arise from a cycle lattice | (in README §2) |
| Th-50 / Th-51 | direction mutual-exclusion (necessary); forbidden slopes `{0,∞,−1}`, slope 1 ⇔ diagonal | (in README §2) |
| Th-54 / Th-55 | displacement 2-cover equivalence; row/column safety | (in README §2) |
| Th-57 / Th-58 | cross-resonance criterion (necessary); s=1 Sidon (necessary) | (in README §2) |
| Container foundation | `2n` configurations = independent sets of the danger hypergraph `Hₙ`; codegree exactly `n−2` | `hypergraph_framework.md` (hard result) |

---

## 🔬 Strong computational evidence (not yet proven)

| Observation | Support |
|------------|---------|
| Structural invariants scale smoothly | source ratio ~0.26 stable (m=7–28 mean 0.2647, σ<0.005); multi-cycle norm; longest cycle ~0.79·m | `structural_scaling_2026-07-12.md` |
| m=37 sits inside the satisfiable regime | no phase transition at m=36→37 (structural + constraint-density continuous) | `m37_satisfiability_window.md` |
| **All m=3..36 have rot4 solutions** | corrected from an earlier false "gap" (was a `.mvr` loader bug) |
| FDR empirical corroboration | rot4/rot2/rct4/dia2/full/dia1 a−b 100%; ort1/iden 0% | `fdr_theorem.md` §empirical |
| R8 computational check | 93 known solutions pass `(X+S)`; 7,500 random templates `miss=0, false=0` | `quadratic_sidon_completeness.md` |

---

## ❓ Conjectures / Open problems

- **m=37 rot4 existence** — the frontier.  ⇔ satisfiability of the R8 quadratic
  CSP.  Exact CP-SAT re-encoding: **1,264,378** per-line at-most-2 constraints
  over 1,369 binary vars (`analysis/results/cpsat_encoding.md`).  OR-Tools
  attack launched (background task `fkB7gu`, 1800 s budget) — result pending.
  **OPEN.**
- **Guy–Kelly D(n)=2n** infinite family — holds for *some* n≤72, not all.
- **C4 necessity** ⇔ no collinearity 2-factor (Th-39 direction).
- **α lower bound ~1.5n** (Hall 1975) for general NTIL.
- **Extremal-config taxonomy** — beyond symmetric families: are there
  asymmetric / sporadic structure types? (left open by the SE revision)
- **Container ⇒ symmetry bridge** — does a high-min-degree container force a
  non-trivial symmetry `G` (hence fall under FDR / Layer 1)?  Unproved.

---

## 🔄 Revised / Abandoned

- **SE (Symmetry Emergence)** — *demoted.* "iden extinct at high n" was a
  sampling artefact: iden is 99.4% of solutions at n=20 and simply has no
  enumeration data for n≥28 (search infeasible, not extinction).  No longer a
  core assumption.
- **"rot4 solutions are a permutation matrix"** — *revised 2026-07-12.*
  The m fundamental-domain cells are an **m-subset** (rows/cols may repeat);
  verified on the full known-solution corpus.  R8's core theorem is untouched
  (it only needs m distinct orbit representatives).
- **"m=23..26 have no rot4 solutions"** — *revised.*  Was a `.mvr` loader bug;
  all m=3..36 have solutions.
- **"container method proves structured"** — *rephrased.*  Container methods
  give a *concentration framework* (few containers hold all configs), **not** a
  proof of symmetry/rigidity of contained configs.  (Per external review.)

---

## 📌 How the pieces connect

```
symmetry G ≤ D₄
   → Layer 1 LINEAR rigidity (FDR, THEOREM): slope±1 preservation ⇒ Sidon
   → R7 gap: linear insufficient (cross-quadrant collinearity is quadratic)
   → Layer 2 GEOMETRIC rigidity (R8, THEOREM): (X)∧(S) quadratic CSP ⇔ NTIL
   → R6: Layer-2 system is 0-dimensional (over-determined)
   → m=37 existence = Layer-2 CP-SAT (1.26M per-line constraints)  [OPEN, attack running]
```
See `two_layer_rigidity.md` for the full synthesis.
