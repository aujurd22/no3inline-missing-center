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
| **★ SIRH** (Symmetry-Induced Rigidity Hierarchy) | **The umbrella theorem** compressing FDR/R7/R8.  Part I (FDR, linear, *necessary*) + Part II (why the quadratic layer is *forced* — proven group-independently; strict insufficiency **proven for C4 (R7) + rot2** (constructive witness)) + Part III (exact quadratic completeness: **proven for ALL six FDR groups = R8-G**, *resolves Hypothesis H*) + **Part IV (reverse: RESOLVED negative — Sidon is one-way necessary, not a symmetry witness)**. | `analysis/results/rigidity_hierarchy_theorem.md`, `analysis/results/part4_reverse.md` |
| **R8** (Quadratic-Sidon Completeness) | rot4 NTIL ⇔ `(X)∧(S)` quadratic CSP on an m-subset of the m×m fundamental quadrant | `analysis/results/r8_proof.md`, `quadratic_sidon_completeness.md` |
| **R8-G** (Generalized Quadratic Completeness) | **For every FDR group** `G∈{C4,C2,D4,dia1,dia2,D2d}`, `G`-symmetric NTIL ⇔ per-line weighted at-most-2 ⇔ finite quadratic CSP (det≠0) on `F_G`.  *Resolves Hypothesis H.*  Constructive proof (`cpsat_symmetric_ntil.py`); 45/45 validate `verify=True` + seed-known C4/C2 reachability. | `analysis/results/r8_generalized.md` |
| **R9** (mod-p descent + 2-regular reformulation) | **R9a:** rot4 NTIL ⇒ (X)+(S) over F_p **only for p>2m** — the natural "mod-m descent for prime m" is **FALSE** (two distinct lifted points can differ by m and coincide mod m, spuriously forcing det≡0 and rejecting *all* valid solutions).  **R9b:** rot4 NTIL ⇔ a permutation π on the m odd vertices + (X)+(S); the current cell-subset model searches a superset — the permutation/assignment CP-SAT model enforces 2-regularity for free and is the structurally correct, faster formulation. | `analysis/results/r9_modp_descent.md` |
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
| **Costas C1–C5** (FDR-transfer) | A Costas array can have exactly **6** D₄ symmetry types; `H`/`V` reflection **impossible** (ort1-boundary analogue); full D₄ impossible; `C4` needs `n ≡ 0,1 (mod 4)`.  Empirically confirmed to n=13 (`D4(full)=0`; `C2=C4=D2=0`). | `analysis/results/costas_symmetry_theorem.md` |
| **Costas C6** (★ BLOCKBUSTER) | **No C4-rotational Costas array exists for any n>1.**  A C4 4-orbit is a square = parallelogram ⇒ two distinct ordered pairs share a displacement ⇒ violates Costas.  Rigorously verified: proof + edge-case stress test + C4 CP-SAT sweep INFEASIBLE to n=32. | `costas_symmetry_theorem.md` §8, §8.1 |
| **Costas R** (★ rotational classification) | **No C2- or C4-symmetric Costas array exists for n ≥ 3.**  Rotational Costas ⇔ n ≤ 2 (the two order-2 arrays, both D2).  Pure-C2 type empty at all orders; D2 only at n=2.  Closes the "do rotational Costas exist?" open problem.  Verified by C2-enumeration (n=2→2, n=4..10→0). | `costas_symmetry_theorem.md` §9 |
| **Costas C** (capstone) | **Complete D₄ symmetry classification.** For n≥3 a Costas array's symmetry `G ≤ D₄` is in `{id, D, AD}` (no rotation, no H/V).  Only diagonal involutions D/AD occur (from order 5); all rotations and axis reflections impossible beyond n≤2.  Full classification resolved. | `costas_symmetry_theorem.md` §10 |
| **LQ** (linear/quadratic classification) | "no-collapse" conditions split into LINEAR (equality of linear forms, e.g. Costas displacement-distinct) vs QUADRATIC (vanishing 2×2 determinant, e.g. NTIL no-3-collinear = R7).  Explains why finite-field constructions work for Costas but fail for rot4-NTIL. | `sidon_costas_unification.md` §2 |
| **R8-C** (Costas analogue of R8) | C4-symmetric Costas of order 4m ⇔ *linear* CSP on m fundamental-quadrant cells (vs R8's *quadratic* CSP for NTIL).  Prescribes a CP-SAT attack on orders 32/33 orthogonal to finite fields. | `sidon_costas_unification.md` §4 |

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

- **★ Hypothesis H (Quadratic Rigidity Theorem) — RESOLVED (2026-07-13, R8-G).**
  For *every* FDR group `G`, `G`-symmetric NTIL ⇔ a finite quadratic CSP
  (det≠0 over orbit pairs in `F_G`).  Proven uniformly by the per-line
  weighted at-most-2 encoding (R8-G), validated computationally (45/45
  `verify=True` across all six groups; seed-known C4/C2 reachability).  Strict
  *insufficiency* of the linear layer was the hard part — proven for C4 (R7) +
  rot2 (constructive witness).  **No longer open; see R8-G.**
- **m=37 rot4 existence (= SIRH Part III, C4 instance)** — ⇔ satisfiability of
  the R8-G quadratic CSP (**1,264,378** per-line at-most-2 constraints over 1,369
  binary vars).  Encoding proven sound + reachable (45/45 `verify=True`;
  `--seed-known` C4/C2 admits cached solutions).  **Under resilient exact CP-SAT
  attack** via the proven R8-G model (background task `HNdAZl`, 7200s, 8 workers,
  checkpoint `results/m37_ckpt.json` → survives container teardown; log
  `results/m37_slice.log`).  **Faster next step (R9b):** re-host the same exact
  (X)+(S) characterization on a *permutation/assignment* CP-SAT model (2-regular
  structure enforced for free) instead of the 1.26M-constraint per-line model.
  **OPEN** (decidability question, not a board search).
- **Guy–Kelly D(n)=2n** infinite family — holds for *some* n≤72, not all.
- **C4 necessity** ⇔ no collinearity 2-factor (Th-39 direction).
- **α lower bound ~1.5n** (Hall 1975) for general NTIL.
- **Extremal-config taxonomy** — beyond symmetric families: are there
  asymmetric / sporadic structure types? (left open by the SE revision)
- **Container ⇒ symmetry bridge** — does a high-min-degree container force a
  non-trivial symmetry `G` (hence fall under FDR / Layer 1)?  Unproved.
- **SIRH Part IV (reverse, Sidon→symmetry?)** — **RESOLVED negative (2026-07-13).**
  Sidon is a strict one-way necessary condition: asymmetric configs satisfy it
  (Lemma E / P4.1), and Sidon carriers are *dominated* by asymmetric configs
  (counting, P4.2); even the full quadratic signature is one-way (P4.3).  SIRH
  is a one-way hierarchy, not an equivalence.  See `part4_reverse.md`.
- **Costas orders 32/33** — famous open frontier; classified by C1–C6/R/C
  (symmetric witness, if any, must be D/AD; rotational & axis symmetry
  impossible).  Exact CP-SAT attacks (`cpsat_costas.py --n 32/33`) were
  **deprioritized 2026-07-13 per user pivot** (stop horizontal expansion; the
  Costas cross-application already yielded its blockbuster C6/R/C and the *same*
  rigidity hierarchy transfers — see SIRH §7).  Left OPEN; not the current focus.

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
★ SIRH (Symmetry-Induced Rigidity Hierarchy — the umbrella theorem)
symmetry G ≤ D₄  (the group action)
   → orbit decomposition → natural fundamental domain F_G
   → Part I LINEAR rigidity (FDR, THEOREM): slope±1 preservation ⇒ Sidon  [necessary]
   → Part II: quadratic layer is FORCED (proven group-independently; strict
              insufficiency for C4 = R7) — FDR sees only slope±1, the rest is quadratic
   → Part III GEOMETRIC rigidity (R8-G, THEOREM for ALL FDR groups): per-line
              weighted at-most-2 ⇔ quadratic CSP (det≠0) ⇔ NTIL  [Hypothesis H RESOLVED]
   → Part IV reverse (Sidon→symmetry?): RESOLVED NEGATIVE (one-way only)
   → m=37 existence = Part III, C4 instance (1.26M per-line constraints)  [OPEN, exact attack running]
```
**Cross-application (NEW 2026-07-13):** the FDR symmetry lens transfers to the
Costas array problem — same D₄ lattice, complementary classification:
```
symmetry G ≤ D₄  (Costas)
   → C1: H/V reflection impossible (permutation-matrix structure)  [= FDR ort1 boundary]
   → C5: only 6 admissible types; full D₄ impossible
   → C6 + R: NO rotational symmetry for n≥3 (C4 empty all n>1; pure C2 empty; D2 only n=2)
   → C (capstone): for n≥3, G ∈ {id, D, AD} — only diagonal involutions occur
```
See `costas_symmetry_theorem.md`, `costas_rigidity.md`.
**Unification (NEW 2026-07-13 night):** Golomb ruler / Sidon set / Costas array
/ difference set are one family of "distinct-difference" lattice designs,
differing only in how densely differences may collide (sparse Sidon ≤1 vs dense
planar difference set =1 ⇒ projective plane).  Bolstered by new theorems
(**LQ** linear/quadratic split; **R8-C** C4-symmetric Costas = linear CSP;
**C6 + R** rotational Costas classification).  Key literature bridge: circular
Costas ⇔ Welch, *proved* via direct-product difference sets → finite projective
planes (Golomb–Moreno, 2015) — so Costas and 2D difference sets are **provably**
linked, not merely analogous.  Computational probes: Welch works to order 18;
C4-symmetric Costas **absent to n=32** (CP-SAT sweep, INFEASIBLE); C2-symmetric
**absent n=4..10** (enumeration) — both now THEOREMS (C6, R), not open.  See
`sidon_costas_unification.md`, `costas_symmetric_search.py`,
`costas_c2_corroborate.py`.
See `two_layer_rigidity.md` for the full synthesis.
