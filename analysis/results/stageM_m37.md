# Stage M — m=37 (n=74) existence attack in the CORRECT (X+S) quadratic space

**Date:** 2026-07-12 (night)
**Engine:** `analysis/constraint_prop_solver.py :: stage_M` (first solver living in the
R8-equivalent `(X) ∧ (S)` quadratic space — *not* the earlier linear-Sidon space).
**Parameters:** `m=37`, `n=74`, Monte-Carlo `mc_trials=2000`, backtracking cap
`bt_timeout=900.0 s`, `node_cap=20_000_000`.

---

## 1. Raw results

| phase | outcome |
|---|---|
| (a) Monte-Carlo | 2000 random m-subsets → **0** passed `(X+S)`  `[0.7 s]` |
| (b) Backtracking | `nodes = 629,861`, `max_depth = 25`, `found = 0`  `[915.1 s, cap 900 s]` |

Interpretation line printed by the engine:

> `max_depth << m` means `(X+S)` prunes the search far from any complete
> assignment — i.e. the Complete-Determination effect is already very strong at m=37.

---

## 2. Why this failure is *informative*, not discouraging

This is the **first** m=37 attack that runs in the space guaranteed (Theorem R8)
to be *exactly equivalent* to rot4 NTIL. Every one of the ~55 earlier searches
(Z3/CP-SAT/SA/Sidon/permutation/…) lived in the **linear** Sidon space, which
R7 proved *cannot in principle* capture cross-quadrant collinearity. Their
failures were therefore "wrong-space" failures and said nothing about existence.

Stage M's failure is of a *different, deeper* kind:

- **It is not a proof of non-existence.** `found=0` only shows that a greedy
  depth-first *construction* could not complete an m=37 assignment. Greedy
  construction was already predicted to fail by the **global-rigidity /
  complete-determination** findings (R6, and the 2026-07-12 20:30 notes):
  - the 16 quadratic forms constitute a **0-dimensional** system — once enough
    cells are fixed, `(X+S)` *uniquely* forces every remaining cell;
  - the *local* determination power **decays** with m (m=10 dist-3 ≈ 97% →
    m=36 dist-3 ≈ 32%), so the system is **globally coupled**, not decomposable
    into sub-problems;
  - consequently a greedy choice at depth *k* cannot be validated until all 37
    cells are placed → DFS dies at `max_depth=25 ≪ 37` and backtracks helplessly.

So the result is perfectly consistent with the standing picture: **m=37 is a
hard, globally-coupled, 0-dimensional CSP** — exactly the regime where
backtracking is the *wrong paradigm* but where **clause-learning SAT** is the
right one.

---

## 3. The correct next attack: a CP-SAT instance in `(X+S)`-clause space

Because collinearity of a *fixed* geometric triple is a **static boolean fact**,
the quadratic determinant conditions linearise perfectly into CNF over the
position-selection bits:

- variables: `sel[p] ∈ {0,1}` for every quadrant position `p ∈ {0..m-1}²`
  (`m² = 1369` vars for m=37);
- cardinality: `Σ_p sel[p] = m` (pick exactly m distinct orbit reps);
- **X-clause**: for every line `L` of the `2m×2m` board and every three distinct
  quadrant positions `p1,p2,p3` whose C4-images land on `L`, forbid co-selection:
  `¬sel[p1] ∨ ¬sel[p2] ∨ ¬sel[p3]`;
- **S-clause**: for every line `L` and every position `p` appearing with
  **≥2 distinct rotations** on `L` together with another position `q` appearing
  ≥1× on `L`, forbid: `¬sel[p] ∨ ¬sel[q]`.

Line-enumeration generates these clauses *completely* (every collinear lifted
triple lies on some line `L'`, and its three source positions are all in `P_{L'}`),
so the instance is an exact re-encoding of `(X)∧(S)`.  This is precisely the
"约束传播在正确二次空间" prescription from the quadratic-gap note.  Clause-learning
(CDCL) + the strong `Σsel=m` propagation is the tool that can pierce a 0-d
globally-coupled system where DFS cannot.

→ implemented in `analysis/cpsat_m37.py` (validated on small m where solutions
are known, then escalated to m=37 with a time budget).

---

## 4. Status

**m=37 rot4 existence remains OPEN** — neither proved nor disproved.
- Evidence *for* existence: all m=3..36 have rot4 solutions; structural invariants
  (source-ratio ≈0.26, longest-cycle ≈0.79·m) and the constraint-density window
  are continuous across m=36→37, so m=37 is *not* at a satisfiability phase
  transition.
- Evidence *against* (weak): naive random/Mc sampling hit rate ≈0 at m=37; DFS
  cannot construct one. Both are algorithmic limitations, not impossibility proofs.

The question is now a clean, finite, exactly-encoded CSP; the answer is decided
by the algebra of 140,304 (dedup) / 529,248 (full) quadratic non-equalities over
1,369 binary variables.
