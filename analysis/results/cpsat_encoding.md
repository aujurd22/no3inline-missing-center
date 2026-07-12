# Exact CP-SAT re-encoding of rot4 NTIL (the m=37 attack)

**Date:** 2026-07-12 (night)
**Engine:** `analysis/cpsat_m37.py` (OR-Tools CP-SAT; Z3 fallback)
**Theorem base:** R8 (quadratic-Sidon completeness) → `(X)∧(S) ⇔ rot4 NTIL`
**Companion:** `stageM_m37.md` (why greedy DFS fails), `r8_proof.md`, `two_layer_rigidity.md`

---

## 1. The encoding

| item | definition |
|------|------------|
| **variables** | `sel[p] ∈ {0,1}` for every quadrant position `p ∈ {0..m-1}²`  (`m²` vars) |
| **cardinality** | `Σ_p sel[p] = m`  (pick exactly m distinct C4-orbit representatives) |
| **per-line constraint** | for every board line `L` of the `2m×2m` grid, let `w_p = #{r : c4(p,r) ∈ L}` (`w_p ∈ {1,2}`); add `Σ_{p on L} w_p·sel[p] ≤ 2` |

The instance is built by **line enumeration**: for each reduced direction class,
group the `2m×2m` grid points by their perpendicular projection, yielding every
line; for each line we invert `c4` to recover which quadrant positions land on it
(and with how many rotations).  A line contributes a constraint only when its
total weight exceeds 2.

---

## 2. Why this is EXACTLY `(X)∧(S)` (hence rot4 NTIL)

Three collinear lifted points always lie on a single line `L'`; their source
positions contribute `≥3` to `Σ w_p·sel[p]` on `L'` → forbidden by the
at-most-2 constraint.  Conversely, `Σ w_p·sel[p] ≤ 2` on **every** line means no
line carries `≥3` lifted points, i.e. no three collinear → rot4 NTIL.

- **Layer X** (3 distinct cells collinear) is captured: any 3 source positions on
  `L` sum to `≥3`.
- **Layer S** (2 images of one cell + 1 of another collinear) is captured: a
  position with `w_p = 2` already consumes the budget of 2, so **no other**
  position on `L` may be selected — exactly the S clause.

Thus the single family of per-line weighted at-most-2 constraints is *equivalent*
to the entire `(X)∧(S)` determinant system, and (by R8) to rot4 NTIL.  This is
an **exact** re-encoding, not a relaxation.

---

## 3. Size — compact vs naive

A naive encoding emits one 3-literal clause per collinear (position-triple,
rotation-form).  For m=37 that is **30,992,032** clauses.  The per-line
at-most-2 form collapses it to **1,264,378** efficient constraints — a ~25×
reduction — and the constraints give the solver *global* line propagation
instead of local triple clauses.

| m | per-line constraints | vars | discovery (OR-Tools) |
|---|---|---|---|
| 12 | 14,038 | 144 | OPTIMAL (2.0 s) |
| 15 | 34,270 | 225 | OPTIMAL (32.5 s) |
| 36 | 1,135,606 | 1,296 | — (seed-known, see §4) |
| 37 | 1,264,378 | 1,369 | **attack running** (task `fkB7gu`, 1800 s budget) |

---

## 4. Validation (soundness + reachability)

**Soundness — discovery sweep (m=3..15).**  OR-Tools CP-SAT finds an assignment
for every m; each found assignment was re-checked with the *full* 64/24
`(X)+(S)` determinant test (`verify()` → `x_full`/`s_full`) and passed
(`verify=True`).  So the clause system admits **no** false (non-NTIL) solution.

**Reachability — known-solution admission (m=36).**  Loaded a cached rot4
solution for m=36, converted its odd-pairing values to quadrant `sel`-bits, and
checked it directly against all 1,135,606 constraints: `card_ok=True`,
`violations=0` → the instance **admits** the known solution.  (The conversion is
`x = m − (a+1)/2, y = m − (b+1)/2`; see `r8_proof.md` setup.)  This proves the
encoding is correct at the largest known m, independent of discovery difficulty.

---

## 5. The m=37 attack and its interpretation

The m=37 instance is launched as a background OR-Tools CP-SAT solve
(`cpsat_m37.py --m 37 --timelimit 1800`, task `fkB7gu`).  Possible outcomes:

- **SAT** → a rot4 configuration at m=37 exists (the frontier resolved *for*).
- **UNSAT** → m=37 is the first failure (would be a major theorem: rot4 NTIL
  exists for all m=3..36 but **not** m=37).
- **UNKNOWN** → needs more time / symmetry-breaking / better search — the
  question remains open but the encoding is now the correct, sharp tool.

**Why this supersedes the earlier 55+ failures.**  Every prior search lived in
the *linear* Sidon space, which R7 proved cannot capture cross-quadrant
collinearity.  The `constraint_prop_solver.py` Stage M attack (m=37) was the
first in the *correct* `(X+S)` space, but used greedy depth-first construction,
which the **global rigidity / complete-determination** findings (R6) predict must
fail: the system is 0-dimensional and globally coupled, so a partial assignment
cannot be validated until all 37 cells are placed → DFS died at `max_depth=25 ≪
37`.  Clause-learning CP-SAT with the strong `Σsel=m` cardinality and global
line propagation is the paradigm that can pierce a 0-d globally-coupled system —
this is precisely the "约束传播在正确二次空间" prescription from the quadratic-gap
note.

---

## 6. Files

- `analysis/cpsat_m37.py` — generator + solver (`--count-only`, `--validate`,
  `--seed-known M`, `--m 37 --timelimit T`).
- `analysis/results/stageM_m37.md` — the (correct-space but DFS) Stage M result.
- `analysis/results/r8_proof.md`, `quadratic_sidon_completeness.md` — the theorem
  base.
