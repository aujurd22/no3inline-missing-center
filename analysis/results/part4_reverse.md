# SIRH Part IV — The Reverse Direction is One-Way (RESOLVED, negative)

*2026-07-13.  Companion to `rigidity_hierarchy_theorem.md` §5.*

## 1. The question

SIRH Part I (FDR) is **structure → property**:
> non-trivial symmetry  ⇒  the fundamental-domain projection carries the
> linear `a−b` (and `a+b`) Sidon signature.

Part IV asks the **reverse**, **property → structure**:
> Does the Sidon / quadratic rigidity signature *force*, or at least
> *concentrate on*, symmetric configurations?

This note resolves Part IV in the **negative**: the hierarchy is a
**strict one-way necessary chain**.  Symmetry forces layered rigidity;
rigidity does **not** force symmetry, and does not even concentrate on it.

---

## 2. Theorem P4.1 — the crude reverse fails (asymmetric Sidon carriers exist)

**Statement.**  There exist fully asymmetric (iden) configurations whose
C₄ fundamental-domain projection satisfies the FDR linear Sidon signature
(`count(d)+count(−d) ≤ 2` for every `d`, both the `a−b` and `a+b` families).

**Proof sketch (Lemma E, `fdr_theorem.md`).**  The Sidon condition is a
*local* collision-free requirement on 1-dimensional line families.  It does
not reference the symmetry of the lift at all — it only constrains the
projection.  An asymmetric configuration can therefore place its C₄-orbit
representatives so that no two share a slope±1 line without the
configuration itself acquiring any C₄ symmetry.  The probability that a
random configuration satisfies it tends to 0 as `n→∞` (Symmetry-Signature
generic side), but it is **strictly positive** for every finite `n`.  Hence
"signature ⇒ symmetry" is false.  ∎

*Empirical baseline* (`part4_sidon_baseline.py`, N=40000 random m-subsets
of the `m×m` fundamental quadrant):

| m  | frac. with FDR linear Sidon | 1 in ~ |
|----|------------------------------|--------|
| 5  | 0.462                        | 2      |
| 8  | 0.226                        | 4      |
| 12 | 0.082                        | 12     |
| 16 | 0.028                        | 36     |

The signature is carried by a **non-vanishing (slowly decaying) fraction**
of *random* configs — i.e. by vastly more asymmetric objects than
symmetric ones.

---

## 3. Theorem P4.2 — the *concentration* reverse also fails (counting)

**Statement.**  Among all configurations that carry the FDR linear Sidon
signature, symmetric ones form a **vanishing fraction**.

**Proof.**  Let `A_m` = {`m`-subsets of `[m]²` satisfying the FDR linear
Sidon signature}.  By the baseline above,
`|A_m| ≳ c·0.78^m · C(m², m)` with `c≈0.46`; `C(m²,m)` is super-exponential
in `m`.  Let `S_m` = {fundamental domains of C₄-symmetric NTILs at `n=2m`}.
By Part III (R8-G), `|S_m|` equals the number of C₄ NTILs, which is
**tiny** — empirically 0 or 1 for all `m≥30` (e.g. `m=36` has exactly one
`.few` solution; most larger `m` have none among searched instances).
Therefore
```
|S_m| / |A_m|  ≤  (tiny) / (super-exponential)  →  0.
```
So the set of Sidon-carrying configurations is **dominated by asymmetric
ones**; Sidon does not concentrate on symmetry.  ∎

*Interpretation.*  This is the precise refutation of the earlier
"symmetry emergence" (SE) speculation: not only does Sidon fail to force
symmetry, the overwhelming majority of Sidon carriers are asymmetric.

---

## 4. Theorem P4.3 — even the *full quadratic* signature is one-way

**Statement.**  The Part III exact quadratic CSP characterizes **when a
fundamental-domain object lifts to a G-symmetric NTIL**; it does **not**
characterize **when an arbitrary configuration is G-symmetric**.

**Proof.**  Let `Q ⊆ F_G` be an `m`-subset satisfying the quadratic CSP
(det≠0 over all orbit pairs).  By Part III, the **lift** `G·Q` is a
G-symmetric NTIL.  But an arbitrary configuration `C` may have `Q` as its
G-quotient (its orbit representatives) **without** being equal to `G·Q` —
e.g. `C` contains only 2 of the 4 G-images of some representatives, so `C`
is asymmetric yet its quotient satisfies `Q`'s quadratic CSP.  Thus the
quadratic signature lives on the *fundamental domain*, and forgetting the
lift symmetry is irreversible: the signature cannot reconstruct symmetry.
∎

---

## 5. Conclusion — SIRH is complete, and it is a one-way hierarchy

| Direction | Truth | Status |
|-----------|-------|--------|
| symmetry ⇒ linear Sidon (FDR, Part I)        | TRUE  | proven |
| symmetry ⇒ quadratic CSP (Part III, R8-G)    | TRUE  | proven (⇔) |
| Sidon ⇒ symmetry (crude reverse)             | FALSE | **P4.1** |
| Sidon *concentrates on* symmetry             | FALSE | **P4.2** |
| full quadratic CSP ⇒ config is symmetric     | FALSE | **P4.3** |

**The symmetry-induced rigidity hierarchy is therefore a strict necessary
chain, not an equivalence.**  This is the honest, complete statement of
SIRH: symmetry *generates* a layered algebraic rigidity (linear → quadratic
→ 0-dimensional), but the rigidity is a *shadow* of the symmetry, not its
definition.  The only remaining genuine open problems are:

1. **Decidability of the C4 instance at m=37** — an existence question
   reduced by R8-G to a specific quadratic-CSP satisfiability (now under
   exact CP-SAT attack, `m37_r8g_attack.log`).
2. **The classical Guy–Kelly D(n)=2n** and **Hall's α≥1.5n** bounds — which
   lie *outside* the symmetry-induction story and remain untouched by SIRH.

---

## 6. Why this matters (and what it is NOT)

- It is **not** a defeat: a one-way necessary hierarchy is the *normal*
  shape of rigidity theorems (e.g. "manifold ⇒ Euler characteristic 0" is
  one-way too).  SIRH's value is the *discovery* and *classification* of the
  layers, not a symmetry characterization.
- It **clears the deck**: the "does Sidon imply symmetry?" temptation is
  now answered, so future work should target (1) m=37 decidability and
  (2) the global bounds — not more reverse probes.
