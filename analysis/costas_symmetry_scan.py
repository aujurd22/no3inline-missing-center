"""
costas_symmetry_scan.py
=======================
Apply the D4 symmetry-classification lens developed for the FDR theorem
(see analysis/results/fdr_theorem.md) to the COSTAS ARRAY problem.

A Costas array of order n is a permutation pi of {0..n-1} such that all
n(n-1)/2 displacement vectors (pi[j]-pi[i], j-i)  (i<j) are distinct --
geometrically, a permutation matrix whose pairwise differences are unique
(a 2-D Sidon set).

We enumerate ALL Costas arrays up to a modest order by backtracking, then
classify each by its symmetry group G <= D4 (acting on the n x n grid
centred at ((n-1)/2,(n-1)/2)).  This transfers our "symmetry-induced
rigidity" methodology from NTIL (FDR/R8) to Costas, and produces a real
table of which symmetry types admit Costas arrays.

Key theoretical point: a permutation matrix can NEVER be invariant under a
horizontal/vertical reflection (such a reflection keeps row i but flips the
column, forcing pi(i) = (n-1)/2 for every i -- impossible).  So only
{id, C2 (180deg), C4 (90deg), diagonal-involution, anti-diagonal} can occur.
This is the Costas analogue of the FDR boundary (ort1, which swaps the
slope+-1 families, is exactly the forbidden horizontal/vertical reflection).

Usage:
    python costas_symmetry_scan.py [max_n]
"""
import sys

def costas_enumerate(n, results):
    """Backtracking enumeration of all Costas arrays of order n."""
    cols_used = [False] * n
    perm = []
    diffs = set()          # (dx, dy) with dy > 0  (each unordered pair once)

    def rec(row):
        if row == n:
            results.append(tuple(perm))
            return
        for c in range(n):
            if cols_used[c]:
                continue
            new = []
            ok = True
            for j in range(row):
                d = (c - perm[j], row - j)
                if d in diffs:
                    ok = False
                    break
                new.append(d)
            if not ok:
                continue
            for d in new:
                diffs.add(d)
            cols_used[c] = True
            perm.append(c)
            rec(row + 1)
            perm.pop()
            cols_used[c] = False
            for d in new:
                diffs.discard(d)
    rec(0)

def symmetry_group(perm, n):
    """Return the set of D4 symmetries satisfied by permutation pi."""
    c2 = all(perm[n - 1 - i] == n - 1 - perm[i] for i in range(n))          # 180 rot
    c4 = all(perm[n - 1 - perm[i]] == i for i in range(n))                  # 90 rot
    dia = all(perm[perm[i]] == i for i in range(n))                         # diagonal refl (involution)
    adia = all(perm[n - 1 - perm[i]] == n - 1 - i for i in range(n))        # anti-diagonal refl
    elems = {'id'}
    if c2:   elems.add('C2')
    if c4:   elems.update(['C2', 'C4', 'C4b'])
    if dia:  elems.add('D')
    if adia: elems.add('AD')
    if dia and adia:                 # D * AD = 180
        elems.update(['C2', 'D', 'AD'])
    if c4 and dia:                   # full D4
        elems.update(['C2', 'C4', 'C4b', 'D', 'AD'])
    return frozenset(elems)

LABEL = {
    frozenset({'id'}): 'trivial',
    frozenset({'id', 'C2'}): 'C2 (180)',
    frozenset({'id', 'C2', 'C4', 'C4b'}): 'C4 (90)',
    frozenset({'id', 'D'}): 'diag-involution',
    frozenset({'id', 'AD'}): 'anti-diag',
    frozenset({'id', 'C2', 'D', 'AD'}): 'D2 (diag refl)',
    frozenset({'id', 'C2', 'C4', 'C4b', 'D', 'AD'}): 'D4 (full)',
}

def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 17
    print(f"{'n':>3} {'total':>8} {'trivial':>9} {'C2':>6} {'C4':>6} "
          f"{'diag':>6} {'anti':>6} {'D2':>5} {'D4':>5}   examples(C4/diag)")
    known_ok = True
    for n in range(1, max_n + 1):
        res = []
        costas_enumerate(n, res)
        if n <= 10:
            # sanity: known Costas counts C(1..10) = 1,2,4,12,40,116,200,444,760,2160
            expected = {1: 1, 2: 2, 3: 4, 4: 12, 5: 40, 6: 116,
                        7: 200, 8: 444, 9: 760, 10: 2160}
            if res and len(res) != expected[n]:
                known_ok = False
                print(f"   !! n={n} count {len(res)} != expected {expected[n]}")
        tally = {}
        examples = []
        for pi in res:
            g = symmetry_group(pi, n)
            tally[g] = tally.get(g, 0) + 1
            if g in (frozenset({'id', 'C2', 'C4', 'C4b'}),
                     frozenset({'id', 'D'}),
                     frozenset({'id', 'C2', 'C4', 'C4b', 'D', 'AD'})) and len(examples) < 2:
                examples.append((LABEL[g], pi))
        triv = tally.get(frozenset({'id'}), 0)
        c2 = tally.get(frozenset({'id', 'C2'}), 0)
        c4 = tally.get(frozenset({'id', 'C2', 'C4', 'C4b'}), 0)
        dia = tally.get(frozenset({'id', 'D'}), 0)
        adia = tally.get(frozenset({'id', 'AD'}), 0)
        d2 = tally.get(frozenset({'id', 'C2', 'D', 'AD'}), 0)
        d4 = tally.get(frozenset({'id', 'C2', 'C4', 'C4b', 'D', 'AD'}), 0)
        ex = "; ".join(f"{l}:{list(p)}" for l, p in examples)
        print(f"{n:>3} {len(res):>8} {triv:>9} {c2:>6} {c4:>6} "
              f"{dia:>6} {adia:>6} {d2:>5} {d4:>5}   {ex}")
    if known_ok:
        print("\nSanity check vs known Costas counts (n<=10): PASS")
    else:
        print("\nSanity check: MISMATCH (enumerator bug?)")

if __name__ == '__main__':
    main()
