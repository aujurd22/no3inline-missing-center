"""
characterize_x.py — Systematic empirical characterization of "X":
the gap between "pairing satisfies Sidon" and "pairing lifts to valid rot4 NTIL".

Three experiments:
1. GENERATE random Sidon-4 pairings, test how many also pass Sidon-8, then C4-lift
2. CLASSIFY all cross-quadrant collinearities in C4-lift failures
3. PRECISION-ID: which extra condition(s) separate Sidon from NTIL?
"""

import os, sys, glob, random, math
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rot4_loader as rot4
from itertools import combinations

CACHE = rot4.CACHE

# ────────────────────────────────────────────
# 0. C4 lift + validation (from pairing_construction)
# ────────────────────────────────────────────

def c4_lift(cells, n):
    pts = set()
    for (x, y) in cells:
        for p in [(x, y), (n-1-y, x), (n-1-x, n-1-y), (y, n-1-x)]:
            pts.add(p)
    return pts

def has_collinear(pts):
    pl = list(pts)
    for i in range(len(pl)):
        x1, y1 = pl[i]
        for j in range(i+1, len(pl)):
            x2, y2 = pl[j]
            dx, dy = x2-x1, y2-y1
            for k in range(j+1, len(pl)):
                x3, y3 = pl[k]
                if dx * (y3 - y1) == dy * (x3 - x1):
                    return True, (pl[i], pl[j], pl[k])
    return False, None

def cells_from_pairing(pairs, m):
    """Convert odd-number pairs (a,b) to quadrant cells (x,y)."""
    return [(m - (a + 1) // 2, m - (b + 1) // 2) for a, b in pairs]

def pairing_from_cells(cells, m):
    """Inverse: cells → pairing."""
    return [(2*(m-x)-1, 2*(m-y)-1) for (x,y) in cells]

def check_2cover(pairs, m):
    all_vals = []
    for a, b in pairs:
        all_vals.append(a); all_vals.append(b)
    expected = sorted([2*i+1 for i in range(m)] * 2)
    return sorted(all_vals) == expected

def validate_rot4_c4lift(pairs, m):
    """Check if a pairing lifts to a valid rot4 NTIL solution."""
    n = 2 * m
    if not check_2cover(pairs, m):
        return False, "2-cover", None
    for a, b in pairs:
        if a == b:
            return False, "self-loop", None
    cells = cells_from_pairing(pairs, m)
    pts = c4_lift(cells, n)
    if len(pts) != 4 * m:
        return False, f"pts={len(pts)}", None
    bad, triple = has_collinear(pts)
    if bad:
        return False, "collinear", triple
    rc = Counter(p[1] for p in pts)
    cc = Counter(p[0] for p in pts)
    if any(c != 2 for c in rc.values()):
        return False, "row", None
    if any(c != 2 for c in cc.values()):
        return False, "col", None
    return True, "VALID", None

# ────────────────────────────────────────────
# 0b. Sidon check on pairings (a,b parameterization)
# ────────────────────────────────────────────

FORMS_4 = [(1, -1), (1, 1)]        # a-b, a+b  (FDR core)
FORMS_8 = [(2, -1), (-1, 2),       # 2a-b, a-2b
           (2, 1), (1, 2),          # 2a+b, a+2b
           (3, -1), (-1, 3)]        # 3a-b, a-3b

def check_forms(pairs, forms, m):
    """Check Sidon-like condition: each |c1*a + c2*b| value appears ≤ 2 times."""
    n = 2 * m
    for c1, c2 in forms:
        vals = [abs(c1*a + c2*b) for a, b in pairs]
        cnt = Counter(vals)
        # Normalize: count(d) + count(-d) ≤ 2
        abs_cnt = Counter()
        for v, c in cnt.items():
            abs_cnt[v] += c
        if any(c > 2 for c in abs_cnt.values()):
            return False
    return True

def check_4form(pairs, m):
    return check_forms(pairs, FORMS_4, m)

def check_8form(pairs, m):
    return check_forms(pairs, FORMS_8, m)

# ────────────────────────────────────────────
# 1. Generate random Sidon-4 pairings
# ────────────────────────────────────────────

def generate_random_2cover(m, seed=None):
    """Generate a random 2-cover pairing of odd numbers {1,3,...,2m-1}.
    Returns list of m unordered pairs (a,b) with a<b.
    Each odd number appears exactly twice across all pairs."""
    if seed is not None:
        random.seed(seed)
    odds = [2*i+1 for i in range(m)] * 2  # each twice
    random.shuffle(odds)
    pairs = []
    for i in range(m):
        a = odds[2*i]
        b = odds[2*i+1]
        if a == b:
            continue  # self-loop, reject
        pairs.append((a, b) if a < b else (b, a))
    # Check 2-cover
    if check_2cover(pairs, m):
        return pairs
    return None

def generate_random_sidon4(m, max_attempts=10000, seed_offset=0):
    """Generate a random pairing that satisfies 2-cover AND 4-form Sidon."""
    for attempt in range(max_attempts):
        pairs = generate_random_2cover(m, seed=attempt + seed_offset)
        if pairs is None:
            continue
        if check_4form(pairs, m):
            return pairs
    return None

def generate_random_sidon8(m, max_attempts=20000, seed_offset=0):
    """Generate random pairing satisfying 2-cover + 4-form + 8-form Sidon."""
    for attempt in range(max_attempts):
        pairs = generate_random_2cover(m, seed=attempt + seed_offset)
        if pairs is None:
            continue
        if check_4form(pairs, m) and check_8form(pairs, m):
            return pairs
    return None

# ────────────────────────────────────────────
# 2. Classify collinearity patterns in failures
# ────────────────────────────────────────────

def classify_collinear_triple(triple, pair_list, m, n):
    """Given a collinear triple of points from different C4 orbits,
    determine which pairings and which quadrants are involved."""
    def classify_point(p, cells_set):
        x, y = p
        if (x, y) in cells_set:
            return "Q0"
        if (n-1-y, x) in cells_set:
            return "Q1"
        if (n-1-x, n-1-y) in cells_set:
            return "Q2"
        if (y, n-1-x) in cells_set:
            return "Q3"
        return "?"
    
    cells = [(m - (a+1)//2, m - (b+1)//2) for a,b in pair_list]
    cells_set = set(cells)
    
    p1, p2, p3 = triple
    q1, q2, q3 = classify_point(p1, cells_set), classify_point(p2, cells_set), classify_point(p3, cells_set)
    
    # Determine which pairings are involved
    involved_pairs = set()
    for p in triple:
        x, y = p
        for idx, (cx, cy) in enumerate(cells):
            if (x, y) in [(cx, cy), (n-1-cy, cx), (n-1-cx, n-1-cy), (cy, n-1-cx)]:
                involved_pairs.add(idx)
                break
    
    qs = tuple(sorted([q1, q2, q3]))
    
    return {
        'quadrants': qs,
        'num_involved_pairs': len(involved_pairs),
        'involved_pairs': sorted(involved_pairs),
    }

# ────────────────────────────────────────────
# 3. Main: Systematic gap measurement
# ────────────────────────────────────────────

def experiment1(m, num_sidon4=200, num_sidon8=200):
    """Experiment 1: Convergence rates Sidon-4 → Sidon-8 → C4-lift valid.
    For a given m, try to generate random pairings at each level."""
    results = {
        'm': m,
        'n': 2*m,
        's4_attempts': 0,
        's4_found': 0,
        's8_attempts': 0,
        's8_found': 0,
        'lift_attempts': 0,
        'lift_valid': 0,
        'lift_fail_reasons': Counter(),
        'failure_types': Counter(),
    }
    
    # Try to find Sidon-4 pairings
    base_seed = m * 10000
    s4_list = []
    for attempt in range(num_sidon4 * 5):  # up to 5x attempts
        pairs = generate_random_sidon4(m, max_attempts=50, seed_offset=base_seed + attempt * 200)
        if pairs is not None:
            s4_list.append(pairs)
            results['s4_found'] += 1
            if len(s4_list) >= num_sidon4:
                break
    results['s4_attempts'] = base_seed + (attempt + 1) * 200
    
    # Check how many also pass Sidon-8
    s8_list = []
    for pairs in s4_list:
        results['s8_attempts'] += 1
        if check_8form(pairs, m):
            s8_list.append(pairs)
            results['s8_found'] += 1
    
    # C4-lift test: test ALL Sidon-8, plus up to some Sidon-4-only for comparison
    for pairs in s8_list[:num_sidon8]:
        results['lift_attempts'] += 1
        status, reason, triple = validate_rot4_c4lift(pairs, m)
        if status:
            results['lift_valid'] += 1
        else:
            results['lift_fail_reasons'][reason] += 1
            if reason == 'collinear' and triple is not None:
                info = classify_collinear_triple(triple, pairs, m, 2*m)
                ft = f"Q={info['quadrants']} P={info['num_involved_pairs']}"
                results['failure_types'][ft] += 1
    
    # Also test Sidon-4-only (non-Sidon-8) for comparison
    s4_only = [p for p in s4_list if not check_8form(p, m)]
    for pairs in s4_only[:50]:
        status, reason, triple = validate_rot4_c4lift(pairs, m)
        if not status:
            results['lift_fail_reasons'][f"s4only_{reason}"] += 1
    
    return results

def experiment2(valid_pairs_list, m):
    """Experiment 2: From known valid solutions, systematically modify
    1-3 pairings and test where it breaks. This pinpoints the EXACT constraints."""
    n = 2 * m
    modified_results = {
        'm': m,
        'tested': 0,
        'still_valid': 0,
        'break_reasons': Counter(),
    }
    
    for base_pairs in valid_pairs_list[:10]:  # limit to 10 base solutions
        # Strategy: pick a pair, swap its b with another pair's b
        for i in range(m):
            for j in range(i+1, m):
                modified_results['tested'] += 1
                mp = list(base_pairs)
                # Swap b values between pair i and pair j
                ai, bi = mp[i]
                aj, bj = mp[j]
                mp[i] = (ai, bj) if ai < bj else (bj, ai)
                mp[j] = (aj, bi) if aj < bi else (bi, aj)
                
                # Check basic validity
                if not check_2cover(mp, m):
                    modified_results['break_reasons']['2cover_after_swap'] += 1
                    continue
                if not check_4form(mp, m):
                    modified_results['break_reasons']['4form_after_swap'] += 1
                    continue
                if not check_8form(mp, m):
                    modified_results['break_reasons']['8form_after_swap'] += 1
                    continue
                    
                status, reason, triple = validate_rot4_c4lift(mp, m)
                if status:
                    modified_results['still_valid'] += 1
                else:
                    modified_results['break_reasons'][f'c4lift_{reason}'] += 1
    
    return modified_results


def main():
    print("=" * 90)
    print("EXPERIMENT 1: Sidon convergence + C4-lift survival rates")
    print("=" * 90)
    
    for m in range(4, 13):  # m=4..12 manageable for random search
        print(f"\n--- m={m} (n={2*m}) ---")
        r = experiment1(m, num_sidon4=100, num_sidon8=100)
        
        s8_rate = 100.0 * r['s8_found'] / r['s4_found'] if r['s4_found'] else 0
        lift_rate = 100.0 * r['lift_valid'] / r['lift_attempts'] if r['lift_attempts'] else 0
        
        print(f"  Sidon-4 found: {r['s4_found']} / {r['s4_attempts']} attempts")
        print(f"  Sidon-8 found: {r['s8_found']} / {r['s4_found']} S4 pairs = {s8_rate:.1f}%")
        print(f"  C4-lift valid: {r['lift_valid']} / {r['lift_attempts']} S8 pairs = {lift_rate:.1f}%")
        print(f"  Fail reasons: {dict(r['lift_fail_reasons'])}")
        
        if r['failure_types']:
            print(f"  Collinear patterns (top 5):")
            for pat, cnt in r['failure_types'].most_common(5):
                print(f"    {pat}: {cnt}")
    
    # Experiment 2: Modify from known valid solutions
    print("\n" + "=" * 90)
    print("EXPERIMENT 2: Systematic modification of known valid solutions")
    print("=" * 90)
    
    for m_target in [5, 6, 7, 8]:
        files = sorted(glob.glob(os.path.join(CACHE, f'n{2*m_target}_rot4*')))
        if not files:
            continue
        path = files[0]
        with open(path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        valid_pairs_list = []
        for line in lines[:20]:
            pts = rot4.decode_line(line, 2*m_target)
            if len(pts) != 4*m_target or not rot4.is_valid(pts, 2*m_target):
                continue
            pairs = [(2*(m_target-x)-1, 2*(m_target-y)-1) for (x,y) in pts if x < m_target and y < m_target]
            if len(pairs) == m_target:
                valid_pairs_list.append(pairs)
        
        if valid_pairs_list:
            r2 = experiment2(valid_pairs_list, m_target)
            print(f"\nm={m_target}: {r2['tested']} swaps tested, {r2['still_valid']} still valid")
            print(f"  Break reasons: {dict(r2['break_reasons'])}")
    
    print("\nDone.")


if __name__ == '__main__':
    main()
