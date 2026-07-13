"""
surgery_m36_to_37_v3.py — Direction 2: BUG-FIXED surgery with correct N=74 lift

v1 bug (v2): incremental check used wrong N=72 base points, producing false positives.
v3 fix: ALL 37 cells lifted with N=74, brute-force collinear check on 148 points.

Strategy:
  Phase 1: 1-edge split (36 attempts)
  Phase 2: 2-switch + split any of the 36 edges
  Phase 3: 3-switch + split any of the 36 edges
  Phase 4: random k-switch + split (50000 attempts)
"""

import sys, os, json, time, random, math
from itertools import permutations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rot4_loader import decode_line

rng = random.Random(20260713)

# ── Geometry ──────────────────────────────────────────────────────────
def c4(p, r, N):
    x, y = p
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    if r == 3: return (y, N - 1 - x)

def check_no_collinear(cells, m):
    """Brute-force check: C4-lift 4m points with N=2m, verify no 3 collinear."""
    N = 2 * m
    pts = []
    for (x, y) in cells:
        for r in range(4):
            pts.append(c4((x, y), r, N))
    n = len(pts)
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            dx1, dy1 = xj - xi, yj - yi
            for k in range(j + 1, n):
                xk, yk = pts[k]
                dx2, dy2 = xk - xi, yk - yi
                if dx1 * dy2 == dy1 * dx2:
                    return False  # collinear
    return True

# ── 2-factor helpers ──────────────────────────────────────────────────
def cells_from_sel(sel, m):
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]

def verify_2factor(sel, m):
    for i in range(m):
        if sum(sel[i]) + sum(sel[r][i] for r in range(m)) != 2:
            return False
    return True

def get_edges(sel, m):
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]


def apply_k_switch(sel_old, m_old, edge_set, bijection):
    """Remove edge_set, add bijection edges. Returns None if conflict."""
    sel = [row[:] for row in sel_old]
    for (i, j) in edge_set:
        sel[i][j] = 0
    for (src, tgt) in bijection:
        if sel[src][tgt] == 1:
            return None
        sel[src][tgt] = 1
    return sel


def try_split(sel_mid, m_old, split_edge):
    """
    Given m_old×m_old selection matrix, try splitting split_edge=(a,b)
    by inserting vertex m_old. Returns (ok, cells) or (False, None).
    ALL 37 cells lifted with N=74, full brute-force check.
    """
    m_new = m_old + 1
    a, b = split_edge

    # Build new selection matrix
    sel = [row[:] for row in sel_mid]
    sel.append([0] * m_new)
    for r in range(m_new):
        sel[r].append(0)

    if sel[a][b] != 1:
        return False, None
    sel[a][b] = 0
    sel[a][m_old] = 1
    sel[m_old][b] = 1

    if not verify_2factor(sel, m_new):
        return False, None

    cells = cells_from_sel(sel, m_new)
    if len(cells) != m_new:
        return False, None

    # Full brute-force check with N=2*m_new
    if not check_no_collinear(cells, m_new):
        return False, None

    return True, cells


# ══════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════

def generate_vertex_disjoint_triples(edges):
    """Generate all 3-edge subsets of edges with 6 distinct vertices."""
    n = len(edges)
    for i in range(n):
        ei = edges[i]
        vi = {ei[0], ei[1]}
        for j in range(i + 1, n):
            ej = edges[j]
            if ej[0] in vi or ej[1] in vi:
                continue
            vj = vi | {ej[0], ej[1]}
            for k in range(j + 1, n):
                ek = edges[k]
                if ek[0] in vj or ek[1] in vj:
                    continue
                yield [ei, ej, ek]


def generate_vertex_disjoint_pairs(edges):
    """Generate all 2-edge subsets with 4 distinct vertices."""
    n = len(edges)
    for i in range(n):
        ei = edges[i]
        vi = {ei[0], ei[1]}
        for j in range(i + 1, n):
            ej = edges[j]
            if ej[0] in vi or ej[1] in vi:
                continue
            yield [ei, ej]


def main():
    m_old = 36
    m_new = 37

    # 1. Load m=36 solution
    line = open('flammenkamp_cache/n72_rot4.few').read().strip()
    pts = decode_line(line, 72)
    cells_old = [(x, y) for (x, y) in pts if x < m_old and y < m_old]
    cells_old.sort()
    sel_old = [[0]*m_old for _ in range(m_old)]
    for (x,y) in cells_old:
        sel_old[x][y] = 1
    assert verify_2factor(sel_old, m_old)

    edges_old = get_edges(sel_old, m_old)
    print(f"[surgery] Loaded m=36: {len(edges_old)} edges, verify={check_no_collinear(cells_old, m_old)}")

    total_attempts = 0
    found = []
    t0 = time.time()

    # ── Phase 1: simple split (36 attempts) ──────────────────────
    print(f"\n{'='*55}\nPhase 1: 1-edge split (36 attempts)")
    for (a,b) in edges_old:
        total_attempts += 1
        ok, cells = try_split(sel_old, m_old, (a,b))
        if ok:
            found.append(('split', (a,b), None, cells))
            print(f"  FOUND: split {a}→{b}")
            break
    if found:
        print(f"  Result: 1 success in {total_attempts} attempts")
    else:
        print(f"  Result: 0/36")

    # ── Phase 2: 2-switch + split ────────────────────────────────
    if not found:
        print(f"\n{'='*55}\nPhase 2: 2-switch + split")
        for edge_pair in generate_vertex_disjoint_pairs(edges_old):
            e1, e2 = edge_pair
            s1, t1 = e1
            s2, t2 = e2
            sources = [s1, s2]
            targets = [t1, t2]
            
            for perm in permutations(targets):
                bijection = list(zip(sources, perm))
                # Skip identity
                if bijection[0][1] == t1 and bijection[1][1] == t2:
                    continue
                sel_mid = apply_k_switch(sel_old, m_old, edge_pair, bijection)
                if sel_mid is None:
                    continue
                if not verify_2factor(sel_mid, m_old):
                    continue
                
                # Try splitting any of the 36 edges in the new 2-factor
                edges_mid = get_edges(sel_mid, m_old)
                for (a,b) in edges_mid:
                    total_attempts += 1
                    ok, cells = try_split(sel_mid, m_old, (a,b))
                    if ok:
                        found.append(('2switch+split', edge_pair, (a,b), cells))
                        break
                if found:
                    break
            if found:
                break
        if found:
            print(f"  FOUND: {found[-1][0]} in {total_attempts} attempts")
        else:
            print(f"  Result: 0 successes")

    # ── Phase 3: 3-switch + split ────────────────────────────────
    if not found:
        print(f"\n{'='*55}\nPhase 3: 3-switch + split")
        n_triples = 0
        for edge_trip in generate_vertex_disjoint_triples(edges_old):
            sources = [e[0] for e in edge_trip]
            targets = [e[1] for e in edge_trip]
            
            for perm in permutations(targets):
                bijection = list(zip(sources, perm))
                if all(b[1]==edge_trip[i][1] for i,b in enumerate(bijection)):
                    continue  # skip identity
                
                sel_mid = apply_k_switch(sel_old, m_old, edge_trip, bijection)
                if sel_mid is None:
                    continue
                if not verify_2factor(sel_mid, m_old):
                    continue
                
                # Check the switched 2-factor has no collinearity
                cells_mid = cells_from_sel(sel_mid, m_old)
                if not check_no_collinear(cells_mid, m_old):
                    continue  # switch itself created collinearity
                
                edges_mid = get_edges(sel_mid, m_old)
                for (a,b) in edges_mid:
                    total_attempts += 1
                    ok, cells = try_split(sel_mid, m_old, (a,b))
                    if ok:
                        found.append(('3switch+split', edge_trip, (a,b), cells))
                        break
                if found:
                    break
                n_triples += 1
            
            if found:
                break
            if n_triples % 500 == 0 and n_triples > 0:
                print(f"  triples={n_triples}, attempts={total_attempts}")
        
        if found:
            print(f"  FOUND: {found[-1][0]} in {total_attempts} attempts")
        else:
            print(f"  Result: 0 successes after {n_triples} triples, {total_attempts} attempts")

    # ── Phase 4: random k-switch + split ─────────────────────────
    if not found:
        print(f"\n{'='*55}\nPhase 4: random k-switch + split (50000 attempts)")
        for _ in range(50000):
            k = rng.randint(2, 4)
            # Pick k vertex-disjoint edges
            shuffled = list(edges_old)
            rng.shuffle(shuffled)
            edge_set = []
            used_v = set()
            for e in shuffled:
                if e[0] not in used_v and e[1] not in used_v:
                    edge_set.append(e)
                    used_v.add(e[0])
                    used_v.add(e[1])
                    if len(edge_set) == k:
                        break
            if len(edge_set) < k:
                continue

            sources = [e[0] for e in edge_set]
            targets = [e[1] for e in edge_set]
            perm = list(targets)
            rng.shuffle(perm)
            bijection = list(zip(sources, perm))

            sel_mid = apply_k_switch(sel_old, m_old, edge_set, bijection)
            if sel_mid is None:
                continue
            if not verify_2factor(sel_mid, m_old):
                continue

            cells_mid = cells_from_sel(sel_mid, m_old)
            if not check_no_collinear(cells_mid, m_old):
                continue

            edges_mid = get_edges(sel_mid, m_old)
            rng.shuffle(edges_mid)
            for (a,b) in edges_mid[:20]:  # try up to 20 split targets
                total_attempts += 1
                ok, cells = try_split(sel_mid, m_old, (a,b))
                if ok:
                    found.append(('random', edge_set, (a,b), cells))
                    break
            if found:
                break

        if found:
            print(f"  FOUND: random in {total_attempts} attempts")
        else:
            print(f"  Result: 0 successes")

    # ── Report ────────────────────────────────────────────────────
    dt = time.time() - t0
    print(f"\n{'='*55}")
    print(f"Total: {total_attempts} attempts, {dt:.1f}s, solutions={len(found)}")

    if found:
        surgery_type, details, split_edge, cells = found[0]
        print(f"\n✅ m=37 SOLUTION FOUND via {surgery_type}")
        
        result = {
            'm': 37, 'surgery': surgery_type,
            'switched_edges': [[e[0],e[1]] for e in details] if details else [],
            'split_edge': list(split_edge) if split_edge else [],
            'cells': sorted(cells),
            'n_cells': len(cells),
            'verify_no_collinear': check_no_collinear(cells, 37),
            'attempts': total_attempts, 'time_sec': dt,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        with open('results/m37_surgery_result.json', 'w') as f:
            json.dump(result, f, indent=1)
        with open('results/m37_surgery_found.txt', 'w') as f:
            f.write(f"m=37 rot4 NTIL solution via direction 2 surgery\n")
            f.write(f"Surgery: {surgery_type}\n")
            f.write(f"Switched edges: {details}\nSplit edge: {split_edge}\n\nCells:\n")
            for (x,y) in sorted(cells):
                f.write(f"  ({x}, {y})\n")
            f.write(f"\nVerify: {check_no_collinear(cells, 37)}\n")
        print(f"Saved to results/m37_surgery_result.json")
    else:
        print(f"\n❌ No m=37 solution found via surgery.")
        with open('results/m37_surgery_log.txt', 'w') as f:
            f.write(f"m=36→37 surgery v3\nAttempts: {total_attempts}\nTime: {dt:.1f}s\n")
            f.write(f"Result: none found\n")
            f.write(f"Methods: split, 2switch+split, 3switch+split, rand\n")


if __name__ == '__main__':
    main()
