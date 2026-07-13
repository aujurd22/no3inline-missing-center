"""
surgery_m36_to_37.py — Direction 2: m=36 → 37 有界2-因子手术

Strategy:
  1. Load the known m=36 rot4 solution (single 36-cycle).
  2. Represent as directed 2-factor sel[m][m] where sel[i][j]=1 means "i→j".
  3. For each edge (i→j), try a 1-edge split: replace (i→j) with (i→36)+(36→j).
     This adds vertex 37 with degree 2 and preserves all other degrees.
  4. For each resulting m=37 2-factor, C4-lift → 4m=148 points → brute-check collinearity.
  5. If found, save. If not, escalate to k-edge switch + split.

Usage:  python surgery_m36_to_37.py
"""

import sys, os, json, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ── Geometry primitives (same as quadratic_sidon_completeness.py) ──────
def c4(p, r, N):
    x, y = p
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    if r == 3: return (y, N - 1 - x)

def lift_cells(cells, m):
    """C4-lift fundamental quadrant cells (x,y) to 4m points on N×N grid, N=2m."""
    N = 2 * m
    pts = []
    for (x, y) in cells:
        for r in range(4):
            pts.append(c4((x, y), r, N))
    return pts

def brute_collinear(pts):
    """Check if ANY 3 of pts are collinear (integer cross-product)."""
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
                    return True  # collinear triple found
    return False

def cells_from_sel(sel, m):
    """Extract list of (x,y) cells from selection matrix."""
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]

def sel_from_cells(cells, m):
    """Build selection matrix from cell list."""
    sel = [[0] * m for _ in range(m)]
    for (i, j) in cells:
        sel[i][j] = 1
    return sel

def verify_2factor(sel, m):
    """Check rowSum[i]+colSum[i]==2 for all i."""
    for i in range(m):
        rs = sum(sel[i])
        cs = sum(sel[r][i] for r in range(m))
        if rs + cs != 2:
            return False
    return True


# ══════════════════════════════════════════════════════════════════════
#  SURGERY PRIMITIVES
# ══════════════════════════════════════════════════════════════════════

def split_edge(sel_old, m_old, edge):
    """
    Split an edge (i→j) by inserting new vertex at index m_old.
    Returns the new selection matrix for m = m_old + 1.
    """
    i, j = edge
    m_new = m_old + 1
    sel = [[0] * m_new for _ in range(m_new)]
    # Copy old edges
    for a in range(m_old):
        for b in range(m_old):
            sel[a][b] = sel_old[a][b]
    # Remove old edge i→j
    sel[i][j] = 0
    # Add i→m_old and m_old→j
    sel[i][m_old] = 1
    sel[m_old][j] = 1
    return sel


def two_switch_and_split(sel_old, m_old, e1, e2, split_target):
    """
    Perform a 2-switch on edges e1=(a,b), e2=(c,d):
      Remove (a,b), (c,d)
      Add (a,c), (b,d)             -- 2-switch (standard)
    
    Then split split_target by inserting vertex m_old.
    split_target is one of 'ac', 'bd' (the two new edges after switch).
    The old edges e1,e2 are given as tuples (i,j).
    """
    a, b = e1
    c, d = e2
    m_new = m_old + 1
    sel = [[0] * m_new for _ in range(m_new)]
    
    # Copy all old edges
    for i in range(m_old):
        for j in range(m_old):
            sel[i][j] = sel_old[i][j]
    
    # Remove old edges
    sel[a][b] = 0
    sel[c][d] = 0
    
    # Add switched edges
    if split_target == 'ac':
        # 2-switch: (a,b)+(c,d)→(a,c)+(b,d), then split (a,c) into (a,m)+(m,c)
        sel[a][c] = 1
        sel[b][d] = 1
        # Now split (a,c):
        sel[a][c] = 0
        sel[a][m_old] = 1
        sel[m_old][c] = 1
    elif split_target == 'bd':
        # 2-switch: (a,b)+(c,d)→(a,c)+(b,d), then split (b,d)
        sel[a][c] = 1
        sel[b][d] = 1
        sel[b][d] = 0
        sel[b][m_old] = 1
        sel[m_old][d] = 1
    elif split_target == 'ca':
        # 2-switch alt: (a,b)+(c,d)→(c,a)+(d,b), then split (c,a)
        sel[c][a] = 1
        sel[d][b] = 1
        sel[c][a] = 0
        sel[c][m_old] = 1
        sel[m_old][a] = 1
    elif split_target == 'db':
        # 2-switch alt: (a,b)+(c,d)→(c,a)+(d,b), then split (d,b)
        sel[c][a] = 1
        sel[d][b] = 1
        sel[d][b] = 0
        sel[d][m_old] = 1
        sel[m_old][b] = 1
    else:
        raise ValueError(f"Unknown split_target: {split_target}")
    
    return sel


# ══════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════

def main():
    m_old = 36
    m_new = 37
    
    # 1. Load m=36 solution
    line = open('flammenkamp_cache/n72_rot4.few').read().strip()
    from rot4_loader import decode_line
    pts = decode_line(line, 72)
    cells_old = [(x, y) for (x, y) in pts if x < m_old and y < m_old]
    cells_old.sort()
    sel_old = sel_from_cells(cells_old, m_old)
    
    assert len(cells_old) == 36
    assert verify_2factor(sel_old, m_old)
    
    edges_old = [(i, j) for i in range(m_old) for j in range(m_old) if sel_old[i][j]]
    print(f"[surgery] Loaded m=36 solution: {len(edges_old)} edges, single 36-cycle")
    
    # ── Phase 1: Simple edge split (36 attempts) ────────────────────
    print(f"\n{'='*60}")
    print(f"Phase 1: 1-edge split (insert vertex 37 into each existing edge)")
    print(f"{'='*60}")
    t0 = time.time()
    found = []
    for idx, (i, j) in enumerate(edges_old):
        sel_new = split_edge(sel_old, m_old, (i, j))
        assert verify_2factor(sel_new, m_new), f"2-factor broken for edge {i}→{j}"
        
        cells_new = cells_from_sel(sel_new, m_new)
        pts_new = lift_cells(cells_new, m_new)
        coll = brute_collinear(pts_new)
        if not coll:
            found.append((i, j, 'split', cells_new))
            print(f"  ✓ edge {i}→{j}: NO 3-collinear! SOLVED!")
        else:
            print(f"  ✗ edge {i}→{j}: collinear", end='')
            if (idx+1) % 12 == 0:
                print()
    
    dt = time.time() - t0
    print(f"\nPhase 1: {len(found)}/{len(edges_old)} successes in {dt:.1f}s")
    
    if found:
        _save_and_exit(found, m_new)
        return
    
    # ── Phase 2: 2-edge switch + split ─────────────────────────────
    print(f"\n{'='*60}")
    print(f"Phase 2: 2-edge switch + split")
    print(f"{'='*60}")
    t0 = time.time()
    
    # Generate edge pairs (order does not matter for unordered pairs)
    npairs = 0
    for ei in range(len(edges_old)):
        for ej in range(ei + 1, len(edges_old)):
            e1 = edges_old[ei]
            e2 = edges_old[ej]
            # Disallow self-pair and adjacent pairs (share a vertex)
            if len(set(e1) & set(e2)) > 0:
                continue  # skip edges sharing a vertex
            for split_target in ['ac', 'bd', 'ca', 'db']:
                try:
                    sel_new = two_switch_and_split(sel_old, m_old, e1, e2, split_target)
                except (AssertionError, ValueError):
                    continue
                
                if not verify_2factor(sel_new, m_new):
                    continue
                
                cells_new = cells_from_sel(sel_new, m_new)
                pts_new = lift_cells(cells_new, m_new)
                coll = brute_collinear(pts_new)
                npairs += 1
                if not coll:
                    found.append((e1, e2, split_target, cells_new))
                    print(f"  ✓ switch({e1},{e2})+split({split_target}): SOLVED!")
    
    dt = time.time() - t0
    print(f"\nPhase 2: {len(found)} successes out of {npairs} attempts in {dt:.1f}s")
    
    if found:
        _save_and_exit(found, m_new)
        return
    
    print(f"\n{'='*60}")
    print(f"RESULT: No surgery found a valid m=37 solution.")
    print(f"Directions not yet tried: 3-edge switch, longer alternating cycles.")
    print(f"{'='*60}")


def _save_and_exit(found, m_new):
    """Save the first found solution."""
    cells = found[0][3]
    result = {
        'm': m_new,
        'surgery_type': str(found[0][:3]),
        'cells': cells,
        'n_cells': len(cells),
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    # Also compute a compressed summary
    N = 2 * m_new
    pts = lift_cells(cells, m_new)
    result['n_points'] = len(pts)
    result['verify_no_collinear'] = not brute_collinear(pts)
    
    path = 'results/m37_surgery_result.json'
    with open(path, 'w') as f:
        json.dump(result, f, indent=1)
    # Also as simple text file for human reading
    with open('results/m37_surgery_result.txt', 'w') as f:
        f.write(f"m=37 rot4 NTIL solution via direction 2 surgery\n")
        f.write(f"Surgery: {found[0][:3]}\n")
        f.write(f"Coefficients (top-left quadrant, {len(cells)} cells):\n")
        for (x,y) in cells:
            f.write(f"  ({x}, {y})\n")
        f.write(f"\nVerify: collinear={result['verify_no_collinear']}\n")
    
    print(f"\nSaved to results/m37_surgery_result.json + .txt")


if __name__ == '__main__':
    main()
