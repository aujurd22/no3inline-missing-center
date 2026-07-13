"""
surgery_m36_to_37_v4.py — Direction 2: OPTIMISED correct-N surgery

Fix: ALL points computed with N=74. Incremental check: only 8 new points vs existing.
~82K triple checks per split attempt vs ~530K for full brute. ~6× faster.
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

def line_of(p, q):
    x1, y1 = p; x2, y2 = q
    A, B, C = y2 - y1, x1 - x2, x2 * y1 - x1 * y2
    g = math.gcd(math.gcd(A, B), C)
    if g: A//=g; B//=g; C//=g
    if A < 0 or (A == 0 and B < 0): A, B, C = -A, -B, -C
    return (A, B, C)

def build_line_pts(pts):
    """dict: line_key -> set of indices. Only keeps lines with ≥2 points."""
    n = len(pts)
    lp = {}
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            if xi == xj and yi == yj: continue
            k = line_of((xi, yi), (xj, yj))
            s = lp.get(k)
            if s is None:
                lp[k] = {i, j}
            else:
                s.add(i); s.add(j)
    # Remove single-point entries
    return {k: v for k, v in lp.items() if len(v) >= 3}

def points_from_cells(cells, m):
    """C4-lift all cells with N=2m. Returns list of (x,y)."""
    N = 2 * m
    pts = []
    for (x, y) in cells:
        for r in range(4):
            pts.append(c4((x, y), r, N))
    return pts

def check_new_points_no_collinear(pts_all, line_pts_old, old_n, new_indices):
    """
    Check collinearity involving new_indices.
    Returns True if NO 3-collinear (safe), else False.
    Checks: (new, new, new), (new, old, old), (new, new, old).
    """
    # Check all pairs involving at least one new point
    n = len(pts_all)
    new_set = set(new_indices)
    
    for i in new_indices:
        xi, yi = pts_all[i]
        for j in range(n):
            if j <= i: continue
            xj, yj = pts_all[j]
            if xi == xj and yi == yj: continue
            k = line_of((xi, yi), (xj, yj))
            # Get all points on this line
            S = line_pts_old.get(k)
            if S is None:
                S = set()
            # Count points on this line
            on_line = len(S)
            if i >= old_n: on_line += 1
            if j >= old_n: on_line += 1
            # Also check if any OLD point we missed is on this line
            # (only needed if both i,j are new — line may have old pts)
            if i >= old_n and j >= old_n and S:
                # At least 2 old points on this line already
                if len(S) >= 1:
                    return False  # new+new+old = collinear
            if on_line >= 3:
                return False
    return True

def check_cells_no_collinear(cells, m):
    """Full brute force check. Returns True if safe."""
    pts = points_from_cells(cells, m)
    n = len(pts)
    lp = build_line_pts(pts)
    return len(lp) == 0  # no lines with ≥3 points

def cells_from_sel(sel, m):
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]

def verify_2factor(sel, m):
    for i in range(m):
        if sum(sel[i]) + sum(sel[r][i] for r in range(m)) != 2:
            return False
    return True

def get_edges(sel, m):
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]

def apply_k_switch(sel_old, edge_set, bijection):
    sel = [row[:] for row in sel_old]
    for (i, j) in edge_set: sel[i][j] = 0
    for (src, tgt) in bijection:
        if sel[src][tgt] == 1: return None
        sel[src][tgt] = 1
    return sel

def try_split_incremental(sel_mid, cells_mid, line_pts_mid, pts_mid, m_old, split_edge):
    """Try split, incremental check on 8 new points. Returns (ok, cells) or (False,None)."""
    m_new = m_old + 1
    a, b = split_edge
    # Build selection matrix
    sel = [row[:] for row in sel_mid]
    sel.append([0] * m_new)
    for r in range(m_new): sel[r].append(0)
    if sel[a][b] != 1: return False, None
    sel[a][b] = 0
    sel[a][m_old] = 1
    sel[m_old][b] = 1
    if not verify_2factor(sel, m_new): return False, None
    
    cells = cells_from_sel(sel, m_new)
    if len(cells) != m_new: return False, None
    
    # New cells: (a, m_old) and (m_old, b)
    old_n = len(pts_mid)
    pts_all = list(pts_mid)
    N = 2 * m_new
    new_cells = [(a, m_old), (m_old, b)]
    new_idx = []
    for (cx, cy) in new_cells:
        idx0 = len(pts_all)
        for r in range(4):
            pts_all.append(c4((cx, cy), r, N))
        new_idx.extend(range(idx0, len(pts_all)))
    
    # Incremental check
    if check_new_points_no_collinear(pts_all, line_pts_mid, old_n, new_idx):
        return True, cells
    return False, None

def generate_vd_pairs(edges):
    n = len(edges)
    for i in range(n):
        ei = edges[i]; vi = {ei[0], ei[1]}
        for j in range(i+1, n):
            ej = edges[j]
            if ej[0] in vi or ej[1] in vi: continue
            yield [ei, ej]

def generate_vd_triples(edges):
    n = len(edges)
    for i in range(n):
        ei = edges[i]; vi = {ei[0], ei[1]}
        for j in range(i+1, n):
            ej = edges[j]
            if ej[0] in vi or ej[1] in vi: continue
            vj = vi | {ej[0], ej[1]}
            for k in range(j+1, n):
                ek = edges[k]
                if ek[0] in vj or ek[1] in vj: continue
                yield [ei, ej, ek]

def surgery_main():
    m_old = 36
    m_new = 37
    
    # Load seed
    line = open('flammenkamp_cache/n72_rot4.few').read().strip()
    pts = decode_line(line, 72)
    cells_old = sorted([(x,y) for (x,y) in pts if x < m_old and y < m_old])
    sel_old = [[0]*m_old for _ in range(m_old)]
    for (x,y) in cells_old: sel_old[x][y] = 1
    assert verify_2factor(sel_old, m_old)
    edges_old = get_edges(sel_old, m_old)
    print(f"[v4] m=36: {len(edges_old)} edges, seed OK", flush=True)
    
    total = 0
    found = []
    t0 = time.time()
    
    # ── Phase 1: simple split ────────────────────────────────────
    # Lift m=36 cells with N=74 (not N=72)
    pts36 = points_from_cells(cells_old, 37)  # N=2*37=74
    lp36 = build_line_pts(pts36)
    assert len(lp36) == 0, "m=36 cells lifted with N=74 should have 0 collinear"
    
    for (a,b) in edges_old:
        total += 1
        ok, cells = try_split_incremental(sel_old, cells_old, lp36, pts36, m_old, (a,b))
        if ok:
            found.append(('split', (a,b), cells))
            break
    if found:
        print(f"Phase 1: FOUND split {found[0][1]}! (attempts={total})", flush=True)
    else:
        print(f"Phase 1: 0/{len(edges_old)} (attempts={total})", flush=True)
    
    # ── Phase 2: 2-switch + split ────────────────────────────────
    if not found:
        print(f"Phase 2: 2-switch+split...", flush=True)
        n2 = 0
        for ep in generate_vd_pairs(edges_old):
            e1, e2 = ep
            sources = [e1[0], e2[0]]
            targets = [e1[1], e2[1]]
            for perm in permutations(targets):
                bi = list(zip(sources, perm))
                if bi[0][1] == targets[0] and bi[1][1] == targets[1]:
                    continue
                sm = apply_k_switch(sel_old, ep, bi)
                if sm is None: continue
                if not verify_2factor(sm, m_old): continue
                cm = cells_from_sel(sm, m_old)
                pm = points_from_cells(cm, 37)
                lm = build_line_pts(pm)
                if lm: continue  # switch broke collinearity
                
                em = get_edges(sm, m_old)
                for (a,b) in em:
                    total += 1
                    ok, cells = try_split_incremental(sm, cm, lm, pm, m_old, (a,b))
                    if ok:
                        found.append(('2switch+split', ep, cells))
                        break
                if found: break
            if found: break
            n2 += 1
            if n2 % 100 == 0:
                print(f"  2-switch: {n2} checked, {total} total attempts", flush=True)
        if found:
            print(f"Phase 2: FOUND! {found[0][1]}", flush=True)
        else:
            print(f"Phase 2: 0 after {n2} switches, {total} attempts", flush=True)
    
    # ── Phase 3: 3-switch + split ────────────────────────────────
    if not found:
        print(f"Phase 3: 3-switch+split...", flush=True)
        n3 = 0
        for et in generate_vd_triples(edges_old):
            sources = [e[0] for e in et]
            targets = [e[1] for e in et]
            for perm in permutations(targets):
                bi = list(zip(sources, perm))
                if all(bi[i][1] == targets[i] for i in range(3)):
                    continue
                sm = apply_k_switch(sel_old, et, bi)
                if sm is None: continue
                if not verify_2factor(sm, m_old): continue
                cm = cells_from_sel(sm, m_old)
                pm = points_from_cells(cm, 37)
                lm = build_line_pts(pm)
                if lm: continue
                
                em = get_edges(sm, m_old)
                for (a,b) in em:
                    total += 1
                    ok, cells_proto = try_split_incremental(sm, cm, lm, pm, m_old, (a,b))
                    if ok:
                        found.append(('3switch+split', et, cells_proto))
                        break
                if found: break
            if found: break
            n3 += 1
            if n3 % 200 == 0:
                print(f"  3-switch: {n3} checked, {total} total attempts", flush=True)
        if found:
            print(f"Phase 3: FOUND! {found[0][1]}", flush=True)
        else:
            print(f"Phase 3: 0 after {n3} triples, {total} attempts", flush=True)
    
    # ── Report ────────────────────────────────────────────────────
    dt = time.time() - t0
    print(f"\nTotal: {total} attempts, {dt:.1f}s, solutions={len(found)}", flush=True)
    
    if found:
        st, details, cells = found[0]
        print(f"\n✅ m=37 SOLUTION via {st}:", flush=True)
        for (x,y) in sorted(cells):
            print(f"  ({x}, {y})", flush=True)
        
        # Double verify with full brute
        v = check_cells_no_collinear(cells, m_new)
        print(f"Full brute verify: {v} (should be True)", flush=True)
        
        result = {
            'm': 37, 'surgery': st,
            'details': str(details),
            'cells': sorted(cells),
            'verify_no_collinear': v,
            'attempts': total, 'time_sec': round(dt, 1),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        json.dump(result, open('results/m37_surgery_result.json', 'w'), indent=1)
        with open('results/m37_surgery_found.txt', 'w') as f:
            f.write(f"m=37 rot4 NTIL via {st}\n{details}\n\n")
            for (x,y) in sorted(cells): f.write(f"({x},{y})\n")
            f.write(f"\nverify={v}\n")
        print(f"Saved to results/", flush=True)
    else:
        print(f"\n❌ No m=37 solution found via surgery.", flush=True)
        json.dump({'attempts': total, 'time_sec': dt, 'found': False,
                   'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')},
                  open('results/m37_surgery_log.txt', 'w'))

if __name__ == '__main__':
    surgery_main()
