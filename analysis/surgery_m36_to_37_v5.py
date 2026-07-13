"""
surgery_m36_to_37_v5.py — Direction 2: CORRECT incremental collinearity check

After k-switch: 36 cells (144 pts, N=74, verified safe).
After split: 35 common cells (140 old pts) + 2 new cells (8 new pts) = 148 pts total.
Check only triples involving ≥1 new point: ~1200 pair-of-points checks, <1ms total.
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
    g = math.gcd(math.gcd(abs(A), abs(B)), abs(C))
    if g: A //= g; B //= g; C //= g
    if A < 0 or (A == 0 and B < 0): A, B, C = -A, -B, -C
    return (A, B, C)

def points_from_cells(cells, m):
    N = 2 * m
    pts = []
    for (x, y) in cells:
        for r in range(4):
            pts.append(c4((x, y), r, N))
    return pts

def check_no_collinear_line(cells, m):
    """Line-based check (O(N²) pairs, not O(N³)). Returns True if safe."""
    pts = points_from_cells(cells, m)
    lp = build_line_pts(pts)
    return len(lp) == 0

def check_no_collinear_brute(cells, m):
    """Full brute force (O(N³)). For backward compat / final verify."""
    return check_no_collinear_line(cells, m)

def build_line_count(pts):
    """Map line_key → count of points on that line (for given pts list)."""
    n = len(pts)
    lc = {}
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            k = line_of((xi, yi), (xj, yj))
            # We only count unique points on this line
            lc[k] = lc.get(k, 0) + 1  # count of PAIRS, not points
    return lc

# Wait — line count of pairs doesn't work. C(t,2) pairs from t points.
# We need: does this line have ≥3 points? NOT how many pairs.
# Let me use a set-based approach but only for lines that matter.

def build_line_pts(pts):
    """Map line_key → set of point indices on that line (≥3)."""
    n = len(pts)
    lp = {}
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            if xi == xj and yi == yj: continue
            k = line_of((xi, yi), (xj, yj))
            s = lp.get(k)
            if s is None: lp[k] = {i, j}
            else: s.add(i); s.add(j)
    return {k: v for k, v in lp.items() if len(v) >= 3}

def safe_split_check(pts_old, new_pts, old_start, n_new):
    """
    Incremental collinearity check.
    pts_old: list of (x,y) for old points (safe among themselves)
    new_pts: list of (x,y) for new points
    old_start: index in combined array where old points start
    n_new: number of new points
    
    Strategy: build line_count for old points only (set-based → exact count),
    then check all new+old and new+new+old triples.
    """
    n_old = len(pts_old)
    
    # Build line→count for OLD points only (exact point count per line)
    old_lc = {}
    for i in range(n_old):
        xi, yi = pts_old[i]
        for j in range(i + 1, n_old):
            xj, yj = pts_old[j]
            k = line_of((xi, yi), (xj, yj))
            s = old_lc.get(k)
            if s is None: old_lc[k] = {i, j}
            else: s.add(i); s.add(j)
    # Convert to count
    old_count = {k: len(v) for k, v in old_lc.items()}
    
    # Check (new, old, old): for each new point, check lines with old points
    for ni in range(n_new):
        nxi, nyi = new_pts[ni]
        for oj in range(n_old):
            oxj, oyj = pts_old[oj]
            k = line_of((nxi, nyi), (oxj, oyj))
            # How many old points are on this line?
            cnt = old_count.get(k, 0)
            if cnt >= 2:  # ni + 2 old = 3 collinear
                return False
    
    # Check (new, new, old) and (new, new, new): for each pair of new points
    for ni in range(n_new):
        nxi, nyi = new_pts[ni]
        for nj in range(ni + 1, n_new):
            nxj, nyj = new_pts[nj]
            k = line_of((nxi, nyi), (nxj, nyj))
            cnt = old_count.get(k, 0)
            if cnt >= 1:  # ni + nj + 1 old = 3 collinear
                return False
            # Check if a third new point is also on this line
            # (rare, but possible: all 3 new points collinear)
            for nk in range(nj + 1, n_new):
                nxk, nyk = new_pts[nk]
                if (nxj - nxi) * (nyk - nyi) == (nyj - nyi) * (nxk - nxi):
                    return False  # 3 new points collinear
    
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

def apply_k_switch(sel_old, edge_set, bijection):
    sel = [row[:] for row in sel_old]
    for (i, j) in edge_set: sel[i][j] = 0
    for (src, tgt) in bijection:
        if sel[src][tgt] == 1: return None
        sel[src][tgt] = 1
    return sel

def try_split_m37(sel_mid, m_old, split_edge):
    """
    Try splitting split_edge=(a,b) by inserting vertex m_old.
    Returns (ok, cells, pts_all) for a valid m=37 solution.
    """
    m_new = m_old + 1
    a, b = split_edge
    
    # Build new selection matrix  
    sel = [row[:] for row in sel_mid]
    sel.append([0] * m_new)
    for r in range(m_new): sel[r].append(0)
    
    if sel[a][b] != 1: return False, None, None
    sel[a][b] = 0
    sel[a][m_old] = 1
    sel[m_old][b] = 1
    if not verify_2factor(sel, m_new): return False, None, None
    
    cells = cells_from_sel(sel, m_new)
    if len(cells) != m_new: return False, None, None
    
    # Lift cells with N=2*m_new
    N = 2 * m_new
    # Gather all 37 cells into a list of 148 points
    pts_all = []
    cell_to_pt_start = {}
    idx = 0
    for (cx, cy) in cells:
        cell_to_pt_start[(cx, cy)] = idx
        for r in range(4):
            pts_all.append(c4((cx, cy), r, N))
        idx += 4
    
    # The old cell (a, b) was removed. The new cells are (a, m_old) and (m_old, b).
    # So the new points are: the 8 pts from these 2 new cells.
    # We need to build the "old points" set: pts_all minus the 8 new pts.
    new_cells = {(a, m_old), (m_old, b)}
    new_indices = set()
    for (cx, cy), start in cell_to_pt_start.items():
        if (cx, cy) in new_cells:
            new_indices.update(range(start, start + 4))
    
    # Build old pts list
    old_pts = [pts_all[i] for i in range(len(pts_all)) if i not in new_indices]
    new_pts = [pts_all[i] for i in sorted(new_indices)]
    
    # Incremental check
    if safe_split_check(old_pts, new_pts, 0, len(new_pts)):
        return True, cells, pts_all
    return False, None, None


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
    print(f"[v5] m=36: {len(edges_old)} edges, full verify: {check_no_collinear_brute(cells_old, 37)}", flush=True)
    
    total = 0
    found = []
    t0 = time.time()
    
    # ── Phase 1: simple split ────────────────────────────────────
    print(f"Phase 1: split (0/{len(edges_old)})", flush=True)
    for idx, (a,b) in enumerate(edges_old):
        total += 1
        ok, cells, _ = try_split_m37(sel_old, m_old, (a,b))
        if ok:
            found.append(('split', (a,b), cells))
            print(f"  FOUND: split {a}→{b}!", flush=True)
            break
        elif (idx+1) % 12 == 0:
            print(f"  {idx+1}/{len(edges_old)}...", flush=True)
    
    if not found: print(f"Phase 1: 0/{len(edges_old)}", flush=True)
    
    # ── Phase 2: 2-switch + split ────────────────────────────────
    if not found:
        print(f"Phase 2: 2-switch+split...", flush=True)
        n2 = 0
        for ep in generate_vd_pairs(edges_old):
            sources = [ep[0][0], ep[1][0]]
            targets = [ep[0][1], ep[1][1]]
            for perm in permutations(targets):
                bi = list(zip(sources, perm))
                if bi[0][1] == targets[0] and bi[1][1] == targets[1]:
                    continue
                sm = apply_k_switch(sel_old, ep, bi)
                if sm is None: continue
                if not verify_2factor(sm, m_old): continue
                cm = cells_from_sel(sm, m_old)
                if not check_no_collinear_brute(cm, 37): continue
                
                em = get_edges(sm, m_old)
                for (a,b) in em:
                    total += 1
                    ok, cells, _ = try_split_m37(sm, m_old, (a,b))
                    if ok:
                        found.append(('2switch+split', ep, cells))
                        break
                if found: break
            if found: break
            n2 += 1
            if n2 % 200 == 0:
                print(f"  2-switch: {n2}, total att={total}", flush=True)
        if found:
            print(f"Phase 2: FOUND! att={total}", flush=True)
        else:
            print(f"Phase 2: 0 after {n2} switches, {total} att", flush=True)
    
    # ── Phase 3: 3-switch + split ────────────────────────────────
    if not found:
        print(f"Phase 3: 3-switch+split...", flush=True)
        n3 = 0
        for et in generate_vd_triples(edges_old):
            sources = [e[0] for e in et]
            targets = [e[1] for e in et]
            for perm in permutations(targets):
                bi = list(zip(sources, perm))
                if all(bi[i][1]==targets[i] for i in range(3)): continue
                sm = apply_k_switch(sel_old, et, bi)
                if sm is None: continue
                if not verify_2factor(sm, m_old): continue
                cm = cells_from_sel(sm, m_old)
                if not check_no_collinear_brute(cm, 37): continue
                
                em = get_edges(sm, m_old)
                for (a,b) in em:
                    total += 1
                    ok, cells, _ = try_split_m37(sm, m_old, (a,b))
                    if ok:
                        found.append(('3switch+split', et, cells))
                        break
                if found: break
            if found: break
            n3 += 1
            if n3 % 500 == 0:
                print(f"  3-switch: {n3}, total att={total}", flush=True)
        if found:
            print(f"Phase 3: FOUND! att={total}", flush=True)
        else:
            print(f"Phase 3: 0 after {n3} triples, {total} att", flush=True)
    
    # ── Phase 4: random k-switch + split ─────────────────────────
    if not found:
        print(f"Phase 4: random k-switch (50000 attempts)...", flush=True)
        for _ in range(50000):
            k = rng.randint(2, 4)
            shuffled = list(edges_old)
            rng.shuffle(shuffled)
            edge_set = []; used_v = set()
            for e in shuffled:
                if e[0] not in used_v and e[1] not in used_v:
                    edge_set.append(e); used_v.add(e[0]); used_v.add(e[1])
                    if len(edge_set) == k: break
            if len(edge_set) < k: continue
            sources = [e[0] for e in edge_set]
            targets = [e[1] for e in edge_set]
            perm = list(targets); rng.shuffle(perm)
            bi = list(zip(sources, perm))
            sm = apply_k_switch(sel_old, edge_set, bi)
            if sm is None: continue
            if not verify_2factor(sm, m_old): continue
            cm = cells_from_sel(sm, m_old)
            if not check_no_collinear_brute(cm, 37): continue
            em = get_edges(sm, m_old)
            rng.shuffle(em)
            for (a,b) in em[:50]:
                total += 1
                ok, cells, _ = try_split_m37(sm, m_old, (a,b))
                if ok:
                    found.append(('random', edge_set, cells))
                    break
            if found: break
        if found:
            print(f"Phase 4: FOUND! att={total}", flush=True)
        else:
            print(f"Phase 4: 0 after {total} attempts", flush=True)
    
    # ── Report ────────────────────────────────────────────────────
    dt = time.time() - t0
    print(f"\nTotal: {total} attempts, {dt:.1f}s, solutions={len(found)}", flush=True)
    
    if found:
        st, details, cells = found[0]
        v = check_no_collinear_brute(cells, m_new)
        print(f"\n✅ m=37  via {st} (verify={v})", flush=True)
        for (x,y) in sorted(cells):
            print(f"  ({x},{y})", flush=True)
        
        json.dump({'m':37,'surgery':st,'details':str(details),
                   'cells':sorted(cells),'verify':v,
                   'attempts':total,'time_sec':round(dt,1),
                   'ts':time.strftime('%Y-%m-%d %H:%M:%S')},
                  open('results/m37_surgery_result.json','w'),indent=1)
        with open('results/m37_surgery_found.txt','w') as f:
            f.write(f"m=37 rot4 NTIL via {st}\n{details}\n\n")
            for (x,y) in sorted(cells): f.write(f"({x},{y})\n")
            f.write(f"\nverify={v}\n")
        print("Saved!", flush=True)
    else:
        print(f"\n❌ No m=37 solution via surgery.", flush=True)
        json.dump({'attempts':total,'time_sec':dt,'found':False,
                   'ts':time.strftime('%Y-%m-%d %H:%M:%S')},
                  open('results/m37_surgery_log.txt','w'))

if __name__ == '__main__':
    surgery_main()
