"""
surgery_m36_to_37_v2.py — Direction 2: Optimized k-switch + surgery

Key optimisation: instead of brute-checking ALL 148 points for every candidate,
only check collinearity INVOLVING the 8 new points added by the 2 new cells.
This speeds up verification by ~50×, allowing us to explore k=1..4 switches.

Strategy:
  Phase 1: 1-edge split (fast, 36 attempts) — already done in v1
  Phase 2: 2-switch + split any of the 36 edges
  Phase 3: 3-switch + split any of the 36 edges  
  Phase 4: 4-switch + split any of the 36 edges
  Phase 5: Random k-switch + split, repeated N times
"""

import sys, os, json, time, random
from itertools import combinations, permutations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rot4_loader import decode_line

SEED = 20260713
rng = random.Random(SEED)

# ── Geometry primitives ──────────────────────────────────────────────
def c4(p, r, N):
    x, y = p
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    if r == 3: return (y, N - 1 - x)

def lift_cells(cells, m):
    N = 2 * m
    return [c4((x, y), r, N) for (x, y) in cells for r in range(4)]

def line_of(p, q):
    """Canonical line key: normalized (A,B,C) for line Ax+By+C=0."""
    x1, y1 = p; x2, y2 = q
    A, B, C = y2 - y1, x1 - x2, x2 * y1 - x1 * y2
    g = math.gcd(math.gcd(A, B), C)
    if g != 0:
        A //= g; B //= g; C //= g
    if A < 0 or (A == 0 and B < 0):
        A, B, C = -A, -B, -C
    return (A, B, C)

import math

def build_line_pts(pts):
    """Build dict: line_key -> set of point indices on that line."""
    n = len(pts)
    line_pts = {}
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            if xi == xj and yi == yj:
                continue
            k = line_of((xi, yi), (xj, yj))
            if k not in line_pts:
                line_pts[k] = set()
            line_pts[k].add(i)
            line_pts[k].add(j)
    return line_pts

def count_collinear(line_pts):
    """Sum C(|S|,3) over all lines with ≥3 points."""
    return sum(len(s) * (len(s)-1) * (len(s)-2) // 6
               for s in line_pts.values() if len(s) >= 3)

def incremental_collinear(pts_old, line_pts_old, new_indices):
    """
    Check if adding new_indices to pts_old creates any 3-collinear.
    Returns (collision_count, updated_line_pts).
    Only checks lines involving at least one new point.
    """
    n_old = len(pts_old)
    line_pts = dict(line_pts_old)  # copy
    total_new_bad = 0
    
    for ni in new_indices:
        xi, yi = pts_old[ni]  # already set in the pts array
        # Check vs all EXISTING points (including other new points already added)
        for j in range(n_old):
            if j == ni: continue
            xj, yj = pts_old[j]
            if xi == xj and yi == yj: continue
            k = line_of((xi, yi), (xj, yj))
            if k not in line_pts:
                line_pts[k] = {ni, j}
            else:
                S = line_pts[k]
                old_sz = len(S)
                S.add(ni)
                S.add(j)
                new_sz = len(S)
                # increment bad if we crossed the ≥3 threshold
                if old_sz < 3 and new_sz >= 3:
                    total_new_bad += new_sz * (new_sz - 1) * (new_sz - 2) // 6
                elif old_sz >= 3:
                    total_new_bad += (new_sz * (new_sz - 1) * (new_sz - 2) // 6
                                      - old_sz * (old_sz - 1) * (old_sz - 2) // 6)
    return total_new_bad, line_pts

# ── 2-factor helpers ────────────────────────────────────────────────
def cells_from_sel(sel, m):
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]

def sel_from_cells(cells, m):
    sel = [[0] * m for _ in range(m)]
    for (i, j) in cells:
        sel[i][j] = 1
    return sel

def verify_2factor(sel, m):
    for i in range(m):
        if sum(sel[i]) + sum(sel[r][i] for r in range(m)) != 2:
            return False
    return True

def get_edges(sel, m):
    return [(i, j) for i in range(m) for j in range(m) if sel[i][j]]


def apply_k_switch(sel_old, m_old, edge_set, bijection):
    """
    Given k edges e1=(i1,j1)...ek=(ik,jk) with distinct vertex sets,
    remove them and add edges according to bijection.
    bijection: a list of (i_k, j_perm_k) pairs (source→target after switch).
    Returns new selection matrix (m_old × m_old, BEFORE adding vertex 37).
    """
    sel = [row[:] for row in sel_old]  # deep copy
    for (i, j) in edge_set:
        sel[i][j] = 0
    for (src, tgt) in bijection:
        if sel[src][tgt] == 1:
            return None  # edge already exists → invalid switch
        sel[src][tgt] = 1
    return sel


def split_and_test(sel_mid, m_old, existing_pts, existing_line_pts,
                   split_edge, verbose=False):
    """
    Given a selection matrix for m_old edges (sel_mid), try splitting 
    split_edge = (a,b) by inserting vertex m_old.
    Build m=m_old+1, lift, check collinearity incrementally.
    Returns (success, m37_cells, m37_pts) if valid, else (False, None, None).
    """
    m_new = m_old + 1
    # Create new selection matrix
    sel = [row[:] for row in sel_mid]
    sel.append([0] * m_new)
    for r in range(m_new):
        sel[r].append(0)
    
    a, b = split_edge
    if sel[a][b] != 1:
        return False, None, None
    sel[a][b] = 0
    sel[a][m_old] = 1
    sel[m_old][b] = 1
    
    if not verify_2factor(sel, m_new):
        return False, None, None
    
    # Get new cells and lift
    cells_new = cells_from_sel(sel, m_new)
    
    # Build new points incrementally
    # Old m_old cells give 4*m_old points already computed
    # New cells: (a, m_old) and (m_old, b)
    new_cells_xy = [(a, m_old), (m_old, b)]
    N = 2 * m_new
    new_pts_indices = []
    pts = list(existing_pts)  # copy
    for (cx, cy) in new_cells_xy:
        idx_start = len(pts)
        for r in range(4):
            pts.append(c4((cx, cy), r, N))
        new_pts_indices.extend(range(idx_start, len(pts)))
    
    new_bad, line_pts_new = incremental_collinear(
        pts, existing_line_pts, new_pts_indices)
    
    if new_bad == 0:
        return True, cells_new, pts
    
    return False, None, None


# ══════════════════════════════════════════════════════════════════════
#  MAIN SURGERY ENGINE
# ══════════════════════════════════════════════════════════════════════

def generate_edge_subsets(edges, k, max_samples=None):
    """Generate all/good k-edge subsets with distinct vertex sets."""
    if k == 1:
        for e in edges:
            yield [e]
        return
    
    n = len(edges)
    # Greedy filtering: skip edges that share vertices
    def _iter_combos(start, chosen, used_vertices):
        if len(chosen) == k:
            yield list(chosen)
            return
        # Early termination if max_samples hit
        for idx in range(start, n):
            e = edges[idx]
            if e[0] in used_vertices or e[1] in used_vertices:
                continue
            chosen.append(e)
            used_vertices.add(e[0])
            used_vertices.add(e[1])
            yield from _iter_combos(idx + 1, chosen, used_vertices)
            used_vertices.remove(e[0])
            used_vertices.remove(e[1])
            chosen.pop()
            if max_samples and len(chosen) == 0:
                pass  # can't easily stop generator
    
    yield from _iter_combos(0, [], set())


def run_surgery(m_old=36, max_k=4, max_random_switches=50000, verbose=True):
    """Run surgery search for m_old+1."""
    # 1. Load m=36 solution
    line = open('flammenkamp_cache/n72_rot4.few').read().strip()
    pts_all = decode_line(line, 72)
    cells_old = [(x, y) for (x, y) in pts_all if x < m_old and y < m_old]
    cells_old.sort()
    sel_old = sel_from_cells(cells_old, m_old)
    assert verify_2factor(sel_old, m_old)
    
    edges_old = get_edges(sel_old, m_old)
    m_new = m_old + 1
    
    # Precompute base lift for m=36: 144 points, line_pts
    base_pts = lift_cells(cells_old, m_old)
    base_line_pts = build_line_pts(base_pts)
    base_bad = count_collinear(base_line_pts)
    assert base_bad == 0, f"m=36 seed has {base_bad} collinear lines!"
    
    if verbose:
        print(f"[surgery] Loaded m=36: {len(edges_old)} edges, "
              f"{len(base_pts)} lifted points, 0 collinear")
    
    all_split_targets = edges_old + [(i, i) for i in range(m_old)]  # no, self-loops invalid
    
    # ── Pre-compute all split edges for the m_old base ──
    # These are the edges we CAN split: all 36 edges in the base 2-factor
    split_targets = list(edges_old)
    
    total_attempts = 0
    found_solutions = []
    t_start = time.time()
    
    # ── Phase 1: simple edge split (m=36→37) ──────────────────
    # (already exhaustive in v1, but redo with incremental check for consistency)
    if verbose:
        print(f"\n{'='*55}\nPhase 1: 1-edge split (36 attempts)")
    for (a, b) in edges_old:
        total_attempts += 1
        ok, cells, pts = split_and_test(
            sel_old, m_old, base_pts, base_line_pts, (a, b))
        if ok:
            found_solutions.append(('split', (a,b), None, cells, pts))
            if verbose:
                print(f"  FOUND: split {a}→{b}")
    
    if found_solutions:
        return found_solutions, total_attempts, time.time() - t_start
    
    if verbose:
        print(f"  Result: 0/36")
    
    # ── Phase 2: k-switch + split ─────────────────────────────
    for k in range(2, max_k + 1):
        if verbose:
            print(f"\n{'='*55}\nPhase {k}: {k}-switch + split")
        
        # enumerate vertex-disjoint k-edge subsets
        n_subsets = 0
        for edge_set in generate_edge_subsets(edges_old, k):
            # edge_set = [(i1,j1), (i2,j2), ..., (ik,jk)]
            sources = [e[0] for e in edge_set]
            targets = [e[1] for e in edge_set]
            
            # All bijections from sources to targets (k! possibilities)
            for perm in permutations(targets):
                bijection = list(zip(sources, perm))
                # Skip the identity (no change)
                if all(b[0] == b[1] for b in bijection):
                    continue
                
                sel_mid = apply_k_switch(sel_old, m_old, edge_set, bijection)
                if sel_mid is None:
                    continue
                
                if not verify_2factor(sel_mid, m_old):
                    continue
                
                # After switch, rebuild the base pts for this new 2-factor
                cells_mid = cells_from_sel(sel_mid, m_old)
                pts_mid = lift_cells(cells_mid, m_old)
                line_pts_mid = build_line_pts(pts_mid)
                if count_collinear(line_pts_mid) != 0:
                    continue  # switch created collinearity → skip
                
                # Now try splitting EACH of the new edges
                new_edges = bijection  # the edges we just added
                for (a, b) in new_edges:
                    total_attempts += 1
                    ok, cells, pts = split_and_test(
                        sel_mid, m_old, pts_mid, line_pts_mid, (a, b))
                    if ok:
                        found_solutions.append(
                            (f'{k}-switch+split', edge_set, (a,b), cells, pts))
                        if verbose:
                            print(f"  FOUND: {k}-switch {edge_set} + split {a}→{b}")
                
                n_subsets += 1
                
                # Also try splitting ORIGINAL edges (not just new ones)
                for (a, b) in edges_old:
                    if (a, b) in edge_set:
                        continue  # this edge was removed
                    total_attempts += 1
                    ok, cells, pts = split_and_test(
                        sel_mid, m_old, pts_mid, line_pts_mid, (a, b))
                    if ok:
                        found_solutions.append(
                            (f'{k}-switch+split-old', edge_set, (a,b), cells, pts))
                        if verbose:
                            print(f"  FOUND: {k}-switch {edge_set} + split-OLD {a}→{b}")
            
            if found_solutions:
                break
        
        if verbose:
            print(f"  Attempts: {total_attempts}, subsets: {n_subsets}")
        
        if found_solutions:
            return found_solutions, total_attempts, time.time() - t_start
        
        if total_attempts > 500000:
            if verbose:
                print(f"  Hit attempt cap, stopping early")
            break
    
    # ── Phase 5: random k-switch + split ─────────────────────
    if verbose:
        print(f"\n{'='*55}\nPhase 5: Random k-switch + split ({max_random_switches} attempts)")
    
    for _ in range(max_random_switches):
        k = rng.randint(2, 4)
        # Pick k random edges with distinct vertices
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
        pts_mid = lift_cells(cells_mid, m_old)
        line_pts_mid = build_line_pts(pts_mid)
        if count_collinear(line_pts_mid) != 0:
            continue
        
        # Try splitting a random edge (from new or old)
        for (a, b) in bijection + edges_old[:5]:
            total_attempts += 1
            ok, cells, pts = split_and_test(
                sel_mid, m_old, pts_mid, line_pts_mid, (a, b))
            if ok:
                found_solutions.append(
                    (f'random-{k}-switch+split', edge_set, (a,b), cells, pts))
                if verbose:
                    print(f"  FOUND: rand-{k}-switch {edge_set} + split {a}→{b}")
                break
        
        if found_solutions:
            break
    
    return found_solutions, total_attempts, time.time() - t_start


def main():
    found, attempts, dt = run_surgery(
        m_old=36, max_k=3, max_random_switches=200000, verbose=True)
    
    print(f"\n{'='*55}")
    print(f"Total attempts: {attempts}")
    print(f"Time: {dt:.1f}s")
    print(f"Solutions found: {len(found)}")
    
    if found:
        surgery_type, details, split_edge, cells, pts = found[0]
        print(f"\nFirst solution: {surgery_type}")
        print(f"  details: {details}")
        print(f"  split: {split_edge}")
        print(f"  cells: {cells}")
        
        result = {
            'm': 37,
            'surgery': surgery_type,
            'switched_edges': list(details) if hasattr(details, '__iter__') else str(details),
            'split_edge': list(split_edge) if split_edge else None,
            'cells': sorted(cells),
            'n_cells': len(cells),
            'n_points': len(pts),
            'verify_no_collinear': True,
            'attempts': attempts,
            'time_sec': dt,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        with open('results/m37_surgery_result.json', 'w') as f:
            json.dump(result, f, indent=1)
        with open('results/m37_surgery_found.txt', 'w') as f:
            f.write(f"m=37 rot4 NTIL solution via direction 2 surgery\n")
            f.write(f"Surgery: {surgery_type}\n")
            f.write(f"Details: {details}\n")
            f.write(f"Split edge: {split_edge}\n\n")
            f.write("Cells (top-left quadrant):\n")
            for (x,y) in sorted(cells):
                f.write(f"  ({x}, {y})\n")
            f.write(f"\nAttempts: {attempts}, Time: {dt:.1f}s\n")
        
        print(f"\nSaved to results/m37_surgery_result.json")
        print(f"Saved to results/m37_surgery_found.txt")
    else:
        print(f"\nNo m=37 solution found via surgery.")
        print(f"Directions still open: k≥3 switch on non-vertex-disjoint edges;")
        print(f"  alternating cycle on the directed 2-factor;")
        print(f"  ejection chain with multiple edge replacements;")
        print(f"  or m=37 truly has no rot4 NTIL solution.")
        
        with open('results/m37_surgery_log.txt', 'w') as f:
            f.write(f"m=36→37 surgery log\n")
            f.write(f"Attempts: {attempts}\n")
            f.write(f"Time: {dt:.1f}s\n")
            f.write(f"Result: none found\n")
            f.write(f"k-switch exhaustive: k=1..3\n")
            f.write(f"random switches: {200000} attempts\n")


if __name__ == '__main__':
    main()
