#!/usr/bin/env python3
"""E3: Direction/slope multiset analysis for known NTIL solutions.

Key insight: In a valid NTIL, each canonical line key (A,B,C) can have at most
2 lifted points (otherwise 3+ collinear). Therefore at most 1 point-pair per key.
So total distinct keys >= C(4m,2) is necessary.

We compare against the total key space available under C4 symmetry.
"""

import json, math, sys, os, time
from collections import defaultdict

# ── Geometry primitives ──────────────────────────────────────────────

def c4(p, r, N):
    x, y = p
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    if r == 3: return (y, N - 1 - x)

def line_key(p, q):
    """Canonical line key (A,B,C) for the line through p and q."""
    x1, y1 = p; x2, y2 = q
    if x1 == x2 and y1 == y2:
        return (0, 0, 0)  # degenerate
    A, B, C = y2 - y1, x1 - x2, x2 * y1 - x1 * y2
    g = math.gcd(math.gcd(abs(A), abs(B)), abs(C))
    if g:
        A //= g; B //= g; C //= g
    if A < 0 or (A == 0 and B < 0):
        A, B, C = -A, -B, -C
    return (A, B, C)

def c4_orbit_key(p, N):
    """C4-orbit key for a pair of points — normalize all 4 rotations."""
    x, y = p
    # Generate all 4 rotations of the point
    rots = [(x, y), (N-1-y, x), (N-1-x, N-1-y), (y, N-1-x)]
    # Canonical: lexicographically smallest
    return min(rots)

def points_from_cells(cells, m):
    N = 2 * m
    return [c4((x, y), r, N) for (x, y) in cells for r in range(4)]

# ── Load known solutions ─────────────────────────────────────────────

def load_solutions():
    solutions = {}
    
    # m05-m19 from results/solutions/
    sol_dir = "results/solutions"
    if os.path.isdir(sol_dir):
        for fname in sorted(os.listdir(sol_dir)):
            if fname.endswith(".json"):
                m_val = int(fname[1:].split(".")[0])  # "m10.json" -> 10
                with open(os.path.join(sol_dir, fname)) as f:
                    d = json.load(f)
                cells = d.get("cells") or d.get("solution", [])
                if cells and d.get("verify", False):
                    solutions[m_val] = [tuple(c) for c in cells]
                    print(f"  Loaded m={m_val}: {len(cells)} cells, verify=True")
    
    # m=36 from Flammenkamp cache
    try:
        sys.path.insert(0, ".")
        from rot4_loader import decode_line
        line = open("flammenkamp_cache/n72_rot4.few").read().strip()
        pts = decode_line(line, 72)
        cells36 = sorted([(x, y) for (x, y) in pts if x < 36 and y < 36])
        # Build in_m36_sol.json if not exists
        if not os.path.exists("results/solutions/m36.json"):
            with open("results/solutions/m36.json", "w") as f:
                json.dump({"m": 36, "n": 72, "cells": cells36, "verify": True}, f)
        solutions[36] = cells36
        print(f"  Loaded m=36: {len(cells36)} cells (Flammenkamp cache)")
    except Exception as e:
        print(f"  m=36 NOT loaded: {e}")
    
    return solutions

# ── Direction key analysis ───────────────────────────────────────────

def analyze_solution(cells, m):
    """Analyze direction key distribution for a known NTIL solution."""
    pts = points_from_cells(cells, m)
    n_pts = len(pts)
    assert n_pts == 4 * m, f"Expected {4*m} points, got {n_pts}"
    
    # Count all point-pair direction keys
    key_counts = defaultdict(int)       # (A,B,C) -> number of point-pairs
    key_pointsets = defaultdict(set)    # (A,B,C) -> set of point indices on that line
    
    for i in range(n_pts):
        xi, yi = pts[i]
        for j in range(i + 1, n_pts):
            xj, yj = pts[j]
            if xi == xj and yi == yj:
                continue
            k = line_key((xi, yi), (xj, yj))
            key_counts[k] += 1
            key_pointsets[k].add(i)
            key_pointsets[k].add(j)
    
    # Stats
    total_pairs = len(key_counts)  # distinct direction keys used
    expected_pairs = n_pts * (n_pts - 1) // 2
    
    # Distribution of key_counts (# of point-pairs per key)
    count_dist = defaultdict(int)
    for c in key_counts.values():
        count_dist[c] += 1
    
    # Distribution of o_ell (# of points per line)
    o_dist = defaultdict(int)
    for s in key_pointsets.values():
        o_dist[len(s)] += 1
    
    # Maximum occupancy
    max_pairs_per_key = max(key_counts.values()) if key_counts else 0
    max_pts_per_line = max(len(s) for s in key_pointsets.values()) if key_pointsets else 0
    
    # Lines with exactly 2 points ("critical lines")
    critical_lines = sum(1 for s in key_pointsets.values() if len(s) == 2)
    
    # C4-orbit analysis: group direction keys by their C4 orbit
    # For each key (A,B,C), the C4 rotation transforms it to another key
    # We track how many distinct C4 orbits the direction keys lie in
    
    return {
        "m": m,
        "n_pts": n_pts,
        "total_pairs_expected": expected_pairs,
        "distinct_keys_used": total_pairs,
        "pairs_ratio": total_pairs / expected_pairs if expected_pairs else 0,
        "max_pairs_per_key": max_pairs_per_key,
        "max_pts_per_line": max_pts_per_line,
        "critical_lines_2pts": critical_lines,
        "lines_with_1pt": o_dist.get(1, 0),
        "lines_with_2pt": o_dist.get(2, 0),
        "lines_with_3pt": o_dist.get(3, 0),
        "lines_with_4pt": o_dist.get(4, 0),
        "lines_with_ge3": sum(v for k, v in o_dist.items() if k >= 3),
        "key_count_dist": dict(count_dist),
        "o_dist": {str(k): v for k, v in sorted(o_dist.items())},
    }

# ── Main ─────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("E3: Slope/Direction Multiset Analysis")
    print("=" * 60)
    
    t0 = time.time()
    solutions = load_solutions()
    print(f"\nLoaded {len(solutions)} solutions in {time.time()-t0:.1f}s")
    print(f"m values: {sorted(solutions.keys())}")
    
    results = {}
    for m in sorted(solutions.keys()):
        cells = solutions[m]
        ti = time.time()
        stats = analyze_solution(cells, m)
        dt = time.time() - ti
        results[m] = stats
        print(f"\n--- m={m:2d} ({dt:.1f}s) ---")
        print(f"  4m={stats['n_pts']} pts, expected pairs={stats['total_pairs_expected']}")
        print(f"  Distinct direction keys used: {stats['distinct_keys_used']}")
        print(f"  Pairs/key ratio (1=ideal): {stats['pairs_ratio']:.4f}")
        print(f"  Max pairs per key: {stats['max_pairs_per_key']}")
        print(f"  Max pts per line:  {stats['max_pts_per_line']}")
        print(f"  Critical lines (o=2): {stats['critical_lines_2pts']}")
        print(f"  Lines with o≥3: {stats['lines_with_ge3']}")
        print(f"  Lines with o=2 count: {stats['lines_with_2pt']}")
        dist_str = "; ".join(f"o={k}:{v}" for k, v in sorted(stats['o_dist'].items()))
        print(f"  o_ell distribution: {{{dist_str}}}")
    
    # Compute scaling trends
    print("\n" + "=" * 60)
    print("SCALING TRENDS")
    print("=" * 60)
    header = f"{'m':>4} {'4m':>4} {'pairs':>10} {'keys':>10} {'ratio':>8} {'max_ppk':>8} {'max_pts':>8} {'crit_ln':>8} {'crit/m':>8}"
    print(header)
    for m in sorted(results.keys()):
        r = results[m]
        print(f"{m:4d} {r['n_pts']:4d} {r['total_pairs_expected']:10d} {r['distinct_keys_used']:10d} "
              f"{r['pairs_ratio']:8.4f} {r['max_pairs_per_key']:8d} {r['max_pts_per_line']:8d} "
              f"{r['critical_lines_2pts']:8d} {r['critical_lines_2pts']/m:8.2f}")
    
    # Save to file
    outfile = "results/e3_slope_stats.json"
    with open(outfile + ".tmp", "w") as f:
        json.dump({"n_ms": len(results), "ms": sorted(results.keys()),
                    "data": {str(k): v for k, v in results.items()}}, f, indent=1)
    os.replace(outfile + ".tmp", outfile)
    print(f"\nSaved to {outfile} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
