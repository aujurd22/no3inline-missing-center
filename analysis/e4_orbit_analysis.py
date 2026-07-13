#!/usr/bin/env python3
"""E4: C4-orbit analysis of direction key space.

For N×N grid (N=2m), count the number of C4-orbits of distinct direction keys.
Compare with requirement: C(4m,2) distinct keys needed for NTIL.
Since each C4-orbit contributes at most 4 keys, need ≥ C(4m,2)/4 orbits.

Also check: does this ratio approach 1 as m grows?
Does m odd/even create different orbit counts?
"""

import math, sys, time
from collections import defaultdict


def line_key(p, q):
    x1, y1 = p; x2, y2 = q
    if x1 == x2 and y1 == y2:
        return (0, 0, 0)
    A = y2 - y1
    B = x1 - x2
    C = x2 * y1 - x1 * y2
    g = math.gcd(math.gcd(abs(A), abs(B)), abs(C))
    if g:
        A //= g; B //= g; C //= g
    if A < 0 or (A == 0 and B < 0):
        A, B, C = -A, -B, -C
    return (A, B, C)


def c4_key_transform(k, N):
    """Apply C4 rotation to a direction key (A,B,C).
    Returns the new key after normalizing.
    """
    A, B, C = k
    # Under 90° rotation: (A,B,C) → (-B, A, C - B*(N-1))
    A2, B2, C2 = -B, A, C - B * (N - 1)
    # Normalize
    g = math.gcd(math.gcd(abs(A2), abs(B2)), abs(C2))
    if g:
        A2 //= g; B2 //= g; C2 //= g
    if A2 < 0 or (A2 == 0 and B2 < 0):
        A2, B2, C2 = -A2, -B2, -C2
    return (A2, B2, C2)


def orbit_of_key(k, N):
    """Return the set of all keys in the C4-orbit of k."""
    orbit = set()
    kk = k
    for _ in range(4):
        orbit.add(kk)
        kk = c4_key_transform(kk, N)
    return frozenset(orbit)


def enumerate_all_keys_and_orbits(N):
    """Enumerate all distinct keys and their C4-orbits in N×N grid."""
    all_pts = [(x, y) for x in range(N) for y in range(N)]
    n_pts = len(all_pts)
    print(f"  N={N}, cells={n_pts}, pairs={n_pts*(n_pts-1)//2}", flush=True)
    
    keys = set()
    orbit_map = {}  # key → orbit_id
    
    for i in range(n_pts):
        x1, y1 = all_pts[i]
        for j in range(i + 1, n_pts):
            x2, y2 = all_pts[j]
            k = line_key((x1, y1), (x2, y2))
            keys.add(k)
        if (i + 1) % max(1, n_pts // 20) == 0:
            print(f"    {i+1}/{n_pts} pts, {len(keys)} keys", flush=True)
    
    # Build orbits
    orbits = {}
    for k in keys:
        orb = orbit_of_key(k, N)
        if orb not in orbits:
            orbits[orb] = len(orb)
    
    # Orbit size distribution
    size_dist = defaultdict(int)
    for sz in orbits.values():
        size_dist[sz] += 1
    
    orbit_list = list(orbits.keys())
    
    return len(keys), len(orbit_list), size_dist, orbit_list[:5]


# More efficient: count by method 1 (brute force) for small m, 
# use formula for large m

def main():
    print("=" * 70)
    print("E4: C4-Orbit Direction Key Space Analysis")
    print("=" * 70)
    
    results = {}
    t0 = time.time()
    
    for m in [5, 6, 7, 8, 9, 10, 11, 12, 13, 14]:
        N = 2 * m
        ti = time.time()
        n_keys, n_orbits, sz_dist, _ = enumerate_all_keys_and_orbits(N)
        dt = time.time() - ti
        
        required_pairs = (4 * m) * (4 * m - 1) // 2
        min_orbits_needed = (required_pairs + 3) // 4  # ceil division by 4
        orbit_ratio = min_orbits_needed / n_orbits if n_orbits else 0
        
        results[m] = {
            'N': N, 'n_keys': n_keys, 'n_orbits': n_orbits,
            'required': required_pairs, 'min_orbits_needed': min_orbits_needed,
            'orbit_ratio': orbit_ratio,
            'size_dist': {str(k): v for k, v in sorted(sz_dist.items())}
        }
        
        print(f"\n--- m={m:2d} N={N:2d} ({dt:.1f}s) ---")
        print(f"  Distinct keys: {n_keys}")
        print(f"  C4-orbits: {n_orbits}")
        print(f"  Required pairs for NTIL: {required_pairs}")
        print(f"  Min orbits needed (÷4): {min_orbits_needed}")
        print(f"  Orbit ratio (need/avail): {orbit_ratio:.4f}")
        print(f"  Orbit size distribution: {dict(sz_dist)}")
    
    print("\n" + "=" * 70)
    print("TREND TABLE")
    print("=" * 70)
    header = f"{'m':>4} {'N':>4} {'keys':>8} {'orbits':>8} {'needed':>8} {'ratio':>8} {'orbit_sizes':>20}"
    print(header)
    for m in sorted(results.keys()):
        r = results[m]
        sz_str = "/".join(f"{k}={v}" for k,v in sorted({int(k):v for k,v in r['size_dist'].items()}.items()))
        print(f"{m:4d} {r['N']:4d} {r['n_keys']:8d} {r['n_orbits']:8d} "
              f"{r['min_orbits_needed']:8d} {r['orbit_ratio']:8.4f} {sz_str}")
    
    print(f"\nTotal: {time.time()-t0:.1f}s")

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("If orbit_ratio → 1: orbit space is the bottleneck.")
    print("If orbit_ratio → 0: plenty of orbital slack.")
    print("If orbit_ratio < 1 for m=37: orbit space NOT the bottleneck.")
    print()
    print("But: NTIL points must be C4-symmetric, so their line keys also form")
    print("C4-orbits. The key question: are there enough ORBIT STARTING-POINTS")
    print("in the C4-invariant subset of the 2m×2m grid?")


if __name__ == "__main__":
    main()
