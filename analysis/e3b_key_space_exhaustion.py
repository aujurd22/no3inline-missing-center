#!/usr/bin/env python3
"""E3b: Direction key space exhaustion analysis.

In a valid NTIL: each of C(4m,2) point-pairs gets a UNIQUE direction key (A,B,C).
If the TOTAL number of possible direction keys in a 2m×2m grid is less than C(4m,2),
then by pigeonhole principle m has no NTIL solution.

Compute: total_distinct_keys(N) for N=2m grid,
         required_keys = C(2N, 2) = C(4m, 2),
         ratio = required / total.

If ratio > 1 for m=37 → proof of non-existence.
"""

import math, sys, time
from collections import defaultdict

def enumerate_all_keys(N):
    """
    Enumerate ALL distinct line keys (A,B,C) from ALL pairs of points in N×N grid.
    N=2m. Total pairs = C(N², 2) ≈ N⁴/2 = O(m⁴).
    For m=10 (N=20): C(400,2) = 79,800 pairs — feasible.
    For m=37 (N=74): C(5476,2) = 14,991,750 pairs — heavy but doable with optimization.
    """
    keys = set()
    # Strategy: iterate over all direction vectors (dx,dy) and count intercepts
    # For a fixed direction (dx,dy), points (x,y) give intercept C = x*dy - y*dx
    # But we need normalized key (A,B,C) where A,B,C have gcd=1.
    # A = dy, B = -dx, so (A,B) normalized by gcd(dy, dx).
    # Then for each such primitive (A,B), count distinct C values that appear.
    
    # Actually simpler approach: iterate over all point pairs for small-m computation,
    # use symmetric key enumeration for the count.
    # For larger m, use the direction enumeration method.
    
    all_pts = [(x, y) for x in range(N) for y in range(N)]
    n_pts = len(all_pts)
    print(f"  N={N}, cells={n_pts}, pairs={n_pts*(n_pts-1)//2}", flush=True)
    
    for i in range(n_pts):
        x1, y1 = all_pts[i]
        for j in range(i + 1, n_pts):
            x2, y2 = all_pts[j]
            A = y2 - y1
            B = x1 - x2
            C = x2*y1 - x1*y2
            g = math.gcd(math.gcd(abs(A), abs(B)), abs(C))
            if g:
                A //= g; B //= g; C //= g
            if A < 0 or (A == 0 and B < 0):
                A, B, C = -A, -B, -C
            keys.add((A, B, C))
        if (i+1) % (n_pts // 10) == 0:
            print(f"    {i+1}/{n_pts} points, {len(keys)} keys so far...", flush=True)
    
    return keys


def count_keys_by_direction(N):
    """
    More efficient: count keys by primitive direction (A,B).
    For each primitive direction, C varies.
    This avoids iterating over all O(N⁴) pairs.
    """
    # For each primitive direction (A,B) with gcd(A,B)=1:
    # A = dy, B = -dx in normalized form
    # intercept C = x2*y1 - x1*y2 = x1*(y1+A) - y1*(x1-B) ??? 
    # Let me think... C = x2*y1 - x1*y2 where (x2,y2) = (x1-B, y1+A) after normalization.
    # Actually C = x1*(y1+A) - (x1-B)*y1 = x1*y1 + x1*A - x1*y1 + B*y1 = A*x1 + B*y1
    # Wait, check: C = x2*y1 - x1*y2. And (x2,y2) = (x1 - B, y1 + A) approximately...
    # No, (x2,y2) = (x1 - B, y1 + A) would mean dx = x2-x1 = -B, dy = y2-y1 = A.
    # C = x2*y1 - x1*y2 = (x1-B)*y1 - x1*(y1+A) = x1*y1 - B*y1 - x1*y1 - x1*A = -(B*y1 + x1*A)
    # = -(A*x1 + B*y1)
    # Since C is sign-normalized, we can track C' = A*x + B*y for each point (x,y).
    
    # So for a fixed primitive direction (A,B), the intercepts C = A*x + B*y (up to sign)
    # range over all distinct values of A*x + B*y for (x,y) in the N×N grid.
    
    # Number of distinct C values for fixed (A,B) = |{A*x + B*y : 0 ≤ x,y < N}|
    
    total = 0
    primitive_dirs = []
    
    # All directions (A,B) with |A|,|B| < N, gcd(A,B)=1, normalized
    for A in range(-N+1, N):
        for B in range(-N+1, N):
            if A == 0 and B == 0:
                continue
            if math.gcd(abs(A), abs(B)) != 1:
                continue
            # Normalize (A,B) like line_key does
            if A < 0 or (A == 0 and B < 0):
                continue  # skip the "negative" version, we'll only take the canonical one
            primitive_dirs.append((A, B))
    
    print(f"  Primitive directions (A,B) in [-{N-1},{N-1}]²: {len(primitive_dirs)}", flush=True)
    
    for A, B in primitive_dirs:
        intercepts = set()
        for x in range(N):
            for y in range(N):
                intercepts.add(A*x + B*y)
        # Count distinct keys = number of distinct C values for this (A,B)
        # But we need to account for normalization of C by gcd(|A|,|B|,|C|)
        # Since gcd(A,B)=1, gcd(A,B,C) = gcd(1,|C|) = 1 always.
        # So no further normalization needed.
        # (This is the key insight: primitive directions don't need GCD reduction on C)
        
        # Wait, is this true? If gcd(A,B)=1, then gcd(A,B,C) = gcd(1,C) = 1 or C mod gcd(1).
        # Since gcd(1, C) = 1 for any C, yes! No further normalization.
        
        total += len(intercepts)
    
    return total, len(primitive_dirs)


def main():
    print("=" * 60)
    print("E3b: Direction Key Space Exhaustion Analysis")
    print("=" * 60)
    
    t0 = time.time()
    
    print("\n--- Method 1: Brute force enumeration of all keys (small N) ---")
    for m in [5, 6, 7, 8, 9, 10, 11, 12]:
        N = 2 * m
        ti = time.time()
        keys = enumerate_all_keys(N)
        dt = time.time() - ti
        total_keys = len(keys)
        required = (4 * m) * (4 * m - 1) // 2  # C(4m, 2)
        ratio = required / total_keys if total_keys > 0 else 0
        print(f"  m={m:2d} N={N:2d}: total_keys={total_keys:8d}, required={required:8d}, ratio={ratio:.4f} [{dt:.1f}s]")
    
    print("\n--- Method 2: Primitive direction enumeration (fast, large N) ---")
    for m in [5, 10, 15, 20, 25, 30, 36, 37, 40]:
        N = 2 * m
        ti = time.time()
        total_keys, n_dir = count_keys_by_direction(N)
        dt = time.time() - ti
        required = (4 * m) * (4 * m - 1) // 2
        ratio = required / total_keys if total_keys > 0 else 0
        # Also compute max possible keys (if each pair had unique key)
        max_possible = N * N * (N * N - 1) // 2
        print(f"  m={m:2d} N={N:2d}: total_keys={total_keys:10d}, required={required:10d}, "
              f"ratio={ratio:.4f}, max_pair_keys={max_possible} [{dt:.1f}s]")
    
    print(f"\nTotal: {time.time()-t0:.1f}s")
    
    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    print("If ratio < 1 for all m: direction key space is NOT the bottleneck.")
    print("If ratio → 1 as m grows: pigeonhole principle may eventually force o_ℓ ≥ 3.")
    print()
    print("Note: ratio=1 would mean ALL possible direction keys are used,")
    print("leaving zero slack. ratio < 0.5 means half the key space is free.")


if __name__ == "__main__":
    main()
