"""
Pattern discovery from Flammenkamp rot4 solutions.

Strategy:
1. Look at the row-by-row permutation structure of known C₄ solutions
2. Each row r has 2 points at columns c_r1, c_r2
3. C₄ symmetry means the permutation structure has special properties
4. Find a formula that generates valid patterns for ANY n
"""
import urllib.request, sys
from collections import Counter

char_to_val = {c:i for i,c in enumerate('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')}

def parse_solution(line):
    data = line.strip()
    if data[0] == 'o': data = data[1:]
    n = len(data)//2
    pts = {}
    for r in range(n):
        c1, c2 = char_to_val[data[2*r]], char_to_val[data[2*r+1]]
        pts[r] = (c1, c2)
    return pts, n

# Get all rot4 solutions for known n
def get_all_solutions(n):
    url = f'https://wwwhomes.uni-bielefeld.de/achim/no3in/download/configurations/n{n}_rot4'
    try:
        lines = urllib.request.urlopen(url, timeout=10).read().decode().strip().split('\n')
        return [l for l in lines if l.strip() and l[0] == 'o']
    except:
        return []

# Approach: For each solution, extract the "pair permutation" 
# C₄ symmetry means: if row r has columns (a,b), then row n-1-r has columns (n-1-b, n-1-a)
# This is a necessary condition for C₄ symmetry.

print("Analyzing C₄ solution structure:")
print("For each solution, check: does row r's pair determine row n-1-r's pair?")
print()

for n in [12, 14, 16, 18, 20, 22, 24]:
    lines = get_all_solutions(n)
    m = n // 2
    print(f"n={n} (m={m}): {len(lines)} solutions")
    
    # Show first solution's row-by-row structure
    pts, _ = parse_solution(lines[0])
    print(f"  Row pairs (first 4 + last 4 rows):")
    for r in list(range(4)) + list(range(n-4, n)):
        c1, c2 = pts[r]
        print(f"    row {r:2d}: ({c1:2d}, {c2:2d})", end="")
        # Check C₄ condition with row n-1-r
        r_flip = n - 1 - r
        if r_flip > r:
            cf1, cf2 = pts[r_flip]
            expected = (n-1-c2, n-1-c1)
            match = "✓" if (cf1, cf2) == expected else f"✗ expected ({expected[0]},{expected[1]})"
            print(f"  ↔ row {r_flip:2d}: ({cf1:2d}, {cf2:2d}) {match}", end="")
        print()
    print()

# Now, the key insight: C₄ rot4 solutions have a VERY specific structure.
# If row r has columns (a, b), then due to C₄ symmetry:
# The point (r, a) rotates to (n-1-a, r), which must be in the solution
# The point (r, b) rotates to (n-1-b, r), which must be in the solution
# 
# So row n-1-a must contain r, and row n-1-b must contain r.
#
# This means the permutation pairs form a very structured pattern.
# In fact, this is equivalent to finding a fixed-point-free involution
# or a matching in a certain graph.

# Let me extract the C₄ cycle structure
print("=" * 60)
print("C₄ cycle decomposition of solutions:")
print()

for n in [12, 14, 16, 18]:
    lines = get_all_solutions(n)
    pts, _ = parse_solution(lines[0])
    m = n // 2
    
    # Build the C₄ permutation: each grid point (r,c) maps to (n-1-c, r)
    # For a rot4 solution, this permutation partitions the 2n points into 4-cycles
    
    pt_set = {(r, c) for r in range(n) for c in pts[r]}
    remaining = set(pt_set)
    cycles = []
    
    while remaining:
        p = remaining.pop()
        cycle = [p]
        np = p
        for _ in range(3):
            r, c = np
            np = (n-1-c, r)
            if np in pt_set:
                cycle.append(np)
        if len(cycle) == 4:
            cycles.append(cycle)
        remaining -= set(cycle)
    
    # The cycles tell us which rows are linked
    print(f"n={n}: {len(cycles)} C₄ cycles")
    
    # Show first few cycles
    for cycle in cycles[:4]:
        rows = [p[0] for p in cycle]
        cols = [p[1] for p in cycle]
        print(f"  cycle: rows {rows}  cols {cols}")
    
    if len(cycles) > 4:
        print(f"  ... and {len(cycles)-4} more cycles")
    print()

# Now let's try to extract the C₄ orbit selection pattern
# Each cycle of 4 points comes from one C₄ orbit (i,j)
# Orbit (i,j) generates 4 points, which is exactly one cycle
# So the number of cycles = number of orbits = m = n/2

# Let me find the orbit representatives
print("=" * 60)
print("Orbit representatives (canonical form):")
print()

def fundamental_rep(pt, n):
    i, j = pt
    positions = [(i,j), (n-1-j,i), (n-1-i,n-1-j), (j,n-1-i)]
    best = None
    for (ci,cj) in positions:
        fd_i, fd_j = min(ci, n-1-ci), min(cj, n-1-cj)
        cand = (fd_i, fd_j)
        if best is None or cand < best: best = cand
    return best

for n in [12, 14, 16, 18, 24, 30, 36]:
    lines = get_all_solutions(n)
    if not lines:
        print(f"n={n}: no data")
        continue
    m = n // 2
    
    # Show first 2 solutions
    for sol_idx in range(min(2, len(lines))):
        pts, _ = parse_solution(lines[sol_idx])
        pt_set = {(r,c) for r in range(n) for c in pts[r]}
        
        remaining = set(pt_set)
        orbits = []
        while remaining:
            p = remaining.pop()
            orb = {p}
            np = p
            for _ in range(3):
                np = (n-1-np[1], np[0])
                orb.add(np)
            remaining -= orb
            orbits.append(fundamental_rep(p, n))
        
        orbits.sort()
        i_dist = dict(sorted(Counter(r[0] for r in orbits).items()))
        
        # Count: how many i values have 1 orbit vs 2 orbits?
        ones = sum(1 for v in i_dist.values() if v == 1)
        twos = sum(1 for v in i_dist.values() if v == 2)
        
        print(f"n={n} sol#{sol_idx}: {len(orbits)} orbits, 1×{ones} 2×{twos}, reps={orbits}")
    print()
