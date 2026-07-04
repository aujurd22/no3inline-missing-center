"""
Deep structural analysis of odd-n rot2 solutions from Flammenkamp.
Key questions:
1. How do solution structures differ between odd n?
2. What ring usage patterns distinguish missing-center from has-center?
3. Can we predict missing-center counts?
4. Is there a recurrence relationship?
"""
import urllib.request
from collections import Counter, defaultdict
import math

ALPHABET = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
char_to_val = {c: i for i, c in enumerate(ALPHABET)}

def decode_solution(line):
    rest = line.strip()[1:]
    n = len(rest) // 2
    points = []
    for row in range(n):
        c1 = char_to_val[rest[2*row]]
        c2 = char_to_val[rest[2*row+1]]
        points.append((row, c1))
        points.append((row, c2))
    return n, points

def compute_rings(n, points):
    ctr = (n-1)/2.0
    rings = Counter()
    for i,j in points:
        d2 = int(round((i-ctr)**2 + (j-ctr)**2))
        rings[d2] += 1
    return rings

def has_center(n, points):
    rings = compute_rings(n, points)
    return max(rings.values()) >= 3

def analyze_solutions(n, symm_class='rot2'):
    url = f'https://wwwhomes.uni-bielefeld.de/achim/no3in/download/configurations/n{n}_{symm_class}'
    try:
        with urllib.request.urlopen(url, timeout=30) as f:
            lines = f.read().decode().strip().split('\n')
    except:
        return None
    
    if not lines or not lines[0].strip():
        return None
    
    total = len(lines)
    missing_lines = []
    center_lines = []
    
    for line in lines:
        if not line.strip():
            continue
        n_actual, pts = decode_solution(line)
        if has_center(n_actual, pts):
            center_lines.append(line)
        else:
            missing_lines.append(line)
    
    missing_count = len(missing_lines)
    
    # Analyze ring usage patterns for missing-center solutions
    missing_ring_usage = Counter()
    center_ring_usage = Counter()
    all_ring = Counter()
    
    for line in missing_lines[:min(50, len(missing_lines))]:
        n_act, pts = decode_solution(line)
        rings = compute_rings(n_act, pts)
        for d2, cnt in rings.items():
            missing_ring_usage[d2] += 1
            all_ring[d2] += 1
    
    for line in center_lines[:min(50, len(center_lines))]:
        n_act, pts = decode_solution(line)
        rings = compute_rings(n_act, pts)
        for d2, cnt in rings.items():
            center_ring_usage[d2] += 1
            all_ring[d2] += 1
    
    # Top rings for missing vs center
    missing_top = missing_ring_usage.most_common(10)
    center_top = center_ring_usage.most_common(10)
    
    # How many unique rings are used?
    missing_unique = len(missing_ring_usage)
    center_unique = len(center_ring_usage)
    all_unique = len(all_ring)
    
    return {
        'n': n,
        'total': total,
        'missing': missing_count,
        'rate': missing_count/total*100 if total > 0 else 0,
        'missing_unique_rings': missing_unique,
        'center_unique_rings': center_unique,
        'all_unique_rings': all_unique,
        'missing_top_rings': missing_top,
        'center_top_rings': center_top,
    }

# ================================================================
print("=" * 100)
print("PART 1: Odd n Rot2 Solution Structure Analysis")
print("=" * 100)

results = {}
for n in [11, 13, 15, 17, 19, 21, 27]:
    r = analyze_solutions(n)
    if r:
        results[n] = r
        print(f"\nn={n}: {r['total']} rot2 solutions, {r['missing']} missing-center ({r['rate']:.1f}%)")
        print(f"  Missing solutions use {r['missing_unique_rings']} distinct ring types")
        print(f"  Center  solutions use {r['center_unique_rings']} distinct ring types")
        print(f"  Top missing rings: {r['missing_top_rings'][:5]}")
        print(f"  Top center rings:  {r['center_top_rings'][:5]}")

# ================================================================
print("\n" + "=" * 100)
print("PART 2: Ring Usage Efficiency")
print("=" * 100)

for n in sorted(results.keys()):
    r = results[n]
    # How many rings are there in the full grid?
    m = (n-1)//2
    all_d2 = set()
    for i in range(n):
        for j in range(n):
            a, b = i-m, j-m
            all_d2.add(a*a + b*b)
    total_rings = len(all_d2)
    
    print(f"\nn={n}: {total_rings} total rings, need {n} rings (2/ring)")
    print(f"  Missing uses {r['missing_unique_rings']}/{total_rings} ({r['missing_unique_rings']/total_rings*100:.1f}%)")
    print(f"  Center  uses {r['center_unique_rings']}/{total_rings} ({r['center_unique_rings']/total_rings*100:.1f}%)")

# ================================================================
print("\n" + "=" * 100)
print("PART 3: Ring Collinearity Load — the Hidden Parameter")
print("=" * 100)
print("""
For each distance ring d^2, we compute its "collinearity load" — the number 
of OTHER rings that share collinear triples with it. This is the true 
computational cost measure.
""")

for n in sorted(results.keys()):
    m = (n-1)//2
    
    # Build ring membership: for each d^2, list all grid points at that distance
    ring_points = defaultdict(list)
    for i in range(n):
        for j in range(n):
            a, b = i-m, j-m
            d2 = a*a + b*b
            ring_points[d2].append((i,j))
    
    # For each ring, count how many OTHER rings it interacts with via collinearity
    # We sample rather than full enumeration
    ring_collision_count = Counter()
    rings_list = list(ring_points.keys())
    
    # Sample-based: for each ring, pick a few points and check lines
    for d2, pts in ring_points.items():
        if len(pts) < 2:
            continue
        # Take first 2 points, find all collinear points in other rings
        p1, p2 = pts[0], pts[1]
        dx, dy = p2[0]-p1[0], p2[1]-p1[1]
        
        affected_rings = set()
        for i in range(n):
            # Check if point on line exists in other rings
            if dx != 0 and (i-p1[0]) % dx == 0:
                step = (i-p1[0]) // dx
                j = p1[1] + step * dy
                if 0 <= j < n:
                    a, b = i-m, j-m
                    d2_other = a*a + b*b
                    if d2_other != d2:
                        affected_rings.add(d2_other)
        
        ring_collision_count[d2] = len(affected_rings)
    
    # Distribution
    collision_dist = Counter()
    for d2, cnt in ring_collision_count.items():
        collision_dist[cnt] += 1
    
    min_c, max_c = min(ring_collision_count.values()), max(ring_collision_count.values())
    avg_c = sum(ring_collision_count.values()) / len(ring_collision_count)
    
    print(f"n={n}: {len(ring_collision_count)} rings")
    print(f"  Collinearity load: min={min_c}, max={max_c}, avg={avg_c:.1f}")
    
    # How many rings have load > n? (heavy rings)
    heavy = sum(1 for v in ring_collision_count.values() if v > n)
    print(f"  Heavy rings (load > n={n}): {heavy}")
