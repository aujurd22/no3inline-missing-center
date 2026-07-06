"""
n=74 C₄ solution search via partition:
  Fix m1 orbits from n=56 DB (vertices 0..27)
  Search m2 orbits for vertices 28..36 using Z3 permutation encoding
  
Key insight: 
  - Row constraint = degree constraint for any 2-regular graph ✓
  - A permutation sigma on m vertices gives each vertex degree 2
    deg(k) = 1 (as i in (k,sigma[k])) + 1 (as j in (sigma^-1[k],k)) = 2
  - So we just search for a permutation on 9 vertices that avoids collinearity
"""

import urllib.request, z3, time, sys
from collections import Counter

char_map = {c:i for i,c in enumerate('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')}

def get_all_unique_patterns(n):
    """Get all unique orbit patterns for n."""
    url = f'https://wwwhomes.uni-bielefeld.de/achim/no3in/download/configurations/n{n}_rot4'
    lines = urllib.request.urlopen(url, timeout=10).read().decode().strip().split(chr(10))
    lines = [l for l in lines if l.strip()]
    all_sols = set()
    for data in lines:
        if data[0] == 'o': data = data[1:]
        pts = {}
        for r in range(n):
            pts[r] = (char_map[data[2*r]], char_map[data[2*r+1]])
        pt_set = {(r,c) for r in range(n) for c in pts[r]}
        remaining = set(pt_set)
        orbits = []
        while remaining:
            p = remaining.pop()
            orb = {p}; np = p
            for _ in range(3):
                np = (n-1-np[1], np[0]); orb.add(np)
            remaining -= orb
            fd = lambda pt: (min(pt[0],n-1-pt[0]), min(pt[1],n-1-pt[1]))
            orbits.append(min(fd(p) for p in orb))
        all_sols.add(tuple(sorted(orbits)))
    return list(all_sols)

def c4_points(orbits, n):
    pts_set = set()
    for (i,j) in orbits:
        pts_set.update([(i,j), (n-1-j,i), (n-1-i,n-1-j), (j,n-1-i)])
    return pts_set

def points_from_orbit(i, j, n):
    return [(i,j), (n-1-j,i), (n-1-i,n-1-j), (j,n-1-i)]

def check_collinear(p1, p2, p3):
    x1,y1=p1; x2,y2=p2; x3,y3=p3
    return (x2-x1)*(y3-y1) == (x3-x1)*(y2-y1)

def solve_partition_z3(fixed_orbits, m2_start, m2_end, N):
    """
    Fixed: orbits covering vertices < m2_start (already have deg=2)
    Search: permutation on vertices [m2_start, m2_end] (m2_count = m2_end-m2_start+1)
    
    Returns list of (i,j) orbits for the variable part if found, else None.
    """
    m2_count = m2_end - m2_start + 1
    
    # Precompute fixed orbit points
    fixed_pts_list = list(c4_points(fixed_orbits, N))
    
    # Z3 permutation variables: s[k] = j value for orbit (k+offset, s[k]+offset)
    offset = m2_start
    sigma = [z3.Int(f's_{k}') for k in range(m2_count)]
    
    solver = z3.Solver()
    solver.set('timeout', 120000)
    
    # Each s[k] ∈ [0, m2_count-1]
    for k in range(m2_count):
        solver.add(sigma[k] >= 0)
        solver.add(sigma[k] < m2_count)
    
    # Permutation: Distinct values
    solver.add(z3.Distinct(sigma))
    
    # Precompute collinear triples involving variable orbits
    # For each choice of 3 variable orbits (or 2 var + 1 fixed), check collinearity
    
    var_orbits_local = [(i, j) for i in range(m2_count) for j in range(i, m2_count)]
    # Actually we'll check on the fly since variable values change
    
    # Precompute: for each pair of variable orbits, which collinear triples exist?
    # Determine if 3 specific VAR orbits (represented as their indices) could be collinear
    
    # For a given set of 3 var orbits {(a,b), (c,d), (e,f)}:
    #   check if any point from each is collinear with points from the other two
    
    # Build ALL possible var orbit point sets
    var_pts = {}
    for i in range(m2_count):
        for j in range(i, m2_count):
            var_pts[(i,j)] = points_from_orbit(i+offset, j+offset, N)
    
    print(f'  Precomputing cross-boundary collinearity constraints...')
    
    # Case 1: All 3 orbits from var set → need permutation constraint
    # For any 3 var orbits (i1,j1), (i2,j2), (i3,j3), if they're collinear, block them
    var_var_collinear = []
    var_keys = list(var_pts.keys())
    for vk1_idx in range(len(var_keys)):
        vk1 = var_keys[vk1_idx]
        for vk2_idx in range(vk1_idx+1, len(var_keys)):
            vk2 = var_keys[vk2_idx]
            for vk3_idx in range(vk2_idx+1, len(var_keys)):
                vk3 = var_keys[vk3_idx]
                coll = False
                for p1 in var_pts[vk1]:
                    for p2 in var_pts[vk2]:
                        for p3 in var_pts[vk3]:
                            if check_collinear(p1, p2, p3):
                                coll = True; break
                        if coll: break
                    if coll: break
                if coll:
                    # Block: NOT (selected or all three)
                    var_var_collinear.append((vk1, vk2, vk3))
    
    print(f'    var×var×var: {len(var_var_collinear)} constraints')
    
    # Case 2: 2 var orbits + 1 fixed orbit  
    fix_var_var_collinear = []
    for fi, fj in fixed_orbits:
        f_pts = points_from_orbit(fi, fj, N)
        for vk1_idx in range(len(var_keys)):
            vk1 = var_keys[vk1_idx]
            for vk2_idx in range(vk1_idx+1, len(var_keys)):
                vk2 = var_keys[vk2_idx]
                coll = False
                for p1 in f_pts:
                    for p2 in var_pts[vk1]:
                        for p3 in var_pts[vk2]:
                            if check_collinear(p1, p2, p3):
                                coll = True; break
                        if coll: break
                    if coll: break
                if coll:
                    fix_var_var_collinear.append((vk1, vk2))
    
    print(f'    fix×var×var: {len(fix_var_var_collinear)} constraints')
    
    # Case 3: 1 var orbit + 2 fixed orbits
    fix_fix_var_collinear = []
    for fk1_idx in range(len(fixed_orbits)):
        fi1, fj1 = fixed_orbits[fk1_idx]
        f_pts1 = points_from_orbit(fi1, fj1, N)
        for fk2_idx in range(fk1_idx+1, len(fixed_orbits)):
            fi2, fj2 = fixed_orbits[fk2_idx]
            f_pts2 = points_from_orbit(fi2, fj2, N)
            for vk in var_keys:
                coll = False
                for p1 in f_pts1:
                    for p2 in f_pts2:
                        for p3 in var_pts[vk]:
                            if check_collinear(p1, p2, p3):
                                coll = True; break
                        if coll: break
                    if coll: break
                if coll:
                    fix_fix_var_collinear.append(vk)
    
    print(f'    fix×fix×var: {len(fix_fix_var_collinear)} constraints')
    
    # Now build Z3 constraints based on permutation
    # For each var-orbit triple that's collinear, block it
    # Block condition: (a, sigma[a]), (b, sigma[b]), (c, sigma[c]) are the three orbits
    # where the orbit is (k+offset, sigma[k]+offset) OR (sigma[k]+offset, k+offset) 
    # depending on whether k < sigma[k]
    
    # Since sigma is a permutation, edges are {(k, sigma[k])}
    # So orbit (k, sigma[k]) where k ranges over 0..m2_count-1
    # But orbits are symmetric, so (k, sigma[k]) and (sigma[k], k) are the same
    
    # To get the orbit key: (min(k, sigma[k]), max(k, sigma[k]))
    def orbit_key(k):
        j = sigma[k].as_long() if hasattr(sigma[k], 'as_long') else sigma[k]
        return (min(k, j), max(k, j))
    
    # We'll add constraints block by block in incremental solving
    # For Z3 incremental: block all possibilities upfront where the 3 selected orbits are collinear
    
    # Convert collinear triples to Z3 constraints
    # The permutation sigma determines which orbits are selected.
    # Orbit (k, sigma[k]) is selected iff (k, sigma[k]).
    # Since we use min/max ordering, orbit (a,b) is selected iff 
    # sigma[a] = b OR sigma[b] = a
    
    print(f'  Adding constraints to Z3...')
    
    # Add all permutation constraints upfront
    # For each var×var×var collinear triple: NOT all 3 selected
    for vk1, vk2, vk3 in var_var_collinear:
        i1,j1 = vk1; i2,j2 = vk2; i3,j3 = vk3
        
        # These orbits are selected if the permutation contains their edges
        # Orbit (i,j) selected iff sigma[i-offset] = j-offset or sigma[j-offset] = i-offset
        
        def orbit_cond(k, j):
            """k is the vertex index (0..m2_count-1), j is the target vertex (0..m2_count-1)."""
            return z3.Or(sigma[k] == j, sigma[j] == k)
        
        cond = z3.Not(z3.And(
            orbit_cond(i1, j1),
            orbit_cond(i2, j2),
            orbit_cond(i3, j3)
        ))
        solver.add(cond)
    
    # For fix×var×var: the var orbits are selected, the fixed orbit is fixed
    for vk1, vk2 in fix_var_var_collinear:
        i1,j1 = vk1; i2,j2 = vk2
        cond = z3.Not(z3.And(
            orbit_cond(i1, j1),
            orbit_cond(i2, j2)
        ))
        solver.add(cond)
    
    # For fix×fix×var: the var orbit is selected, the fixed orbits are fixed
    for vk in fix_fix_var_collinear:
        i,j = vk
        cond = z3.Not(orbit_cond(i, j))
        solver.add(cond)
    
    print(f'  Solving...')
    t0 = time.time()
    result = solver.check()
    dt = time.time() - t0
    
    if result == z3.sat:
        model = solver.model()
        sigma_vals = [model[sigma[k]].as_long() for k in range(m2_count)]
        
        # Build orbits from permutation
        var_orbits = []
        for k in range(m2_count):
            i, j = k + offset, sigma_vals[k] + offset
            var_orbits.append((min(i,j), max(i,j)))
        
        print(f'  ✅ Solution found in {dt:.1f}s!')
        print(f'  Permutation: {sigma_vals}')
        print(f'  Variable orbits: {sorted(var_orbits)}')
        return var_orbits
    else:
        print(f'  Z3 result: {result} [{dt:.1f}s]')
        return None


# ===== Main =====
N = 74
M = 37
print('=' * 60)
print(f'n={N} C₄ Solution via partition search')
print('=' * 60)

# Get n=56 patterns
print('\nFetching n=56 patterns...')
patterns_56 = get_all_unique_patterns(56)
print(f'Unique patterns: {len(patterns_56)}')

# Filter to those using only vertices 0..27
valid_patterns = [p for p in patterns_56 if max(max(i,j) for i,j in p) <= 27]
print(f'Patterns with vertices ⊆ [0,27]: {len(valid_patterns)}')

# Also try smaller n patterns for the fixed part  
# n=54 (m=27) has 7696 solutions, so m=27 fixed + m=10 search = 37
print('\nFetching n=54 patterns...')
patterns_54 = get_all_unique_patterns(54)
valid_54 = [p for p in patterns_54 if max(max(i,j) for i,j in p) <= 26]
print(f'n=54 patterns with vertices ⊆ [0,26]: {len(valid_54)}')

# n=52 (m=26) + m=11 = 37
print('\nFetching n=52 patterns...')
patterns_52 = get_all_unique_patterns(52)
valid_52 = [p for p in patterns_52 if max(max(i,j) for i,j in p) <= 25]
print(f'n=52 patterns with vertices ⊆ [0,25]: {len(valid_52)}')

# Try different partition strategies
partitions = [
    # (fixed_orbit_list, m2_start, m2_end)
]

# Strategy 1: n=56 (m=28) + remaining 9 vertices
print('\n\n=== Strategy 1: 28 + 9 partition ===')
for pat_idx, pat in enumerate(valid_patterns[:10]):
    print(f'\n  Pattern {pat_idx}:')
    result = solve_partition_z3(list(pat), 28, 36, N)
    if result:
        all_orbits = list(pat) + result
        # Verify
        deg = Counter()
        for i,j in all_orbits: deg[i] += 1; deg[j] += 1
        rows = Counter()
        for i,j in all_orbits: rows[i] += 1; rows[j] += 1; rows[N-1-i] += 1; rows[N-1-j] += 1
        pts = list(c4_points(all_orbits, N))
        coll = False
        for a in range(len(pts)):
            for b in range(a+1,len(pts)):
                for c in range(b+1,len(pts)):
                    if check_collinear(pts[a],pts[b],pts[c]): coll=True; break
                if coll: break
            if coll: break
        
        print(f'  Degree: {"✅" if all(d==2 for d in deg.values()) else "❌"}')
        print(f'  Rows: {"✅" if all(r==2 for r in rows.values()) else "❌"}')
        print(f'  Collinear: {"✅ NONE" if not coll else "❌ FOUND"}')
        
        if all(d==2 for d in deg.values()) and all(r==2 for r in rows.values()) and not coll:
            print(f'\n🎉🎉🎉 VALID n={N} C₄ SOLUTION FOUND!')
            print(f'All 37 orbits: {sorted(all_orbits)}')
            sys.exit(0)

# Strategy 2: n=54 (m=27) + remaining 10 vertices
print('\n\n=== Strategy 2: 27 + 10 partition ===')
for pat_idx, pat in enumerate(valid_54[:10]):
    print(f'\n  Pattern {pat_idx}:')
    result = solve_partition_z3(list(pat), 27, 36, N)
    if result:
        all_orbits = list(pat) + result
        deg = Counter()
        for i,j in all_orbits: deg[i] += 1; deg[j] += 1
        rows = Counter()
        for i,j in all_orbits: rows[i] += 1; rows[j] += 1; rows[N-1-i] += 1; rows[N-1-j] += 1
        pts = list(c4_points(all_orbits, N))
        coll = False
        for a in range(len(pts)):
            for b in range(a+1,len(pts)):
                for c in range(b+1,len(pts)):
                    if check_collinear(pts[a],pts[b],pts[c]): coll=True; break
                if coll: break
            if coll: break
        print(f'  Degree: {"✅" if all(d==2 for d in deg.values()) else "❌"}')
        print(f'  Rows: {"✅" if all(r==2 for r in rows.values()) else "❌"}')
        print(f'  Collinear: {"✅ NONE" if not coll else "❌ FOUND"}')
        if all(d==2 for d in deg.values()) and all(r==2 for r in rows.values()) and not coll:
            print(f'\n🎉🎉🎉 VALID n={N} C₄ SOLUTION FOUND!')
            print(f'All 37 orbits: {sorted(all_orbits)}')
            sys.exit(0)

print('\n\nNo solution found with these partition strategies.')
