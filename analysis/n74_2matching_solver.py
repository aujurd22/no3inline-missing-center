"""
2-matching incremental solver for n=74 C₄ solution.
Same approach as before but with cleaner code and more iterations.
"""

import z3, time
from collections import Counter

N = 74
M = 37

# Precompute orbit points
print('Precomputing orbit points...')
orbit_pts = {}
for i in range(M):
    for j in range(M):
        pts = [(i,j), (N-1-j,i), (N-1-i,N-1-j), (j,N-1-i)]
        orbit_pts[(i,j)] = pts

# Z3 2-matching variables
x = [[z3.Bool(f'x_{i}_{j}') for j in range(M)] for i in range(M)]
solver = z3.Solver()
solver.set('timeout', 600000)

# Exactly M orbits
solver.add(z3.Sum([z3.If(x[i][j], 1, 0) for i in range(M) for j in range(M)]) == M)

# Degree constraint: each vertex appears in exactly 2 orbits
for i in range(M):
    deg = z3.Sum([z3.If(x[i][j], 1, 0) for j in range(M)]) + \
          z3.Sum([z3.If(x[j][i], 1, 0) for j in range(M)])
    solver.add(deg == 2)

# Row constraint: each row gets exactly 2 points  
for r in range(N):
    solver.add(z3.Sum([
        z3.If(x[i][j], 1, 0)
        for i in range(M) for j in range(M)
        for ri in [i, j, N-1-i, N-1-j]
        if ri == r
    ]) == 2)

print(f'Starting incremental search...')
print(f'Variables: {M*M} boolean, constraints: degree={M}, row={N}')
print()

for iteration in range(10000):
    t0 = time.time()
    result = solver.check()
    dt = time.time() - t0
    
    if result != z3.sat:
        print(f'Iter {iteration}: {result} [{dt:.1f}s]')
        break
    
    model = solver.model()
    selected = [(i,j) for i in range(M) for j in range(M) if model[x[i][j]] == True]
    
    # Build points
    pts_set = set()
    for (i,j) in selected:
        for pt in orbit_pts[(i,j)]:
            pts_set.add(pt)
    
    if len(pts_set) != 2*N:
        print(f'  iter {iteration}: wrong point count {len(pts_set)}')
        # Add blocking constraint for duplicate-inducing orbits
        duplicate_inducing = []
        pt_counts = Counter()
        for (i,j) in selected:
            for pt in orbit_pts[(i,j)]:
                pt_counts[pt] += 1
        for pt, cnt in pt_counts.items():
            if cnt > 1:
                # Find which orbits produce this point multiple times
                producers = [(i,j) for (i,j) in selected if pt in orbit_pts[(i,j)]]
                if len(producers) >= 2:
                    duplicate_inducing.append(producers[:2])
        for (i1,j1),(i2,j2) in duplicate_inducing:
            solver.add(z3.Not(z3.And(x[i1][j1], x[i2][j2])))
        continue
    
    pts_list = list(pts_set)
    
    # Find collinear triples
    violations = []
    for a in range(len(pts_list)):
        for b in range(a+1, len(pts_list)):
            x1,y1 = pts_list[a]; x2,y2 = pts_list[b]
            for c in range(b+1, len(pts_list)):
                x3,y3 = pts_list[c]
                if (x2-x1)*(y3-y1) == (x3-x1)*(y2-y1):
                    # Find which orbits contain these 3 points
                    orb_keys = []
                    for pt in [pts_list[a], pts_list[b], pts_list[c]]:
                        for (i,j) in selected:
                            if pt in orbit_pts[(i,j)]:
                                orb_keys.append((i,j))
                                break
                    violations.append(tuple(orb_keys[:3]))
                    break
            if violations and len(violations) > 0 and violations[-1] is not None:
                if len(violations[-1]) == 3:
                    pass
        if len(violations) > 200:  # cap for speed
            break
    
    if not violations:
        print(f'\n🎉🎉🎉 n=74 C₄ SOLUTION FOUND at iteration {iteration}!')
        print(f'Orbits ({len(selected)}): {sorted(selected)}')
        break
    
    # Block violations (up to 200 to keep solving fast)
    for orb1, orb2, orb3 in violations[:200]:
        solver.add(z3.Not(z3.And(x[orb1[0]][orb1[1]], x[orb2[0]][orb2[1]], x[orb3[0]][orb3[1]])))
    
    if iteration % 10 == 0:
        constraint_count = len(solver.assertions())
        print(f'  iter {iteration}: {len(violations)} collinear blocked [{dt:.1f}s, {constraint_count} constraints]', flush=True)

print(f'\nDone. Iterations: {iteration+1}')
