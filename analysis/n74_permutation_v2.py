"""
n=74 C₄: permutation search with 2-cycle blocking
Fix: sigma has no 2-cycles → all M orbits are distinct.
"""

import z3, time

N = 74; M = 37

orbit_pts = {}
for i in range(M):
    for j in range(i, M):
        orbit_pts[(i,j)] = [(i,j),(N-1-j,i),(N-1-i,N-1-j),(j,N-1-i)]

sigma = [z3.Int(f's_{k}') for k in range(M)]
solver = z3.Solver()
solver.set('timeout', 300000)

for k in range(M):
    solver.add(sigma[k] >= 0)
    solver.add(sigma[k] < M)
solver.add(z3.Distinct(sigma))
solver.add(sigma[0] < M//2)  # symmetry breaking

# Block 2-cycles: NOT(sigma[a]==b AND sigma[b]==a) for all a<b
print('Adding 2-cycle constraints...')
for a in range(M):
    for b in range(a+1, M):
        solver.add(z3.Not(z3.And(sigma[a] == b, sigma[b] == a)))
print(f'  Added {M*(M-1)//2} constraints')

print(f'\nStarting search...')
for iteration in range(10000):
    t0 = time.time()
    result = solver.check()
    dt = time.time() - t0
    if result != z3.sat:
        print(f'Iter {iteration}: {result} [{dt:.1f}s]'); break
    
    model = solver.model()
    sigma_vals = [model[sigma[k]].as_long() for k in range(M)]
    
    # DEDUPLICATE orbits
    orbits_set = set()
    for k in range(M):
        orbits_set.add((min(k, sigma_vals[k]), max(k, sigma_vals[k])))
    orbits = list(orbits_set)
    
    pts_set = set()
    for (i,j) in orbits:
        for pt in orbit_pts[(i,j)]:
            pts_set.add(pt)
    
    if len(pts_set) != 2*N or len(orbits) != M:
        # Block this permutation (has duplicates still)
        solver.add(z3.Not(z3.And([sigma[k] == sigma_vals[k] for k in range(M)])))
        if iteration % 500 == 0:
            print(f'  iter {iteration}: dup (pts={len(pts_set)},orb={len(orbits)}) [{dt:.1f}s]', flush=True)
        continue
    
    pts_list = list(pts_set)
    coll_found = False
    coll_triple_orbits = None
    for a in range(len(pts_list)):
        for b in range(a+1, len(pts_list)):
            x1,y1=pts_list[a]; x2,y2=pts_list[b]
            for c in range(b+1, len(pts_list)):
                x3,y3=pts_list[c]
                if (x2-x1)*(y3-y1) == (x3-x1)*(y2-y1):
                    coll_found = True
                    # Find orbits
                    orb_keys = []
                    for pt in [pts_list[a],pts_list[b],pts_list[c]]:
                        for (oi,oj) in orbits:
                            if pt in orbit_pts[(oi,oj)]:
                                orb_keys.append((oi,oj)); break
                    coll_triple_orbits = orb_keys[:3]
                    break
            if coll_found: break
        if coll_found: break
    
    if not coll_found:
        print(f'\n🎉🎉🎉 n=74 C₄ SOLUTION FOUND at iteration {iteration}!')
        print(f'Permutation: {sigma_vals}')
        print(f'Orbits ({len(orbits)}): {sorted(orbits)}')
        break
    
    # Block this specific collinear triple of orbital values
    if coll_triple_orbits and len(coll_triple_orbits) == 3:
        cond_parts = []
        for (oi, oj) in coll_triple_orbits:
            # Which k produces this orbit? Uniquely determined by sigma
            # orbit (oi,oj) is selected iff sigma[oi]==oj OR sigma[oj]==oi
            cond_parts.append(z3.Or(sigma[oi] == oj, sigma[oj] == oi))
        solver.add(z3.Not(z3.And(cond_parts)))
    
    if iteration % 100 == 0:
        print(f'  iter {iteration}: collinear [{dt:.1f}s, {len(solver.assertions())} constraints]', flush=True)

print(f'\nDone: {iteration+1} iterations')
