"""
n=74 C₄ solution using PySAT (CaDiCaL/Glucose) with ALL constraints upfront.
Optimized precomputation using line-hash method.
"""

import math, time
from collections import defaultdict

N = 74; M = 37

print('=== n=74 C₄ SAT Solver (PySAT + CaDiCaL) ===\n')

# === Step 1: Precompute collinear orbit triples (optimized) ===
print('Step 1: Precomputing collinear orbit triples...')

t0 = time.time()

# Generate all orbit points with their orbit labels
pt_to_orbits = defaultdict(set)
orbit_pts_map = {}  # orbit -> its 4 points
for i in range(M):
    for j in range(i, M):
        pts = [(i,j),(N-1-j,i),(N-1-i,N-1-j),(j,N-1-i)]
        orbit_pts_map[(i,j)] = pts
        for pt in pts:
            pt_to_orbits[pt].add((i,j))

all_pts = list(set(pt_to_orbits.keys()))
n_pts = len(all_pts)
print(f'  Unique points: {n_pts}')

# Build pt_idx for fast lookup
pt_to_idx = {pt: idx for idx, pt in enumerate(all_pts)}
idx_to_pt = {idx: pt for pt, idx in pt_to_idx.items()}

# Step 1a: Map each pair of points to its line
# For line key, use (dx, dy, c) where dx,dy reduced, c = x1*dy - y1*dx
line_to_pts = defaultdict(set)  # line key -> set of point indices

for a in range(n_pts):
    x1,y1 = all_pts[a]
    for b in range(a+1, n_pts):
        x2,y2 = all_pts[b]
        dx, dy = x2-x1, y2-y1
        g = math.gcd(abs(dx), abs(dy))
        dx_r, dy_r = dx//g, dy//g
        
        if dx_r < 0 or (dx_r == 0 and dy_r < 0):
            dx_r, dy_r = -dx_r, -dy_r
        
        # Line constant: c = cross product
        c = x1*dy_r - y1*dx_r
        key = (dx_r, dy_r, c)
        line_to_pts[key].add(a)
        line_to_pts[key].add(b)

dt1 = time.time() - t0
print(f'  Line pairs: {sum(len(v)//2 for v in line_to_pts.values())} [{dt1:.1f}s]')

# Step 1b: For each line with 3+ points, extract collinear orbit triples
collinear_triples = set()
lines_with_3plus = sum(1 for v in line_to_pts.values() if len(v) >= 3)
print(f'  Lines with 3+ points: {lines_with_3plus}')

processed = 0
for key, pt_indices in line_to_pts.items():
    if len(pt_indices) < 3:
        continue
    
    pt_indices = sorted(pt_indices)
    n_line = len(pt_indices)
    
    # Collect orbits for each point on this line
    line_orbit_sets = [list(pt_to_orbits[all_pts[idx]]) for idx in pt_indices]
    
    # For each triple of points on this line
    for a in range(n_line):
        orbits_a = line_orbit_sets[a]
        for b in range(a+1, n_line):
            orbits_b = line_orbit_sets[b]
            for c in range(b+1, n_line):
                orbits_c = line_orbit_sets[c]
                # Check orbit combinations
                for oa in orbits_a:
                    for ob in orbits_b:
                        for oc in orbits_c:
                            if oa != ob and ob != oc and oa != oc:
                                key_t = tuple(sorted([oa, ob, oc]))
                                collinear_triples.add(key_t)
    
    processed += 1
    if processed % 50000 == 0:
        print(f'  Processed {processed} lines, found {len(collinear_triples)} triples...', end='\r')

dt2 = time.time() - t0
print(f'  Done: {len(collinear_triples)} collinear triples in {dt2:.1f}s')

# Filter: remove triples where a vertex appears 3+ times (blocked by degree)
filtered = set()
for (i1,j1),(i2,j2),(i3,j3) in collinear_triples:
    vc = defaultdict(int)
    for v in [i1,j1,i2,j2,i3,j3]: vc[v] += 1
    if max(vc.values()) <= 2:
        filtered.add(((i1,j1),(i2,j2),(i3,j3)))

print(f'  After filtering: {len(filtered)} triples')

# === Step 2: SAT encoding ===
print('\nStep 2: Encoding as SAT...')

# Map orbit (i,j) to a SAT variable (1-indexed)
var = {}
for i in range(M):
    for j in range(i, M):
        var[(i,j)] = len(var) + 1

n_vars = len(var)
print(f'  Variables: {n_vars}')

# Build clauses
clauses = []

# Degree constraint: each vertex appears in exactly 2 orbits
def get_incident_vars(i):
    """Get all orbit variables incident to vertex i."""
    result = []
    for j in range(0, i):
        result.append(var[(j,i)])
    for j in range(i, M):
        result.append(var[(i,j)])
    return list(set(result))

print('  Adding degree constraints...')
for i in range(M):
    inc_vars = get_incident_vars(i)
    k = len(inc_vars)
    # At-most-2: no triple all true
    for a in range(k):
        for b in range(a+1, k):
            for c in range(b+1, k):
                clauses.append([-inc_vars[a], -inc_vars[b], -inc_vars[c]])
    # At-least-2: for each set of k-1 vars, at least 1 true
    for a in range(k):
        clause = [inc_vars[x] for x in range(k) if x != a]
        clauses.append(clause)

# Row constraint is AUTOMATICALLY satisfied by degree constraint for C₄ solutions.
# Proof: orbit (i,j) contributes to rows {i, j, N-1-i, N-1-j}.
# Low rows [0,M-1]: row r gets deg(r) = 2 points ✓
# High rows [M,N-1]: row N-1-k gets deg(k) = 2 points ✓
# So we can skip explicit row clauses (saves ~490M clauses).

# Collinear constraints: ternary clauses
print(f'  Adding {len(filtered)} collinear clauses...')
for (i1,j1),(i2,j2),(i3,j3) in filtered:
    v1 = var[(min(i1,j1), max(i1,j1))]
    v2 = var[(min(i2,j2), max(i2,j2))]
    v3 = var[(min(i3,j3), max(i3,j3))]
    clauses.append([-v1, -v2, -v3])

total_clauses = len(clauses)
print(f'  Total clauses: {total_clauses}')

# === Step 3: Solve with CaDiCaL ===
print('\nStep 3: Solving with CaDiCaL...')

try:
    from pysat.solvers import Solver
    
    # Try CaDiCaL first, fall back to Glucose
    solver_name = None
    for s in ['cadical153', 'glucose4', 'glucose42', 'minisat22']:
        try:
            test = Solver(name=s)
            test.delete()
            solver_name = s
            break
        except:
            continue
    
    print(f'  Using: {solver_name}')
    
    s = Solver(name=solver_name)
    
    t0 = time.time()
    for clause in clauses:
        s.add_clause(clause)
    dt_add = time.time() - t0
    print(f'  Added clauses in {dt_add:.1f}s')
    
    t0 = time.time()
    result = s.solve()
    dt_solve = time.time() - t0
    
    print(f'  Result: {result} [{dt_solve:.1f}s]')
    
    if result:
        model = s.get_model()
        selected = []
        for (i,j), v in var.items():
            if model[v-1] > 0:
                selected.append((i,j))
        
        print(f'\n🎯 Solution: {len(selected)} orbits')
        print(f'Orbits: {sorted(selected)}')
        
        # Verify
        pts_set = set()
        for (i,j) in selected:
            pts_set.update([(i,j),(N-1-j,i),(N-1-i,N-1-j),(j,N-1-i)])
        print(f'Points: {len(pts_set)} (expected {2*N})')
        
        coll = False
        pts_list = list(pts_set)
        for a in range(len(pts_list)):
            for b in range(a+1,len(pts_list)):
                for c in range(b+1,len(pts_list)):
                    x1,y1=pts_list[a]; x2,y2=pts_list[b]; x3,y3=pts_list[c]
                    if (x2-x1)*(y3-y1)==(x3-x1)*(y2-y1): coll=True; break
                if coll: break
            if coll: break
        print(f'Collinear: {"✅ NONE" if not coll else "❌ FOUND"}')
        
        if not coll:
            print('\n🎉🎉🎉 n=74 C₄ SOLUTION FOUND!')
    else:
        print('  UNSAT')
        
except ImportError as e:
    print(f'PySAT not available: {e}')
    print('Writing DIMACS CNF file...')
    
    # Deduplicate clauses
    clause_set = set()
    for c in clauses:
        clause_set.add(tuple(sorted(c)))
    
    filename = 'n74_c4.cnf'
    with open(filename, 'w') as f:
        f.write(f'p cnf {n_vars} {len(clause_set)}\n')
        for clause in sorted(clause_set):
            f.write(' '.join(str(l) for l in clause) + ' 0\n')
    print(f'  Written {len(clause_set)} clauses to {filename}')

print('\nDone.')
