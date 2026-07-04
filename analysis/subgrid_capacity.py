"""
Level 2: Missing-center solvability criterion via subgrid parity
Level 1: Refined regression with per-subgrid 4k+1 prime features
"""

from collections import Counter
import numpy as np

# ── r2(d) computation ──
def factorize(n_):
    factors = {}
    d = 2
    while d * d <= n_:
        while n_ % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n_ //= d
        d += 1 if d == 2 else 2
    if n_ > 1: factors[n_] = 1
    return factors

def r2(d):
    if d == 0: return 1
    factors = factorize(d)
    for p, e in factors.items():
        if p % 4 == 3 and e % 2 == 1: return 0
    result = 4
    for p, e in factors.items():
        if p % 4 == 1: result *= (e + 1)
    return result

def count_4k1(d):
    return sum(1 for p in factorize(d) if p % 4 == 1)

def analyze_rings(n):
    ctr = (n-1)/2.0
    rings = Counter()
    for i in range(n):
        for j in range(n):
            d2 = int(round((i-ctr)**2 + (j-ctr)**2))
            rings[d2] += 1
    
    mod_cnt = {0:0,1:0,2:0}
    mod_cap = {0:0,1:0,2:0}
    max_4k1_mod = {0:0,1:0,2:0}
    
    for d2, pop in rings.items():
        m = d2 % 4
        if m == 3: continue
        mod_cnt[m] += 1
        mod_cap[m] += min(2, pop)
        nk = count_4k1(d2)
        max_4k1_mod[m] = max(max_4k1_mod[m], nk)
    
    total_cap = sum(mod_cap.values())
    return mod_cnt, mod_cap, max_4k1_mod, total_cap, len(rings)

# ── Data ──
missing_data = {
    7:(132,4),8:(380,0),9:(365,8),10:(1135,0),11:(1120,36),
    12:(4348,52),13:(3622,292),14:(10568,84),15:(30634,2716),
    16:(46304,1392),17:(55573,3872),18:(152210,24),19:(258170,10280),
    20:(941580,112),21:(9701,768),22:(5082,52),23:(16007,796),
    24:(11630,456),25:(36147,1804),26:(19724,584),27:(69513,1388),
    28:(48696,4),29:(179534,0),30:(99514,4960)
}

# ════════════════════════════════════════════
# LEVEL 2
# ════════════════════════════════════════════
print('LEVEL 2: MISSING-CENTER SOLVABILITY CRITERION')
print('=' * 70)
h = '{:>3} {:>5} {:>5} {:>4} {:>4} {:>4} {:>5} {:>5} {:>5} {:>8} {:>6} {:>4} {:>8}'
print(h.format('n','type','Rings','R0','R1','R2','Cap0','Cap1','Cap2','TotCap','Need','OK?','Missing'))
print('-' * 70)

for n in range(6, 31):
    mc, mcap, m4k1, total_cap, nrings = analyze_rings(n)
    need = 2 * n
    ptype = '4k+{}'.format(n%4) + ('E' if n%2==0 else 'O')
    feasible = 'YES' if total_cap >= need else 'NO'
    has_m = 'YES' if n in missing_data and missing_data[n][1] > 0 else 'NO'
    print(h.format(n, ptype, nrings, mc[0], mc[1], mc[2],
          mcap[0], mcap[1], mcap[2], total_cap, need, feasible, has_m))

# Even n threshold
print('\nEven n threshold:')
for n in range(6, 31, 2):
    _, _, _, total_cap, _ = analyze_rings(n)
    ratio = total_cap / (2*n)
    has_m = 'YES' if n in missing_data and missing_data[n][1] > 0 else 'NO'
    flag = ' <<<' if n == 12 else ''
    print('  n={:>2}: cap/need={:.3f} missing={}{}'.format(n, ratio, has_m, flag))

# ════════════════════════════════════════════
# LEVEL 1
# ════════════════════════════════════════════
print('\n\nLEVEL 1: REFINED REGRESSION')
print('=' * 70)

X, y, w = [], [], []
for n in sorted(missing_data.keys()):
    ft, fm = missing_data[n]
    mc, mcap, m4k1, total_cap, nrings = analyze_rings(n)
    X.append([nrings, mc[0], mc[1], mc[2], m4k1[0], m4k1[1], m4k1[2],
              mcap[0], mcap[1], mcap[2]])
    y.append(fm / ft * 100)
    w.append(ft)

X = np.array(X); y = np.array(y); w = np.array(w)
W = np.diag(w / w.mean())

# Model A: aggregate
XA = np.column_stack([np.ones(len(X)), X[:,0], np.max(X[:,4:7], axis=1)])
betaA = np.linalg.solve(XA.T @ W @ XA, XA.T @ W @ y)
r2_A = 1 - np.sum(W@((y-XA@betaA)**2))/np.sum(W@((y-np.average(y,weights=w))**2))
print('Model A (aggregate):')
for i,nm in enumerate(['intercept','rings','max4k1']):
    print('  {}: {:.4f}'.format(nm, betaA[i]))
print('  R2 = {:.4f}'.format(r2_A))

# Model B: per-mod4 rings + per-mod4 max_4k1
XB = np.column_stack([np.ones(len(X)), X[:,1], X[:,2], X[:,3], X[:,4], X[:,5], X[:,6]])
betaB = np.linalg.solve(XB.T @ W @ XB, XB.T @ W @ y)
r2_B = 1 - np.sum(W@((y-XB@betaB)**2))/np.sum(W@((y-np.average(y,weights=w))**2))
print('\nModel B (per-mod4 rings + per-mod4 max4k1):')
for i,nm in enumerate(['intercept','R0','R1','R2','4k1_0','4k1_1','4k1_2']):
    print('  {}: {:.4f}'.format(nm, betaB[i]))
print('  R2 = {:.4f}'.format(r2_B))

# Model C: per-mod4 capacities + per-mod4 max_4k1
XC = np.column_stack([np.ones(len(X)), X[:,7], X[:,8], X[:,9], X[:,4], X[:,5], X[:,6]])
betaC = np.linalg.solve(XC.T @ W @ XC, XC.T @ W @ y)
r2_C = 1 - np.sum(W@((y-XC@betaC)**2))/np.sum(W@((y-np.average(y,weights=w))**2))
print('\nModel C (per-mod4 capacities + per-mod4 max4k1):')
for i,nm in enumerate(['intercept','cap0','cap1','cap2','4k1_0','4k1_1','4k1_2']):
    print('  {}: {:.4f}'.format(nm, betaC[i]))
print('  R2 = {:.4f}'.format(r2_C))

print('\n{:>10} {:>8}'.format('Model','R2'))
print('{:>10} {:>8.4f}'.format('Old',0.620))
print('{:>10} {:>8.4f}'.format('A',r2_A))
print('{:>10} {:>8.4f}'.format('B',r2_B))
print('{:>10} {:>8.4f}'.format('C',r2_C))

# Best model predictions
best = 'B' if r2_B >= r2_C else 'C'
Xbest = XB if r2_B >= r2_C else XC
beta = betaB if r2_B >= r2_C else betaC
y_pred = Xbest @ beta
print('\nBest model {} predictions:'.format(best))
print('{:>3} {:>7} {:>7} {:>7}'.format('n','actual','pred','diff'))
for i, n in enumerate(sorted(missing_data.keys())):
    a, p = y[i], y_pred[i]
    print('{:>3} {:>6.2f}% {:>6.2f}% {:>+6.2f}%'.format(n, a, p, a-p))
