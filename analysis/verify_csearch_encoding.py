import json
from collections import defaultdict

d = json.load(open('results/solutions/m10.json'))
cells = d['cells']; m = 10; n = 2 * m

def c4(x, y, r, N):
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    return (y, N - 1 - x)

def igcd(a, b):
    a = abs(a); b = abs(b)
    while b: a, b = b, a % b
    return a

# replicate csearch's line_of exactly
def line_of(p, q):
    a = -(q[1] - p[1]); b = (q[0] - p[0]); c = (q[1] - p[1]) * p[0] - (q[0] - p[0]) * p[1]
    g = igcd(abs(a), abs(b)); g = igcd(g, abs(c))
    if g != 0: a //= g; b //= g; c //= g
    if a < 0 or (a == 0 and b < 0) or (a == 0 and b == 0 and c < 0): a = -a; b = -b; c = -c
    return (a, b, c)

# csearch-style lift (list, no dedup)
pts = []
for (x, y) in cells:
    for r in range(4): pts.append(c4(x, y, r, n))
print('npts', len(pts), 'unique', len(set(pts)))

g = defaultdict(int)
for i in range(len(pts)):
    for j in range(i):
        if pts[i] == pts[j]: continue
        g[line_of(pts[i], pts[j])] += 1
bad = 0
for k, c in g.items():
    if c >= 3: bad += c * (c - 1) * (c - 2) // 6
print('bad(line_of, csearch logic) =', bad)

def brute(pl):
    for i in range(len(pl)):
        x1, y1 = pl[i]
        for j in range(i + 1, len(pl)):
            x2, y2 = pl[j]; dx, dy = x2 - x1, y2 - y1
            for k in range(j + 1, len(pl)):
                x3, y3 = pl[k]
                if dx * (y3 - y1) == dy * (x3 - x1): return True
    return False

print('brute(pts, csearch lift) has-bad =', brute(pts))

# true verify_cells path: pairs -> cell_xy -> c4 (uses set, dedup)
def cell_xy(a, b, mm): return (mm - (a + 1) // 2, mm - (b + 1) // 2)
pairs = [(2 * (m - x) - 1, 2 * (m - y) - 1) for (x, y) in cells]
lifted = set()
for (a, b) in pairs:
    x, y = cell_xy(a, b, m)
    for r in range(4): lifted.add(c4((x, y), r, n))
print('true-verify lifted npts', len(lifted), 'brute-has-bad =', brute(list(lifted)))
