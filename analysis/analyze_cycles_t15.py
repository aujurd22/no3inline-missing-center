"""T15 support: cycle-type statistics + mechanical verification of two claims.
Reads results/cycle_analysis.json (all known rot4-NTIL solutions, m=3.. up to 36).
Outputs results/cycle_type_stats.json + prints a summary.
"""
import json, itertools, random, os

ROOT = "D:/djr82/Documents/workbuddy/2026-07-03-16-29-36/no3inline-rigidity/analysis"
DATA = os.path.join(ROOT, "results/cycle_analysis.json")
STATS = os.path.join(ROOT, "results/cycle_type_stats.json")

with open(DATA) as f:
    blob = json.load(f)
data = blob["data"]

rows = []
all_lengths = set()
len2_present = False
m_with_m = 0
m_with_1m1 = 0
share_num = 0
share_den = 0

for m_str, info in sorted(data.items(), key=lambda kv: int(kv[0])):
    m = int(m_str)
    nsol = info["nsol"]
    dt = info["distinct_types"]
    types = info["types"]                     # list of [signature, count]
    sigs = [tuple(t[0]) for t in types]
    counts = [t[1] for t in types]

    has2 = any(2 in sig for sig in sigs)
    if has2:
        len2_present = True
    has_m = any(sig == (m,) for sig in sigs)
    if has_m:
        m_with_m += 1
    has_1m1 = any(sig in ((1, m - 1), (m - 1, 1)) for sig in sigs)
    if has_1m1:
        m_with_1m1 += 1

    top_share = 0
    for sig, c in zip(sigs, counts):
        if sig == (m,) or sig == (1, m - 1):
            top_share += c
    share_num += top_share
    share_den += nsol

    for sig in sigs:
        for L in sig:
            all_lengths.add(L)

    maxL = max(max(sig) for sig in sigs)
    rows.append({
        "m": m, "nsol": nsol, "distinct_types": dt, "n_sigs": len(sigs),
        "has_len2": has2, "has_single_m": has_m, "has_1m1": has_1m1,
        "top2_share": round(top_share / nsol, 3) if nsol else None,
        "max_cycle": maxL,
        "sigs_sample": [list(s) for s in sigs[:6]],
    })

summary = {
    "m_range": [min(int(k) for k in data), max(int(k) for k in data)],
    "len2_present_anywhere": len2_present,
    "m_with_single_m_cycle": m_with_m,
    "m_with_1m1": m_with_1m1,
    "total_m": len(data),
    "overall_top2_share": round(share_num / share_den, 4) if share_den else None,
    "all_cycle_lengths_observed": sorted(all_lengths),
}

# ---- mechanical verification 1: 3-cycle is always baseframe-safe ----
# cells (x,y),(y,z),(z,x); baseframe coords a=2(m-u)-1.
# det = -1/2[(x-y)^2+(y-z)^2+(z-x)^2]; zero iff x=y=z (degenerate loop).
def det3(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1])

three_cycle_ok = True
for _ in range(5000):
    x, y, z = random.sample(range(20), 3)
    m = 20
    def bf(u, v):
        return (2 * (m - u) - 1, 2 * (m - v) - 1)
    d = det3(bf(x, y), bf(y, z), bf(z, x))
    if d == 0:                       # distinct vertices must NEVER be collinear
        three_cycle_ok = False
        break

# ---- mechanical verification 2: mutual edge (2-cycle) self-collinearity ----
# a mutual edge (x,y),(y,x) lifts to 8 points (4 orientations each cell).
# check whether ANY triple among these 8 is collinear for x != y.
def lifts(m, x, y):
    a = 2 * (m - x) - 1
    b = 2 * (m - y) - 1
    return [(a, b), (-b, a), (-a, -b), (b, -a),
            (b, a), (-a, b), (-b, -a), (a, -b)]

mutual_self_collinear = False
random.seed(0)
for _ in range(20000):
    m = random.randint(3, 40)
    x = random.randint(0, m - 1)
    y = random.randint(0, m - 1)
    if x == y:
        continue
    if x + y >= 2 * m - 1:          # a=-b impossible inside grid -> skip
        continue
    pts = lifts(m, x, y)
    found = False
    for trio in itertools.combinations(range(8), 3):
        if det3(pts[trio[0]], pts[trio[1]], pts[trio[2]]) == 0:
            found = True
            break
    if found:
        mutual_self_collinear = True
        break

# sanity: det3 detects an obviously collinear triple
sanity = det3((0, 0), (1, 1), (2, 2)) == 0

out = {
    "summary": summary,
    "per_m": rows,
    "verify_3cycle_distinct_vertices_never_collinear": three_cycle_ok,
    "verify_mutual_edge_no_self_collinear_triple": (not mutual_self_collinear),
    "verify_det3_sanity": sanity,
}
with open(STATS, "w") as f:
    json.dump(out, f, indent=1)

print(json.dumps(summary, indent=1))
print("3cycle_safe =", three_cycle_ok,
      "| mutual_edge_self_safe =", not mutual_self_collinear,
      "| det3_sanity =", sanity)
