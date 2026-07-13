"""
Dir5 (data-driven invariant mining), Phase-1: extract invariants from KNOWN solutions.

For every known solution m=05..19 (results/solutions/mXX.json), lift to 4m points,
classify every line through >=2 points by its D4-canonical direction, and record the
occupancy (number of points on that line). Since these are true solutions, all
occupancies are <=2. We aggregate across m to find:
  (a) which canonical directions carry the most occ=2 lines (the "critical directions"
      that real solutions must thread between),
  (b) the slope/det multiset structure,
  (c) stability of these patterns across m.

This is lightweight read+analysis only (no random sampling -> no CPU contention with
the running gating experiment). Phase-2 (random near-miss control) is deferred until
after gating finishes, reusing the same lift/line machinery.

Output: results/dir5_invariants.json
"""
import json, math, glob, os
from collections import Counter, defaultdict

D4 = [
    lambda dx, dy: (dx, dy),
    lambda dx, dy: (-dy, dx),
    lambda dx, dy: (-dx, -dy),
    lambda dx, dy: (dy, -dx),
    lambda dx, dy: (dx, -dy),
    lambda dx, dy: (-dx, dy),
    lambda dx, dy: (dy, dx),
    lambda dx, dy: (-dy, -dx),
]

def canon_dir(dx, dy):
    if dx == 0 and dy == 0:
        return None
    g = math.gcd(abs(dx), abs(dy))
    dx //= g; dy //= g
    if dx < 0 or (dx == 0 and dy < 0):
        dx = -dx; dy = -dy
    best = (dx, dy)
    for T in D4:
        a, b = T(dx, dy)
        ga = math.gcd(abs(a), abs(b))
        a //= ga; b //= ga
        if a < 0 or (a == 0 and b < 0):
            a = -a; b = -b
        if (a, b) < best:
            best = (a, b)
    return best

def c4(x, y, r, N):
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    return (y, N - 1 - x)

def lift(cells, N):
    pts = []
    for x, y in cells:
        for r in range(4):
            pts.append(c4(x, y, r, N))
    return pts

def line_id(p, q):
    dx = q[0] - p[0]; dy = q[1] - p[1]
    g = math.gcd(abs(dx), abs(dy)); a = dx // g; b = dy // g
    if a < 0 or (a == 0 and b < 0):
        a = -a; b = -b
    c = b * p[0] - a * p[1]
    return (a, b, c)

def analyze(m, cells):
    N = 2 * m
    pts = lift(cells, N)
    lp = defaultdict(list)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            lp[line_id(pts[i], pts[j])].append(i)
    dir_occ = defaultdict(Counter)   # canon_dir -> Counter(occupancy)
    for (a, b, c), idxs in lp.items():
        occ = len(set(idxs))
        d = canon_dir(a, b)
        dir_occ[d][occ] += 1
    return dir_occ, N

def main():
    sol_dir = "results/solutions"
    files = sorted(glob.glob(os.path.join(sol_dir, "m*.json")))
    per_m = {}
    global_dir_occ = defaultdict(Counter)   # canon_dir -> Counter(occ) summed over all m
    global_dir_lines = Counter()            # canon_dir -> total #lines (occ>=2)
    for f in files:
        m = int(os.path.basename(f)[1:3])
        data = json.load(open(f))
        cells = data["cells"]
        # some files use list-of-lists, some list-of-[x,y]
        cells = [(c[0], c[1]) for c in cells]
        dir_occ, N = analyze(m, cells)
        per_m[m] = {f"{d}": dict(cnt) for d, cnt in dir_occ.items()}
        for d, cnt in dir_occ.items():
            global_dir_occ[d].update(cnt)
            global_dir_lines[d] += sum(cnt.values())
    # summarize critical directions: occ=2 lines per canon_dir (across all m)
    crit = {f"{d}": global_dir_occ[d].get(2, 0) for d in global_dir_occ}
    crit_sorted = sorted(crit.items(), key=lambda kv: -kv[1])
    out = {
        "per_m": per_m,
        "global_dir_occ": {f"{d}": dict(cnt) for d, cnt in global_dir_occ.items()},
        "global_dir_total_lines": dict(global_dir_lines),
        "critical_directions_occ2_ranked": crit_sorted,
        "n_ms": len(per_m),
    }
    os.makedirs("results", exist_ok=True)
    tmp = "results/dir5_invariants.json.tmp"
    json.dump(out, open(tmp, "w"), indent=1)
    os.replace(tmp, "results/dir5_invariants.json")
    # console brief
    print(f"[dir5] analyzed {len(per_m)} solutions (m={sorted(per_m)})")
    print(f"[dir5] distinct canonical directions seen: {len(global_dir_occ)}")
    print("[dir5] top-10 critical directions (occ=2 line count across all m):")
    for d, c in crit_sorted[:10]:
        print(f"    {d}: {c}")
    print("[dir5] wrote results/dir5_invariants.json")

if __name__ == "__main__":
    main()
