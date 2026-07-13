"""
gating_lll_r1.py -- Direction-1 LLL / Switching 可行性 gating 实验.

模型 (R9b / Th-44):
  rot4-NTIL 的 m 个 cell 构成一个有向 2-正则图 (顶点 0..m-1),
  indeg(i)+outdeg(i)==2, 即 `rowSum[i]+colSum[i]==2`.
  随机 2-因子 = 配置模型: 2m 个 stub (每顶点 2 个) 随机配对 -> m 条有向边
  (cell (u,v) = 边 u->v). 允许 loop (u,u) 与 double-edge (2-环).

  每个 cell (x,y) 提升为 4 个 C4 旋转点 -> 共 4m 点. 坏事件 A_L = "几何线 L
  含 >=3 个提升点".

测量 (评审的止损标准):
  (1) E[B]             期望坏线数 (对拍 R9e 的 ~10.3m)
  (2) N, pbar=E[B]/N   候选坏事件数 & 平均逐事件概率
  (3) d (共享 2-因子顶点 的依赖度) & e*pbar*(d+1)  (vanilla LLL 量及其失败判据)
  (4) switch-reducibility  含坏构型中存在可严格减少 B 的 2-switch 的比例
                          = switching / 熵压缩 LLL 的真正可行性判据

输出: results/gating_lll_r1.json + 打印摘要.
"""
import os, sys, json, math, random
from collections import defaultdict
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from quadratic_sidon_completeness import c4

# ───────────────────────── geometry ─────────────────────────
def lift_points(cells, m):
    n = 2 * m
    pts = set()
    for (x, y) in cells:
        for r in range(4):
            pts.add(c4((x, y), r, n))
    return list(pts)

def line_map(pts):
    """key=(rdx,rdy,intercept) -> pair-count on that line. count>=3 <=> >=3 points."""
    mm = {}
    L = len(pts)
    for i in range(L):
        x1, y1 = pts[i]
        for j in range(i + 1, L):
            x2, y2 = pts[j]
            dx, dy = x2 - x1, y2 - y1
            g = gcd(abs(dx), abs(dy)) or 1
            rdx, rdy = dx // g, dy // g
            if rdx < 0 or (rdx == 0 and rdy < 0):
                rdx, rdy = -rdx, -rdy
            perp = (-rdy, rdx)
            inter = perp[0] * x1 + perp[1] * y1
            k = (rdx, rdy, inter)
            mm[k] = mm.get(k, 0) + 1
    return mm

def bad_lines(lm):
    return [k for k, c in lm.items() if c >= 3]

def bad_count(lm):
    return sum(1 for c in lm.values() if c >= 3)

def on_line(P, key):
    rdx, rdy, inter = key
    return (-rdy) * P[0] + rdx * P[1] == inter

def cell_of_point(P, m):
    for r in range(4):
        c = c4(P, r, 2 * m)
        if 0 <= c[0] < m and 0 <= c[1] < m:
            return c
    return None

# ───────────────────────── sampler ─────────────────────────
def random_2factor(m, rng):
    stubs = []
    for v in range(m):
        stubs.append(v)
        stubs.append(v)
    rng.shuffle(stubs)
    return [(stubs[k], stubs[k + 1]) for k in range(0, 2 * m, 2)]

# ───────────────────────── gating ─────────────────────────
def gating_one(m, rng, n_cfg, n_red, K):
    # ---- Phase 1: candidate universe + E[B] + dependency d ----
    candidates = {}          # key -> #configs where bad
    Bs = []
    dsum = 0
    dcount = 0
    dmax = 0
    pLs = []                 # per-line bad probability estimate (fraction)
    for s in range(n_cfg):
        cells = random_2factor(m, rng)
        pts = lift_points(cells, m)
        lm = line_map(pts)
        bl = bad_lines(lm)
        Bs.append(len(bl))
        # touched vertices per bad line (in this config)
        line_verts = {}
        for k in bl:
            verts = set()
            for P in pts:
                if on_line(P, k):
                    c = cell_of_point(P, m)
                    if c:
                        verts.add(c[0])
                        verts.add(c[1])
            line_verts[k] = verts
        for k in bl:
            candidates[k] = candidates.get(k, 0) + 1
            vk = line_verts[k]
            deg = 0
            for k2 in bl:
                if k2 != k and line_verts[k2] & vk:
                    deg += 1
            dsum += deg
            dcount += 1
            if deg > dmax:
                dmax = deg
    E_B = sum(Bs) / len(Bs)
    N = len(candidates)
    pbar = E_B / N if N else 0.0
    d_mean = dsum / dcount if dcount else 0.0
    # per-line p_L (= fraction of configs where bad), aggregated later via candidates
    # e * pbar * (d_mean+1)  (vanilla LLL quantity)
    e_p_d = math.e * pbar * (d_mean + 1) if d_mean else 0.0

    # ---- Phase 2: switch-reducibility (decisive for switching/entropy LLL) ----
    red_configs = 0
    total_Bgt0 = 0
    red_fracs = []
    max_red_seen = 0
    for s in range(n_red):
        cells = random_2factor(m, rng)
        pts = lift_points(cells, m)
        B0 = bad_count(line_map(pts))
        if B0 == 0:
            continue
        total_Bgt0 += 1
        mc = len(cells)
        reducing = 0
        maxred = 0
        for _ in range(K):
            i = rng.randrange(mc)
            j = rng.randrange(mc)
            if i == j:
                continue
            a, b = cells[i]
            c, d = cells[j]
            if a == c and b == d:
                continue
            new_cells = list(cells)
            new_cells[i] = (a, d)
            new_cells[j] = (c, b)
            B1 = bad_count(line_map(lift_points(new_cells, m)))
            red = B0 - B1
            if red > 0:
                reducing += 1
                if red > maxred:
                    maxred = red
        if reducing > 0:
            red_configs += 1
        red_fracs.append(reducing / K)
        if maxred > max_red_seen:
            max_red_seen = maxred

    red_frac_configs = (red_configs / total_Bgt0) if total_Bgt0 else 0.0
    mean_red_frac = (sum(red_fracs) / len(red_fracs)) if red_fracs else 0.0
    return {
        "m": m,
        "n_cfg": n_cfg,
        "n_red": n_red,
        "K": K,
        "E_B": E_B,
        "N_candidates": N,
        "pbar": pbar,
        "d_mean": d_mean,
        "d_max": dmax,
        "e_p_d_vanilla": e_p_d,
        "red_config_frac": red_frac_configs,
        "mean_red_frac": mean_red_frac,
        "max_red_seen": max_red_seen,
    }

def main():
    rng = random.Random(20260713)
    # (m, n_cfg, n_red, K) adaptive by size
    plan = [
        (10, 600, 400, 600),
        (14, 500, 300, 500),
        (18, 400, 250, 400),
        (22, 300, 200, 300),
        (26, 250, 150, 250),
        (30, 200, 120, 200),
    ]
    results = []
    for (m, nc, nr, K) in plan:
        print(f"[gating] m={m}  n_cfg={nc} n_red={nr} K={K} ...", flush=True)
        t0 = __import__("time").time()
        r = gating_one(m, rng, nc, nr, K)
        r["sec"] = round(__import__("time").time() - t0, 1)
        results.append(r)
        print(f"  E_B={r['E_B']:.2f}  N={r['N_candidates']}  pbar={r['pbar']:.4g}  "
              f"d_mean={r['d_mean']:.1f}  e*p*(d+1)={r['e_p_d_vanilla']:.3g}  "
              f"red_config_frac={r['red_config_frac']:.3f}  mean_red_frac={r['mean_red_frac']:.3f}  "
              f"[{r['sec']}s]", flush=True)
    out = {
        "model": "config-model random 2-factor (2m stubs paired)",
        "notes": "e_p_d_vanilla = e*pbar*(d_mean+1) using vertex-sharing dependency; "
                 "vanilla LLL needs <=1 (fails by orders). red_config_frac = fraction of "
                 "bad configs having >=1 reducing 2-switch (switching/entropy-compression viability).",
        "results": results,
    }
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "gating_lll_r1.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[done] wrote {path}")

if __name__ == "__main__":
    main()
