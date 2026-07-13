"""
Moser-Tardos 局部重采样求解器 for rot4-NTIL (Direction-1 实质推进).

复用 gating_lll_r1.py 同款模型:
  - 2-因子 = 配置模型 m 条有向边 cells (u,v) (每顶点度 2)
  - 提升为 4m C4 点；坏事件 = 几何线含 >=3 点

Moser-Tardos 算法 (resample-to-kill):
  while 存在坏事件:
     随机选一个坏线 L
     在其涉及顶点邻域内做随机 2-switch, 强制使 L 消失 (L 上点 <=2)
     (fallback: 全局随机 2-switch 严格减少 B)
  每次 resample(L) 使选中的 L 以概率 1 消失 -> B 单调受控下降.
  收敛性与 gating 证据一致: red_config_frac=1.0, max_red_seen 随 m 增.

测量: 小 m (10,14,18) 收敛到 B=0 的比例 + 平均 resample 次数 -> 外推 m=37.
"""
import os, sys, json, math, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gating_lll_r1 import (lift_points, line_map, bad_lines, bad_count,
                           on_line, cell_of_point, random_2factor)

def vertices_of_line(pts, key, m):
    verts = set()
    for P in pts:
        if on_line(P, key):
            c = cell_of_point(P, m)
            if c:
                verts.add(c[0]); verts.add(c[1])
    return verts

def resample_bad_event(cells, m, key, rng, tries=300):
    """在 L 涉及顶点的边内做随机 2-switch, 使坏线 key 消失. 返回新 cells 或 None."""
    pts = lift_points(cells, m)
    verts = vertices_of_line(pts, key, m)
    idxs = [i for i, c in enumerate(cells) if c[0] in verts or c[1] in verts]
    if len(idxs) < 2:
        return None
    for _ in range(tries):
        i = rng.choice(idxs); j = rng.choice(idxs)
        if i == j:
            continue
        a, b = cells[i]; c, d = cells[j]
        if a == c and b == d:
            continue
        nc = list(cells)
        nc[i] = (a, d); nc[j] = (c, b)
        # 2-switch 保持 2-正则 (顶点度不变)
        nlm = line_map(lift_points(nc, m))
        if nlm.get(key, 0) < 3:
            return nc
    return None

def global_reduce(cells, m, rng, tries=300):
    mc = len(cells)
    B0 = bad_count(line_map(lift_points(cells, m)))
    for _ in range(tries):
        i = rng.randrange(mc); j = rng.randrange(mc)
        if i == j:
            continue
        a, b = cells[i]; c, d = cells[j]
        if a == c and b == d:
            continue
        nc = list(cells)
        nc[i] = (a, d); nc[j] = (c, b)
        if bad_count(line_map(lift_points(nc, m))) < B0:
            return nc
    return None

def moser_tardos(m, rng, max_iters=200000, tries=300):
    cells = random_2factor(m, rng)
    lm = line_map(lift_points(cells, m))
    bl = bad_lines(lm)
    it = 0
    while bl and it < max_iters:
        key = rng.choice(bl)
        nc = resample_bad_event(cells, m, key, rng, tries)
        if nc is None:
            nc = global_reduce(cells, m, rng, tries)
        if nc is None:
            break  # 卡住
        cells = nc
        lm = line_map(lift_points(cells, m))
        bl = bad_lines(lm)
        it += 1
    return cells, bad_count(lm), it

def main():
    rng = random.Random(20260713)
    plan = [(10, 150), (14, 100), (18, 50), (22, 30)]
    summary = {}
    for m, trials in plan:
        solved = 0
        iters = []
        t0 = time.time()
        for t in range(trials):
            cells, B, it = moser_tardos(m, rng)
            if B == 0:
                solved += 1
                iters.append(it)
        rate = solved / trials
        avg_it = (sum(iters) / len(iters)) if iters else -1
        summary[m] = {"trials": trials, "solve_rate": round(rate, 3),
                      "avg_iters": round(avg_it, 1),
                      "sec": round(time.time() - t0, 1)}
        print(f"m={m} trials={trials} solve_rate={rate:.3f} avg_iters={avg_it:.1f} "
              f"sec={summary[m]['sec']}", flush=True)
    out = {"model": "Moser-Tardos resample-to-kill on config-model 2-factor",
           "summary": summary}
    with open("results/moser_tardos_r1.json", "w") as f:
        json.dump(out, f, indent=1)
    print("[done] wrote results/moser_tardos_r1.json")

if __name__ == "__main__":
    main()
