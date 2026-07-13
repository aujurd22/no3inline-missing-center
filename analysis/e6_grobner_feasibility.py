"""
E6: 代数几何路线可行性实验
----------------------
目标：验证 sympy 的 Gröbner 基能否处理 rot4-NTIL 的小规模多项式系统。
步骤：
  Phase 0: m=5, GF(7), 仅 2-因子约束（线性等式）→ Gröbner 基有多大/多快
  Phase 1: m=5, 加一条共线约束 → 系统有无解？
  Phase 2: m=5, 加全部共线约束 → Gröbner 基能否在合理时间内找到解？
  Phase 3: m=7/11, 观察缩放
"""

import sympy as sp
import time, sys, json, math

# === 几何辅助 ===
def c4(p, r, N):
    x, y = p
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    if r == 3: return (y, N - 1 - x)

def make_det_poly(i, j, k, cells, N, vars_dict):
    """生成三元组 (i,j,k) 是否共线的判定式（行列式）。
    返回多项式列表: [det0, det1, ...] 对应各旋转组合。
    若任何 det = 0 则三点共线。
    """
    polys = []
    for r in range(4):
        for s in range(4):
            for t in range(4):
                # (x1,y1) = c4(cells[i], r, N)
                xi, yi = cells[i]
                x1 = sp.sympify(c4((xi, yi), r, N)[0])
                y1 = sp.sympify(c4((xi, yi), r, N)[1])
                if isinstance(cells[j], tuple) and len(cells[j]) == 2 and all(isinstance(v, int) for v in cells[j]):
                    xj_val, yj_val = cells[j]
                    x2 = sp.sympify(c4((xj_val, yj_val), s, N)[0])
                    y2 = sp.sympify(c4((xj_val, yj_val), s, N)[1])
                else:
                    # variable-based
                    xj_val, yj_val = cells[j]
                    x2 = sp.sympify(c4((sp.sympify(xj_val), sp.sympify(yj_val)), s, N)[0])
                    y2 = sp.sympify(c4((sp.sympify(xj_val), sp.sympify(yj_val)), s, N)[1])
                # (x3,y3) = c4(cells[k], t, N)
                if isinstance(cells[k], tuple) and len(cells[k]) == 2 and all(isinstance(v, int) for v in cells[k]):
                    xk_val, yk_val = cells[k]
                    x3 = sp.sympify(c4((xk_val, yk_val), t, N)[0])
                    y3 = sp.sympify(c4((xk_val, yk_val), t, N)[1])
                else:
                    xk_val, yk_val = cells[k]
                    x3 = sp.sympify(c4((sp.sympify(xk_val), sp.sympify(yk_val)), t, N)[0])
                    y3 = sp.sympify(c4((sp.sympify(xk_val), sp.sympify(yk_val)), t, N)[1])
                # det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
                det = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
                polys.append(sp.simplify(det))
    return polys

# === 测试管理器 ===
results = []

def run_phase(name, m, domain='ZZ', timelimit=120):
    """通用相位测试：建变量、建约束、算 Gröbner 基，报告时间与规模"""
    print(f"\n{'='*60}")
    print(f"Phase: {name}  (m={m}, domain={domain})")
    print(f"{'='*60}")
    sys.stdout.flush()
    
    N = 2 * m
    # 变量：对于 permutation [m] 单环，cell (i, yi) 即 xi=i, yi=variable
    # 用 permutation matrix sel 的布尔变量太贵，这里用简化版：
    # 37×37=1369 布尔变量 → 不可行。改用 permutation 编码（m 个变量）
    if domain == 'ZZ':
        # 整数环上 Gröbner→变量是 ZZ
        y = sp.symbols(f'y0:{m}', integer=True)
    else:
        y = sp.symbols(f'y0:{m}')
    
    # 2-因子约束（permutation）：每个 row 恰好选一个，每个 col 恰好被选一个
    # 即 y_i 构成一个双射。用 1:1 约束：
    # 对于 perm: 每个 y_i 取 0..m-1 各一次
    # 多项式约束：∏_{i<j} (y_i - y_j) = 0 太贵 (degree m*(m-1)/2)
    # 改用弱约束：每个 y_i 在 {0..m-1} 中：∏_{t=0}^{m-1} (y_i - t) = 0
    constraints = []
    for i in range(m):
        # y_i 取值约束
        poly = 1
        for t in range(m):
            poly *= (y[i] - t)
        constraints.append(sp.expand(poly))
    
    n_var = m
    n_con = len(constraints)
    
    print(f"  Variables: {n_var}")
    print(f"  Constraints: {n_con} (domain constraints)")
    sys.stdout.flush()
    
    t0 = time.time()
    try:
        # 尝试计算 Gröbner 基
        G = sp.groebner(constraints, *y, order='lex', modulus=None if domain == 'ZZ' else int(domain.lstrip('GF')))
        dt = time.time() - t0
        
        n_terms = sum(len(g.free_symbols) + g.count_ops() for g in G)
        
        result = {
            'phase': name,
            'm': m,
            'domain': domain,
            'n_var': n_var,
            'n_con': n_con,
            'time_s': round(dt, 2),
            'gbasis_len': len(G),
            'gbasis_terms': n_terms,
            'status': 'OK',
        }
        print(f"  Gröbner basis: {len(G)} polys, ~{n_terms} terms, {dt:.2f}s")
        for i, g in enumerate(G[:5]):
            print(f"    g{i}: {g}")
        if len(G) > 5:
            print(f"    ... and {len(G)-5} more")
        
    except Exception as e:
        dt = time.time() - t0
        result = {
            'phase': name,
            'm': m,
            'domain': domain,
            'n_var': n_var,
            'n_con': n_con,
            'time_s': round(dt, 2),
            'gbasis_len': 0,
            'gbasis_terms': 0,
            'status': f'ERROR: {e}',
        }
        print(f"  ERROR ({dt:.1f}s): {e}")
    
    results.append(result)
    sys.stdout.flush()
    return result


def run_e6_scoped(m, timelimit=120):
    """简化版 E6 实验：只检查 collinear-free 条件的代数可行性"""
    print(f"\n{'='*60}")
    print(f"E6 scoped: m={m}, 单 [m] 环型, 枚举共线约束")
    print(f"{'='*60}")
    sys.stdout.flush()
    
    N = 2 * m
    # 变量
    y = sp.symbols(f'y0:{m}')
    
    # 已知 m=5 的一个解：[0->0, 1->2, 2->4, 3->1, 4->3] 
    # 把已知解代入验证
    if m == 5:
        known_sol = [0, 2, 4, 1, 3]  # π(0)=0, π(1)=2, ...
        cells = [(i, known_sol[i]) for i in range(m)]
        
        # 构建所有共线判定式，代入已知解，验证是否无零点
        bad_count = 0
        t0 = time.time()
        for i in range(m):
            for j in range(i + 1, m):
                for k in range(j + 1, m):
                    for r in range(4):
                        for s in range(4):
                            for t_rot in range(4):
                                x1, y1 = c4(cells[i], r, N)
                                x2, y2 = c4(cells[j], s, N)
                                x3, y3 = c4(cells[k], t_rot, N)
                                det = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
                                if det == 0:
                                    bad_count += 1
                                    if bad_count <= 3:
                                        print(f"  BAD: (i={i},j={j},k={k},r={r},s={s},t={t_rot}): det=0")
        dt = time.time() - t0
        print(f"  m={m} known solution: bad_events={bad_count}, {dt:.3f}s")
        print(f"  (expect 0 for a valid solution)")
        sys.stdout.flush()
    
    # 现在测试代数编码：对于一个未知的 π，collinear-free 条件如何表示？
    # 每个三元组 (i,j,k) 需要避免所有 64 种旋转组合的共线
    # 约束：∏_{r,s,t} det(i,j,k,r,s,t) ≠ 0 (但对所有三元组同时成立)
    # 这会产生高次多项式，无法直接做 Gröbner
    
    # 折中：只取一个三元组 (0,1,2) 的 64 个行列式，看 Gröbner 能做什么
    print(f"\n  --- 微型实验: 3 个 cell, 64 个 det 多项式 ---")
    
    # 固定前 3 行的 π 值 (k=0,1,2)
    vars = [y[0], y[1], y[2]]
    cells_v = [(0, y[0]), (1, y[1]), (2, y[2])]
    
    # 生成所有 64 个 det 多项式
    polys = []
    for r in range(4):
        for s in range(4):
            for t_rot in range(4):
                x1_val, y1_val = c4((0, y[0]), r, N)
                x2_val, y2_val = c4((1, y[1]), s, N)
                x3_val, y3_val = c4((2, y[2]), t_rot, N)
                det = x1_val*(y2_val - y3_val) + x2_val*(y3_val - y1_val) + x3_val*(y1_val - y2_val)
                det_s = sp.simplify(det)
                polys.append(det_s)
    
    # 去重
    unique_polys = list(set(str(p) for p in polys))
    print(f"  Unique det polynomials from 64 combos: {len(unique_polys)}")
    
    # 计算这 3 变量的 Gröbner 基
    for modulus in [None, 7, 13, 17]:
        label = 'ZZ' if modulus is None else f'GF({modulus})'
        print(f"\n  --- Gröbner on {label} ({len(vars)} vars, {len(polys)} eqs) ---")
        sys.stdout.flush()
        t0 = time.time()
        try:
            if modulus is None:
                # Show the polynomials first
                pass
            G = sp.groebner(polys[:8], *vars, order='lex', modulus=modulus)
            dt = time.time() - t0
            print(f"    Gröbner basis: {len(G)} polys, {dt:.3f}s")
            for i, g in enumerate(G):
                print(f"      g{i}: {g}")
        except Exception as e:
            dt = time.time() - t0
            print(f"    ERROR ({dt:.1f}s): {e}")
        sys.stdout.flush()
    
    print()
    return

if __name__ == '__main__':
    # Phase 0: m=5 域约束 Gröbner
    run_phase("Phase 0: domain constraints only", m=5, domain='GF(7)', timelimit=60)
    
    # Phase 0b: 整数环
    # run_phase("Phase 0b: domain on ZZ", m=5, domain='ZZ', timelimit=30)
    
    # E6 scoped: 微共线实验
    run_e6_scoped(m=5, timelimit=30)

    # 汇总
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    for r in results:
        print(f"  {r['phase']}: status={r['status']}, time={r['time_s']}s, "
              f"gbasis={r['gbasis_len']} polys, ~{r['gbasis_terms']} terms")
    
    with open('results/e6_feasibility.json', 'w') as f:
        json.dump(results, f, indent=1)
    print("\nSaved to results/e6_feasibility.json")
    sys.stdout.flush()
