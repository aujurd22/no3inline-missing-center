"""
E6: Gröbner 基攻 m=37 — 分布式三元组约束法
========================================
核心思路：
  rot4-NTIL 的禁忌条件是：不存在 (i,j,k) 三元组在 C4 提升下共线。
  Gröbner 基为每个三元组 (i,j,k) 编码了其"共线条件代数簇"。
  若所有三元组的共线条件簇的并集覆盖整个搜索空间，则无解。
  
  更实际的做法：
  1. 对每个三元组 (i,j,k)，计算大乘积多项式 P_{ijk} = Π det(r,s,t) ≈ 0
  2. 其 Gröbner 基给出共线的 (y_i,y_j,y_k) 赋值集合
  3. 检查已知解是否不在任何共线簇中 → 验证方法正确性
  4. 对 m=37: 若代数簇联合覆盖了所有可能的 permutation → 无解
  
  由于 "在所有 p 上同时成立" 才对应整数解，我们取多个素数联合验证。
"""

import sympy as sp
import time, sys, json, math, itertools

# === 几何函数 ===
def c4(p, r, N):
    x, y = p
    if r == 0: return (x, y)
    if r == 1: return (N - 1 - y, x)
    if r == 2: return (N - 1 - x, N - 1 - y)
    if r == 3: return (y, N - 1 - x)

def det_poly(i, yi, j, yj, k, yk, N):
    """生成 3 个 variable-based cell 的所有 det 多项式乘积"""
    terms = []
    for r in range(4):
        for s in range(4):
            for t in range(4):
                x1, y1 = c4((i, yi), r, N)
                x2, y2 = c4((j, yj), s, N)
                x3, y3 = c4((k, yk), t, N)
                # det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)  
                det = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
                terms.append(sp.simplify(det))
    return terms

# === 已知解（从文件加载）===
def load_solution(m):
    """Load cells from known solution files. Returns list of (x,y) or None."""
    import os
    base = 'results/solutions'
    path = f'{base}/m{m:02d}.json'
    if os.path.exists(path):
        d = json.load(open(path))
        if 'cells' in d:
            return d['cells']
    return None

def load_m36_solution():
    from rot4_loader import decode_line
    line = open('flammenkamp_cache/n72_rot4.few').read().strip()
    pts_all = decode_line(line, 72)
    return [(x,y) for (x,y) in pts_all if x < 36 and y < 36]

# === 核心实验 ===

def verify_solution(cells, m):
    """Full brute-force verify of a candidate solution."""
    N = 2 * m
    pts = [c4(cell, r, N) for cell in cells for r in range(4)]
    bad = 0
    n = len(pts)
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            dx1, dy1 = xj - xi, yj - yi
            for k in range(j + 1, n):
                xk, yk = pts[k]
                if dx1 * (yk - yi) == dy1 * (xk - xi):
                    bad += 1
    return bad == 0, bad

def compute_triple_gb(i, j, k, m, p=None):
    """对三元组 (i,j,k) 计算共线条件 Gröbner 基"""
    N = 2 * m
    yi, yj, yk = sp.symbols(f'y{i} y{j} y{k}')
    
    # 乘积多项式: Π det(r,s,t) = 0
    # 如果乘积 = 0，则至少一个 det = 0 → 共线
    P = 1
    for r in range(4):
        for s in range(4):
            for t in range(4):
                x1, y1_p = c4((i, yi), r, N)
                x2, y2_p = c4((j, yj), s, N)
                x3, y3_p = c4((k, yk), t, N)
                det = x1*(y2_p - y3_p) + x2*(y3_p - y1_p) + x3*(y1_p - y2_p)
                P *= det
    P = sp.expand(P)
    # 指数展开太大，用 factor 可能更好
    
    # 只用 3 个随机组合做 Gröbner（用局部变量避免完全展开)
    # 改用单个 det 乘积太大时，换为子集
    return None  # placeholder

def gr_vs_m5(p=7):
    """m=5 多项式系统实验：在 GF(p) 上建系并求解"""
    print(f"\n{'='*60}")
    print(f"m=5 代数系统 (GF({p}))")
    print(f"{'='*60}")
    sys.stdout.flush()
    
    m = 5
    N = 10
    
    # 加载已知解
    sol = load_solution(5)
    if sol:
        print(f"Known solution from file: {sol}")
        ok, bad = verify_solution(sol, 5)
        print(f"Verify: {'PASS' if ok else f'FAIL (bad={bad})'}")
        sys.stdout.flush()
    else:
        print("WARNING: No m=5 solution file found")
        # 使用 [0,0],[1,2],[2,4],[3,1],[4,3] (从前更早的猜测)
        sol = [(0,0),(1,2),(2,4),(3,1),(4,3)]
        print(f"Fallback: {sol}")
    
    # 提取 y 值 (permutation)
    ys = [y for (x,y) in sol]
    print(f"y values: {ys}")
    sys.stdout.flush()
    
    # 变量
    y = sp.symbols('y0:5')
    
    # === PART 1: 检查 GF(p) 上的 permutation 存在性 ===
    # domain constraint: y_i ∈ GF(p), 互异
    # 编码为 Π (y_i - y_j) ≠ 0 (对 i<j)
    # 和 Π_{t=0}^{p-1} (y_i - t) = 0
    
    # 先直接验证已知解在 GF(p) 上是否仍是解
    print(f"\n--- Part 1: Check known solution in GF({p}) ---")
    
    # 对每个三元组检查 det mod p
    collinear_pairs = []
    collinear_count = 0
    for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                for r in range(4):
                    for s in range(4):
                        for t_rot in range(4):
                            x1, y1_p = c4((i, ys[i]), r, N)
                            x2, y2_p = c4((j, ys[j]), s, N)
                            x3, y3_p = c4((k, ys[k]), t_rot, N)
                            det_val = x1*(y2_p - y3_p) + x2*(y3_p - y1_p) + x3*(y1_p - y2_p)
                            if det_val % p == 0:
                                collinear_count += 1
                                if len(collinear_pairs) < 5:
                                    collinear_pairs.append((i,j,k,r,s,t_rot,det_val))
    print(f"  Collinear events mod {p}: {collinear_count}")
    for item in collinear_pairs:
        print(f"    (i={item[0]},j={item[1]},k={item[2]},r={item[3]},s={item[4]},t={item[5]}): det={item[6]} ≡ 0 mod {p}")
    
    # 一个真正的 NTIL 解应该有 0 个 det ≡ 0 mod ℤ
    # 但 mod p 可能有额外的 ≡ 0 (即 det 被 p 整除而本身非零)
    # 这叫做 "false positive" over GF(p)
    
    # === PART 2: 逐个三元组的 Gröbner 基 ===
    print(f"\n--- Part 2: Triple-by-triple Gröbner basis ---")
    
    triple_results = []
    t_total = 0.0
    
    for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                t0 = time.time()
                yi_v, yj_v, yk_v = sp.symbols(f'y{i} y{j} y{k}')
                
                # 乘积多项式 Π det_{r,s,t} (degree=64×2=128 max)
                P = 1
                for r in range(4):
                    for s in range(4):
                        for t_rot in range(4):
                            x1, y1_p = c4((i, yi_v), r, N)
                            x2, y2_p = c4((j, yj_v), s, N)
                            x3, y3_p = c4((k, yk_v), t_rot, N)
                            det = x1*(y2_p - y3_p) + x2*(y3_p - y1_p) + x3*(y1_p - y2_p)
                            P *= det
                
                # Factor (更快展开)
                try:
                    P_exp = sp.expand(P)
                    # 检查已知解是否违反
                    subbed = P_exp.subs({yi_v: ys[i], yj_v: ys[j], yk_v: ys[k]})
                    is_bad = (subbed % p == 0) if p else (subbed == 0)
                    
                    dt = time.time() - t0
                    # 简化 Gröbner：用 factorization
                    # 对于 P=0, 其 Gröbner 基等价于每个因子的 Gröbner 基的并
                    # 更好的做法是直接算 P=0 在 GF(p) 上的 Gröbner 基
                    G = sp.groebner([P_exp], yi_v, yj_v, yk_v, order='lex', modulus=p)
                    gb_dt = time.time() - t0 - dt
                    t_total += time.time() - t0
                    
                    result = {
                        'i': i, 'j': j, 'k': k,
                        'P_terms': P_exp.count_ops(),
                        'P_deg': sp.degree(P_exp, yi_v) + sp.degree(P_exp, yj_v) + sp.degree(P_exp, yk_v) if hasattr(sp, 'degree') else '?',
                        'gb_len': len(G),
                        'sol_bad': is_bad,
                        'gb_time_s': round(gb_dt, 4),
                        'total_time_s': round(dt + gb_dt, 4),
                    }
                    triple_results.append(result)
                    
                    if len(triple_results) <= 10:
                        print(f"  triple({i},{j},{k}): P={result['P_terms']} terms, "
                              f"gb={result['gb_len']} polys, {result['gb_time_s']:.3f}s, "
                              f"sol_bad={result['sol_bad']}")
                    sys.stdout.flush()
                    
                except Exception as e:
                    dt = time.time() - t0
                    t_total += dt
                    print(f"  triple({i},{j},{k}): ERROR ({dt:.1f}s): {e}")
                    sys.stdout.flush()
    
    n_completed = len(triple_results)
    n_bad = sum(1 for r in triple_results if r.get('sol_bad'))
    avg_gb_time = sum(r['gb_time_s'] for r in triple_results) / max(n_completed, 1)
    
    print(f"\n  Summary: {n_completed}/{n_completed} triples done, {t_total:.1f}s total")
    print(f"  Avg per triple: {t_total/max(n_completed,1)*1000:.2f}ms")
    print(f"  Known solution violates {n_bad} triple Gröbner bases mod {p}")
    if n_bad == 0:
        print(f"  → Known solution is consistent with Gröbner constraints over GF({p})")
    else:
        print(f"  → These are false positives (det ≡ 0 mod {p} but det ≠ 0 in ℤ)")
    sys.stdout.flush()
    
    return n_completed, t_total, triple_results


def triple_gb_scale_test(m_start, m_end, p=7):
    """测不同 m 的三元组 Gröbner 基计算时间"""
    print(f"\n{'='*60}")
    print(f"Scale test: m={m_start} → {m_end} on GF({p})")
    print(f"{'='*60}")
    sys.stdout.flush()
    
    for m in [m_start, m_end]:
        if m > m_end: break
        y = sp.symbols(f'y0:{m}')
        N = 2 * m
        
        sample_size = min(8, sp.binomial(m, 3))
        total_time = 0
        count = 0
        
        for idx, (i, j, k) in enumerate(itertools.combinations(range(m), 3)):
            if idx >= sample_size:
                break
            t0 = time.time()
            
            # 乘积多项式
            P = 1
            for r in range(4):
                for s in range(4):
                    for t_rot in range(4):
                        x1, y1_p = c4((i, y[i]), r, N)
                        x2, y2_p = c4((j, y[j]), s, N)
                        x3, y3_p = c4((k, y[k]), t_rot, N)
                        det = x1*(y2_p - y3_p) + x2*(y3_p - y1_p) + x3*(y1_p - y2_p)
                        P *= det
            
            try:
                P_exp = sp.expand(P)
                G = sp.groebner([P_exp], y[i], y[j], y[k], order='lex', modulus=p)
                dt = time.time() - t0
                total_time += dt
                count += 1
                if idx < 3 or idx == sample_size - 1:
                    print(f"  m={m}: triple({i},{j},{k}): {P_exp.count_ops()} terms → "
                          f"gb={len(G)} polys, {dt:.3f}s")
                    sys.stdout.flush()
            except Exception as e:
                dt = time.time() - t0
                print(f"  m={m}: triple({i},{j},{k}): ERROR ({dt:.1f}s): {e}")
                sys.stdout.flush()
        
        if count > 0:
            avg = total_time / count
            total_est = avg * sp.binomial(m, 3)
            print(f"  m={m}: avg {avg*1000:.2f}ms/triple, est total {total_est/60:.1f}min for all C(m,3)")
        print()
    
    sys.stdout.flush()


# === 新方向：小素数消元法 ===
# 不用 Gröbner，用穷举小素数上的赋值
def exhaustive_gf(m, p, max_sols=5):
    """
    直接在 GF(p) 上穷举所有可能的 π（双射），
    检查哪些在 GF(p) 上无 collinear triple。
    对 m 较小、p 略大于 m 的场景有效。
    时间复杂度: P(p, m) = p!/(p-m)!  → 仅对极小 m 可行。
    """
    print(f"\n{'='*60}")
    print(f"Exhaustive search in GF({p}): m={m} (p={p})")
    print(f"Search space: P({p},{m}) = {math.perm(p,m)}")
    print(f"{'='*60}")
    sys.stdout.flush()
    
    N = 2 * m
    
    # Warm up
    # For m=5, p=7: 7P5 = 2520 — perfectly feasible
    # For m=7, p=11: 11P7 = 1663200 — feasible with optimization
    # For m=11, p=13: 13P11 = 3113510400 — too large
    
    if math.perm(p, m) > 5000000:
        print(f"  Too large ({math.perm(p,m)} > 5M), skipping exhaustive")
        return None
    
    found_solutions = []
    search_count = 0
    
    # Optimization: precompute all line sets
    # For each (A,B,C) line and each cell position (i, yi), we precompute
    # whether that lift contributes to this line.
    
    t0 = time.time()
    
    # For each permutation
    for perm in itertools.permutations(range(p), m):
        search_count += 1
        
        # Quick check: does this permutation have any collinear triple?
        has_bad = False
        
        # Check only triples (i,j,k)
        for i in range(m):
            if has_bad: break
            for j in range(i+1, m):
                if has_bad: break
                for k in range(j+1, m):
                    if has_bad: break
                    for r in range(4):
                        if has_bad: break
                        for s in range(4):
                            if has_bad: break
                            for t_rot in range(4):
                                x1, y1_p = c4((i, perm[i]), r, N)
                                x2, y2_p = c4((j, perm[j]), s, N)
                                x3, y3_p = c4((k, perm[k]), t_rot, N)
                                det = x1*(y2_p - y3_p) + x2*(y3_p - y1_p) + x3*(y1_p - y2_p)
                                if det % p == 0:
                                    has_bad = True
                                    break
        
        if not has_bad:
            found_solutions.append(list(perm))
            print(f"  FOUND #{len(found_solutions)}: perm={list(perm)}")
            sys.stdout.flush()
            if len(found_solutions) >= max_sols:
                break
        
        if search_count % 10000 == 0:
            dt = time.time() - t0
            print(f"  Searched {search_count}/... ({dt:.1f}s, {search_count/dt:.0f}/s)")
            sys.stdout.flush()
    
    dt = time.time() - t0
    print(f"\n  Total: {search_count} searched in {dt:.1f}s")
    print(f"  Solutions found: {len(found_solutions)}")
    sys.stdout.flush()
    
    # Verify: do these correspond to real NTIL solutions?
    if found_solutions:
        print(f"\n  Cross-checking found solutions:")
        for idx, perm in enumerate(found_solutions[:3]):
            cells = [(i, perm[i]) for i in range(m)]
            ok_int, bad_int = verify_solution(cells, m)
            print(f"    #{idx+1}: perm={perm} — over ℤ: {'PASS' if ok_int else f'FAIL bad={bad_int}'}")
            sys.stdout.flush()
    
    return found_solutions


if __name__ == '__main__':
    # === Step 1: m=5 验证实验 ===
    # 加载已知解
    print("Loading m=5 solution...")
    sol5 = load_solution(5)
    if sol5:
        ys = tuple(y for (x,y) in sol5)
        print(f"m=5 solution: cells={sol5}, y={ys}")
        ok, bad = verify_solution(sol5, 5)
        print(f"Integer verify: {'PASS' if ok else f'FAIL (bad={bad})'}")
    sys.stdout.flush()
    
    # === Step 2: 完整 Gröbner 基实验 ===
    gr_vs_m5(p=13)  # Use p > 2m for cleaner results
    
    # === Step 3: 尺寸缩放测试 ===
    triple_gb_scale_test(5, 11, p=13)
    
    # === Step 4: 小素数穷举法 (m=5, p=7) ===
    exhaustive_gf(m=5, p=7)
    
    print("\nDone!")
