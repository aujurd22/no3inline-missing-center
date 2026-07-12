"""
verify_quadratic_distance.py — Check whether known valid solutions
"barely avoid" or "comfortably avoid" the quadratic collinearity conditions.

For each valid solution, compute ALL 3-pair × (0,1,2) determinants.
If most are large (>>0), the solution "comfortably" avoids collinearity.
If many are close to zero, it "barely" avoids collinearity.
"""

import os, sys, glob
from collections import Counter
from itertools import combinations, product
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rot4_loader as rot4

CACHE = rot4.CACHE

def rotation_coords(cx, cy, r, N):
    if r == 0:
        return (cx, cy)
    elif r == 1:
        return (N-1-cy, cx)
    elif r == 2:
        return (N-1-cx, N-1-cy)
    elif r == 3:
        return (cy, N-1-cx)

def collinear_det(cells, rot_indices, N):
    pts = [rotation_coords(cx, cy, r, N) for (cx,cy), r in zip(cells, rot_indices)]
    (x1,y1),(x2,y2),(x3,y3) = pts
    return (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

# Key rotation patterns that produce quadratic conditions
KEY_ROTS = [(0,1,2), (0,1,3), (0,2,3), (1,2,3), (0,0,1), (0,0,3)]

def analyze_solution(pairs, m):
    n = 2 * m
    cells = [(m - (a+1)//2, m - (b+1)//2) for a,b in pairs]
    results = {}
    
    for rot in KEY_ROTS:
        dets = []
        for i, j, k in combinations(range(m), 3):
            cell_triple = [cells[i], cells[j], cells[k]]
            det = collinear_det(cell_triple, rot, n)
            dets.append(det)
        results[rot] = {
            'min': min(dets),
            'max': max(dets),
            'mean': sum(dets) / len(dets),
            'near_zero': sum(1 for d in dets if abs(d) <= 1),  # "barely avoided"
            'total': len(dets),
        }
    return results

def main():
    print("=" * 90)
    print("QUADRATIC CONDITION DISTANCE ANALYSIS")
    print("For each valid solution: how 'near' is it to collinearity?")
    print("=" * 90)
    
    for m_target in [4, 5, 8, 10, 12, 16, 20, 28, 36]:
        files = sorted(glob.glob(os.path.join(CACHE, f'n{2*m_target}_rot4*')))
        if not files:
            continue
        path = files[0]
        with open(path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        
        valid_pairs_list = []
        for line in lines[:10]:
            pts = rot4.decode_line(line, 2*m_target)
            if len(pts) != 4*m_target or not rot4.is_valid(pts, 2*m_target):
                continue
            pairs = [(2*(m_target-x)-1, 2*(m_target-y)-1) for (x,y) in pts 
                     if x < m_target and y < m_target]
            if len(pairs) == m_target:
                valid_pairs_list.append(pairs)
        
        if not valid_pairs_list:
            continue
        
        aggregated = {}
        for rot in KEY_ROTS:
            all_dets = []
            near_zero_counts = []
            for pairs in valid_pairs_list[:5]:
                r = analyze_solution(pairs, m_target)[rot]
                all_dets.extend([0])  # placeholder
                near_zero_counts.append(r['near_zero'])
            
            # Better: compute aggregate across all sols
            all_near_zero = 0
            all_total = 0
            for pairs in valid_pairs_list[:5]:
                r = analyze_solution(pairs, m_target)[rot]
                all_near_zero += r['near_zero']
                all_total += r['total']
            
            aggregated[rot] = {
                'near_zero': all_near_zero,
                'total': all_total,
                'pct': 100.0 * all_near_zero / all_total if all_total else 0,
            }
        
        print(f"\nm={m_target}: {len(valid_pairs_list)} solutions")
        for rot in KEY_ROTS:
            a = aggregated[rot]
            print(f"  r={rot}: {a['near_zero']}/{a['total']} near-zero = {a['pct']:.2f}%")

    print("\n\nAnalysis complete.")
    print("'near-zero' = |det| ≤ 1 means the solution is 'barely avoiding' collinearity")
    print("If near-zero % is high, valid solutions critically depend on tiny margin.")

if __name__ == '__main__':
    main()
