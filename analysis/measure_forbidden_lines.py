"""
measure_forbidden_lines.py — For known valid rot4 solutions,
measure how many distinct forbidden (α,β) positions each pair
actually faces from the 16 quadratic forms.

This determines the EFFECTIVE constraint density vs the RAW line count.
"""

import os, sys, glob
from itertools import combinations, product
from collections import defaultdict
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

def det_value(cells, rot_indices, N):
    pts = [rotation_coords(cx, cy, r, N) for (cx,cy), r in zip(cells, rot_indices)]
    (x1,y1),(x2,y2),(x3,y3) = pts
    return (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

# 16 canonical rotation triples (one per equivalence class)
CANONICAL_ROTS = [
    (0,0,0), (0,0,1), (0,0,2), (0,0,3),
    (0,1,0), (0,1,1), (0,1,2), (0,1,3),
    (0,2,0), (0,2,1), (0,2,2), (0,2,3),
    (0,3,0), (0,3,1), (0,3,2), (0,3,3),
]

def count_forbidden_lines_per_pair(pairs, m, pair_idx):
    """
    For a given pair (pair_idx), count how many distinct forbidden
    (α,β) positions exist based on all C(m-1,2) × 16 form evaluations.
    A position (α,β) is "forbidden" if some (j,k) and rotation triple
    makes det = 0 when pair_idx takes that position.
    
    Actually, we need to be smarter. For each (j,k) pair-of-pairs and
    each rotation triple, the condition det=0 defines a LINE in (α,β).
    But I'll count the actual forbidden POINTS on the m×m grid.
    """
    n = 2 * m
    
    # Known values of all pairs except pair_idx
    other_pairs = [p for i, p in enumerate(pairs) if i != pair_idx]
    other_cells = [(m - (a+1)//2, m - (b+1)//2) for a,b in other_pairs]
    
    # Get actual (a,b) of the target pair
    a0, b0 = pairs[pair_idx]
    alpha0 = (a0 + 1) // 2
    beta0 = (b0 + 1) // 2
    
    # For each combination of 2 other pairs, check all 16 rotation triples
    # When does det = 0 assuming pair_idx takes position (α,β)?
    forbidden_positions = set()
    line_count = 0
    
    for j, k in combinations(range(len(other_pairs)), 2):
        cj = other_cells[j]  # (x_j, y_j)
        ck = other_cells[k]  # (x_k, y_k)
        
        for rot in CANONICAL_ROTS:
            r1, r2, r3 = rot
            
            # The determinant for triple (target, j, k) with rotations (r1,r2,r3)
            # We need det = 0 as a function of (α,β) for target
            # 
            # For rotation r1 applied to target cell at (x,y) = (m-α, m-β):
            #   if r1=0: P1 = (m-α, m-β)
            #   if r1=1: P1 = (N-1-(m-β), m-α) = (m-1+β, m-α)
            #   if r1=2: P1 = (N-1-(m-α), N-1-(m-β)) = (m-1+α, m-1+β)
            #   if r1=3: P1 = (m-β, N-1-(m-α)) = (m-β, m-1+α)
            #
            # where α,β ∈ {1,...,m}
            
            # Same for P2, P3 but with fixed (j,k) cells
            
            # Actually let me just brute-force all m² positions
            for alpha in range(1, m+1):
                for beta in range(1, m+1):
                    # Skip the actual position (it must NOT be forbidden)
                    if alpha == alpha0 and beta == beta0:
                        continue
                    
                    # Build the 3 cells with target at (alpha,beta)
                    cells = [(m-alpha, m-beta), cj, ck]
                    
                    det = det_value(cells, [r1, r2, r3], n)
                    if det == 0:
                        forbidden_positions.add((alpha, beta))
        
        line_count += len(CANONICAL_ROTS)
    
    return forbidden_positions, line_count

def main():
    print("=" * 100)
    print("EFFECTIVE FORBIDDEN POSITIONS PER PAIR (known valid solutions)")
    print("=" * 100)
    
    for m_target in [5, 8, 10, 12, 16, 20, 28, 36]:
        files = sorted(glob.glob(os.path.join(CACHE, f'n{2*m_target}_rot4*')))
        if not files:
            continue
        path = files[0]
        with open(path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        
        for line in lines[:3]:  # Check first 3 solutions
            pts = rot4.decode_line(line, 2*m_target)
            if len(pts) != 4*m_target or not rot4.is_valid(pts, 2*m_target):
                continue
            pairs = [(2*(m_target-x)-1, 2*(m_target-y)-1) for (x,y) in pts 
                     if x < m_target and y < m_target]
            if len(pairs) != m_target:
                continue
            
            # For pairs 0, m_target//4, m_target//2 (spread across the solution)
            sample_indices = [0, m_target//4, m_target//2, m_target-1] if m_target > 3 else [0, m_target-1]
            
            for pair_idx in sample_indices:
                if pair_idx >= m_target:
                    continue
                forbidden, total_lines = count_forbidden_lines_per_pair(pairs, m_target, pair_idx)
                
                grid_positions = m_target * m_target
                raw_ratio = total_lines / grid_positions
                effective_ratio = len(forbidden) / grid_positions
                
                print(f"  m={m_target}, pair#{pair_idx}: "
                      f"{len(forbidden)}/{grid_positions} forbidden "
                      f"(raw: {total_lines} lines, {raw_ratio:.1f}x grid)")
            
            break  # Just first solution per m
        
        print()

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)
    print("""
If effective forbidden count << raw line count:
  → Lines are highly degenerate (many pass through same points)
  → The quadratic constraints have a small algebraic basis
  → The "real" constraint density is much lower than raw count

If effective forbidden count ≈ grid size:
  → Each pair is nearly forced to a specific position
  → The solution is "tight" with little flexibility
  → m=37 may be impossible due to over-constraint
""")


if __name__ == '__main__':
    main()
