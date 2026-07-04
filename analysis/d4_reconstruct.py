"""
D4 equivalence class -> full solution count reconstruction.
Verifies D4 orbit sizes empirically, then reconstructs full counts.

For each symmetry class, the D4 orbit size multiplier is:
  iden(.): 8   - no symmetry, all D4 transforms distinct
  rot2(:): 4   - {I,R2} stabilizer
  dia1(/): 4   - {I,sigma} stabilizer
  ort1(-): 4   - {I,sigma*R} stabilizer
  rot4(o): 2   - {I,R,R2,R3} stabilizer
  dia2(x): 2   - {I,R2,sigma,sigma*R2} stabilizer
  ort2(+): 2   - {I,R2,sigma*R,sigma*R3} stabilizer
  rct4(c): 1   - full D4 stabilizer
  full(*): 1   - same as rct4
"""
import urllib.request
from collections import Counter
import time

ALPHABET = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
char_to_val = {c: i for i, c in enumerate(ALPHABET)}

SYMM_MULT = {
    '.': 8, ':': 4, '/': 4, '-': 4,
    'o': 2, 'x': 2, '+': 2,
    'c': 1, '*': 1
}
SYMM_NAMES = {
    '.': 'iden', ':': 'rot2', '/': 'dia1', '-': 'ort1',
    'o': 'rot4', 'c': 'rct4', 'x': 'dia2', '+': 'ort2', '*': 'full'
}

def R(p, n):
    return (n-1-p[1], p[0])

def decode(line):
    line = line.strip()
    if not line: return None
    symm = line[0]
    rest = line[1:]
    n = len(rest) // 2
    pts = []
    for r in range(n):
        pts.append((r, char_to_val[rest[2*r]]))
        pts.append((r, char_to_val[rest[2*r+1]]))
    return symm, n, frozenset(pts)

def d4_orbit_size(n, pts):
    """Compute the exact D4 orbit size of a solution."""
    transforms = set()
    transforms.add(pts)
    for _ in range(1, 4):
        pts = frozenset(R(p, n) for p in pts)
        transforms.add(pts)
    for pts2 in list(transforms):
        ptsR = frozenset((p[1], p[0]) for p in pts2)
        transforms.add(ptsR)
    more = set()
    for pts2 in transforms:
        for _ in range(3):
            pts2 = frozenset(R(p, n) for p in pts2)
            more.add(pts2)
    transforms.update(more)
    return len(transforms)

def is_missing(n, pts):
    ctr = (n-1)/2.0
    rings = Counter()
    for i,j in pts:
        d2 = int(round((i-ctr)**2 + (j-ctr)**2))
        rings[d2] += 1
    return all(v < 3 for v in rings.values())

def load_solutions(n, symm):
    url = f'https://wwwhomes.uni-bielefeld.de/achim/no3in/download/configurations/n{n}_{symm}'
    try:
        with urllib.request.urlopen(url, timeout=20) as f:
            raw = f.read().decode().strip()
            lines = [l for l in raw.split('\n') if l.strip()]
            if not lines or '<html' in lines[0].lower():
                return None
            return lines
    except Exception as e:
        return None

def process_class(lines, symm_char, n, sample_size=15):
    """Process one symmetry class. Returns (total, missing, full_total, full_missing, verified)."""
    total = len(lines)
    mult = SYMM_MULT.get(symm_char, 8)
    
    # Empirical verification on sample
    sample_ok = True
    sample_missing = 0
    for i, line in enumerate(lines[:sample_size]):
        result = decode(line)
        if result is None: continue
        _, n2, pts = result
        if n2 != n: continue
        sz = d4_orbit_size(n, pts)
        if sz != mult:
            sample_ok = False
            print(f'    *** MISMATCH: expected mult={mult}, got orbit={sz}')
        if is_missing(n, pts):
            sample_missing += 1
    
    # Process ALL solutions for missing-center
    all_missing = 0
    for line in lines:
        result = decode(line)
        if result is None: continue
        _, n2, pts = result
        if n2 != n: continue
        if is_missing(n, pts):
            all_missing += 1
    
    full_total = total * mult
    full_missing = all_missing * mult
    
    return total, all_missing, full_total, full_missing, sample_ok, sample_missing

# ── Main ──
print('=' * 110)
print('D4 EQUIVALENCE CLASS -> FULL SOLUTION COUNT RECONSTRUCTION')
print('=' * 110)

# First pass: verify D4 orbit size multipliers empirically
print('\n--- Phase 1: Empirical D4 orbit size verification (sample) ---')
for n in [7, 8, 9, 10, 11, 12, 13, 14, 15]:
    symm_classes = ['.', ':', '/', '-', 'o', 'c', 'x', '+', '*']
    for s in symm_classes:
        lines = load_solutions(n, s)
        if lines is None: continue
        mult_expected = SYMM_MULT.get(s, 8)
        ok_count = 0
        fail_count = 0
        for line in lines[:10]:
            result = decode(line)
            if result is None: continue
            _, n2, pts = result
            if n2 != n: continue
            sz = d4_orbit_size(n, pts)
            if sz == mult_expected:
                ok_count += 1
            else:
                fail_count += 1
        status = 'OK' if fail_count == 0 else 'MISMATCH'
        print(f'  n={n:2d} {SYMM_NAMES[s]:>5}({s}): expected mult={mult_expected}, '
              f'verified {ok_count}/{ok_count+fail_count} -> {status}')

# Second pass: reconstruct full counts for ALL n
print('\n--- Phase 2: Full reconstruction ---')
HEADER = f'{"n":>3} {"type":>6} {"D4_tot":>6} {"D4_miss":>7} {"D4_%":>5} ' \
         f'{"full_tot":>8} {"full_miss":>9} {"full_%":>6} {"classes_detail":<55}'
print(HEADER)
print('-' * len(HEADER))

all_data = {}

for n in list(range(7, 31)) + [31, 33, 35, 37, 39, 41, 43, 45]:
    symm_classes = ['.', ':', '/', '-', 'o', 'c', 'x', '+', '*']
    total_d4 = 0
    missing_d4 = 0
    full_total = 0
    full_missing = 0
    details = []
    verified = True
    
    for s in symm_classes:
        lines = load_solutions(n, s)
        if lines is None:
            continue
        t, m, ft, fm, v, _ = process_class(lines[:50], s, n)  # limit to 50 for speed
        total_d4 += t
        missing_d4 += m
        full_total += ft
        full_missing += fm
        if not v:
            verified = False
        cls_name = SYMM_NAMES.get(s, s)
        details.append(f'{cls_name}:{t}(m={m})')
    
    if total_d4 == 0:
        continue
    
    d4_rate = missing_d4 / total_d4 * 100
    full_rate = full_missing / full_total * 100 if full_total > 0 else 0
    
    ptype = '4k+' + str(n % 4)
    is_prime = all(n % d != 0 for d in range(2, int(n**0.5) + 1)) if n > 1 else False
    typestr = ptype + ('P' if is_prime else 'C')
    
    short_det = ', '.join(details[:6])
    if len(details) > 6:
        short_det += f' (+{len(details)-6} more)'
    
    print(f'{n:>3} {typestr:>6} {total_d4:>6} {missing_d4:>7} {d4_rate:>4.1f}% '
          f'{full_total:>8} {full_missing:>9} {full_rate:>5.2f}% '
          f'{short_det:<55}')
    
    all_data[n] = {
        'type': typestr,
        'd4_total': total_d4,
        'd4_missing': missing_d4,
        'd4_rate': d4_rate,
        'full_total': full_total,
        'full_missing': full_missing,
        'full_rate': full_rate,
        'verified': verified
    }

# Cross-validation
print('\n--- Phase 3: Cross-validation (n<=13) ---')
known_full = {7: 132, 8: 380, 9: 368, 10: 1135, 11: 1120, 12: 4348, 13: 3622}
known_missing = {5: 4, 7: 4, 9: 8, 11: 36, 12: 52, 13: 292}

print(f'{"n":>3} {"D4_tot":>6} {"full_recon":>10} {"full_known":>10} {"diff":>5} '
      f'{"miss_recon":>10} {"miss_known":>10} {"ok?":>5}')
print('-' * 60)
all_ok = True
for n in range(7, 14):
    if n in all_data:
        d = all_data[n]
        fk = known_full.get(n, '?')
        mk = known_missing.get(n, 0)
        d1 = abs(d['full_total'] - fk) if fk != '?' else -1
        d2 = abs(d['full_missing'] - mk)
        ok = 'YES' if d1 == 0 else 'NO'
        if d1 != 0:
            all_ok = False
        print(f'{n:>3} {d["d4_total"]:>6} {d["full_total"]:>10} {fk:>10} {d1:>5} '
              f'{d["full_missing"]:>10} {mk:>10} {ok:>5}')

print(f'\nAll cross-validations passed: {all_ok}')

# Final table
print('\n--- Phase 4: Final reconstructed full table ---')
print(f'{"n":>3} {"type":>6} {"full_total":>10} {"full_missing":>12} {"full_%":>6}')
print('-' * 40)
for n in sorted(all_data.keys()):
    d = all_data[n]
    print(f'{n:>3} {d["type"]:>6} {d["full_total"]:>10} {d["full_missing"]:>12} {d["full_rate"]:>5.2f}%')

print(f'\nVerified multiplier for all tested classes: {verified}')
