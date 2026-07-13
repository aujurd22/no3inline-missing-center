import json, random, sys
# Perturb a known rot4-NTIL solution by K random swaps of ys (keeps 2-factor).
# Usage: perturb.py <sol.json> <K> <out.json> <seed>
f, K, out, seed = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
d = json.load(open(f))
cells = [[x, y] for x, y in d['cells']]
random.seed(seed)
m = len(cells)
for _ in range(K):
    i, j = random.sample(range(m), 2)
    cells[i][1], cells[j][1] = cells[j][1], cells[i][1]
json.dump({"m": m, "n": 2 * m, "cells": cells}, open(out, 'w'))
print("perturbed K=%d -> %s (m=%d)" % (K, out, m))
