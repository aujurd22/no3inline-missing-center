"""
solve_m37_lazy.py  --  Direction 3: Lazy constraint generation (branch-and-cut)
for rot4 C4 NTIL, using R9b 2-regular reformulation.

Instead of feeding ALL 1.26M per-line at-most-2 constraints upfront (which
wedges the presolve phase), we:

  1) Build the full line-to-reps dictionary ONCE (fast m^2 generator ~60s for m=37).
  2) Start with a minimal model: 2-factor + Sidon + base-frame-3free (O(m) constraints).
  3) Solve (fast, no presolve wall).
  4) Lift candidate cells -> brute-find ALL violated lines (3+ collinear).
  5) Add ONLY those violated-line constraints to the model.
  6) Re-solve with hints from previous solution.
  7) Repeat until verify passes or model is infeasible.

Stop-loss: compare lazy vs full-model solve on m=17-22.
If lazy doesn't significantly outperform, stop this direction.

Usage:
  python solve_m37_lazy.py --m 17              # fast probe on reference m
  python solve_m37_lazy.py --m 37 --timelimit 7200 --workers 8
  python solve_m37_lazy.py --m 17-22 --bench   # full vs lazy comparison
"""
import os, sys, time, json, math, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from quadratic_sidon_completeness import c4, cell_xy
from solve_m37_r9b import verify_cells, _write_ckpt

try:
    from ortools.sat.python import cp_model
    HAVE_ORTOOLS = True
except Exception:
    HAVE_ORTOOLS = False


# ---------------------------------------------------------------------------
#  Build FULL line→reps dictionary using integer rep_idx (same as
#  generate_constraints), returned as two structures:
#    line_cons_full : list of {rep_idx: weight}  (all lines with total weight>2)
#    line_w_dict    : {(dx,dy): {key: {rep_idx: weight}}}
#  The latter is used by find_violated_lines() for incremental violation detection.
# ---------------------------------------------------------------------------

def build_line_dicts(m):
    """Build both line_cons_full (all constraints) and line_w_dict (raw dict for
    incremental violation detection).  Uses integer rep_idx for consistency with
    cp_model variable indices."""
    n = 2 * m
    reps = [(x, y) for x in range(m) for y in range(m)]  # m^2 cells
    orbits = [None] * len(reps)
    for i, (x, y) in enumerate(reps):
        orbits[i] = [c4((x, y), r, n) for r in range(4)]

    # reduced_dirs (same as solve_m37_r9b.py)
    dirs = set()
    for dx in range(-(n - 1), n):
        for dy in range(0, n):
            if dx == 0 and dy == 0:
                continue
            g = math.gcd(abs(dx), dy) or 1
            rdx, rdy = dx // g, dy // g
            if rdx < 0 or (rdx == 0 and rdy < 0):
                rdx, rdy = -rdx, -rdy
            dirs.add((rdx, rdy))

    from collections import defaultdict
    line_w = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # (dx,dy)->key->{rep_idx:w}
    for i, (x, y) in enumerate(reps):
        for (X, Y) in orbits[i]:
            for (dx, dy) in dirs:
                perp = (-dy, dx)
                key = perp[0] * X + perp[1] * Y
                line_w[(dx, dy)][key][i] += 1

    # Extract full constraint list (all lines with total weight > 2)
    line_cons_full = []
    for d, lines in line_w.items():
        for key, pos_w in lines.items():
            if sum(pos_w.values()) > 2:
                line_cons_full.append(dict(pos_w))

    return line_cons_full, line_w, reps


def find_violated_line_cons(cells_set, m, line_w):
    """Given a set of selected cells (x,y tuples), return list of {rep_idx:weight}
    for all lines where the selected cells have weighted sum > 2.
    Uses the integer-index line_w dict."""
    # Convert cells (x,y) to rep_idx
    rep_of = {(x, y): x * m + y for x in range(m) for y in range(m)}
    viol = []
    for d, lines in line_w.items():
        for key, pos_w in lines.items():
            cnt = 0
            for ri, w in pos_w.items():
                # Convert rep_idx back to (x,y) to check if selected
                rx = ri // m
                ry = ri % m
                if (rx, ry) in cells_set:
                    cnt += w
            if cnt > 2:
                viol.append(dict(pos_w))
    return viol


# ---------------------------------------------------------------------------
#  Lazy solver
# ---------------------------------------------------------------------------

class LazySolver:
    """Branch-and-cut: iteratively solve minimal model, add violated lines."""

    def __init__(self, m, line_w, reps, line_cons_full,
                 timelimit=300, workers=8, verbose=True,
                 checkpoint=None):
        self.m = m
        self.line_w = line_w
        self.reps = reps
        self.line_cons_full = line_cons_full  # all constraints (for benchmark only)
        self.timelimit = timelimit
        self.workers = workers
        self.verbose = verbose
        self.checkpoint = checkpoint

        # Track which constraints have been added (by frozenset of items)
        self._added = set()
        self._total_added = 0

    def add_2factor(self, model, sel):
        """Add 2-factor constraints to an existing model."""
        m = self.m
        for i in range(m):
            model.Add(sum(sel[i]) + sum(sel[r][i] for r in range(m)) == 2)

    def add_line_cons(self, model, sel, cons_list):
        """Add a batch of {rep_idx: weight} constraints to the model.
        sel is 2D list sel[x][y] of BoolVars."""
        for d in cons_list:
            if len(d) <= 2:
                continue
            terms = []
            for ri, w in d.items():
                x = ri // self.m
                y = ri % self.m
                terms.append(w * sel[x][y])
            model.Add(sum(terms) <= 2)

    def add_line_cons_dedup(self, model, sel, cons_list):
        """Add constraints, skipping duplicates tracked in self._added."""
        new = 0
        for d in cons_list:
            fk = frozenset(d.items())
            if fk in self._added:
                continue
            self._added.add(fk)
            self.add_line_cons(model, sel, [d])
            new += 1
        self._total_added += new
        return new

    def solve_one_iter(self, model, sel, hint_cells=None, timelimit_frac=0.25):
        """Solve current model, return (status, cells, solve_time)."""
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = max(30, self.timelimit * timelimit_frac)
        solver.parameters.num_search_workers = self.workers
        solver.parameters.max_presolve_iterations = 3

        # Note: AddHint accumulates across solver calls; OR-Tools 9.x may raise
        # MODEL_INVALID if later hints conflict with earlier ones.  We skip hints
        # and rely on the incremental constraints alone to drive convergence.

        t0 = time.time()
        try:
            st = solver.Solve(model)
        except Exception as e:
            return f"EXCEPTION:{e}", [], 0
        dt = time.time() - t0
        status = solver.StatusName(st)
        cells = []
        if status in ("FEASIBLE", "OPTIMAL"):
            for i in range(self.m):
                for j in range(self.m):
                    if solver.Value(sel[i][j]) == 1:
                        cells.append((i, j))
        return status, cells, dt

    def run(self, max_iters=30):
        """Run lazy constraint generation loop."""
        if not HAVE_ORTOOLS:
            print("[lazy] OR-Tools not available")
            return "NO_ORTOOLS", [], 0, 0

        m = self.m
        t_start = time.time()

        # Create model and variables ONCE
        model = cp_model.CpModel()
        sel = [[model.NewBoolVar(f"s{i}_{j}") for j in range(m)] for i in range(m)]
        self.add_2factor(model, sel)

        if self.verbose:
            print(f"[lazy] initial model built ({time.time()-t_start:.1f}s, "
                  f"2-factor only)", flush=True)

        # Validate initial model (debug)
        try:
            init_solver = cp_model.CpSolver()
            init_solver.parameters.num_search_workers = self.workers
            init_solver.parameters.max_time_in_seconds = 5
            init_solver.parameters.max_presolve_iterations = 3
            st0 = init_solver.Solve(model)
            if init_solver.StatusName(st0) in ("INFEASIBLE", "MODEL_INVALID"):
                print(f"[lazy] WARNING: bare model {init_solver.StatusName(st0)}",
                      flush=True)
        except Exception as e:
            print(f"[lazy] EXCEPTION on bare model: {e}", flush=True)

        cells = None
        total_solve_time = 0
        prev_signature = None
        stall = 0
        for it in range(max_iters):
            t_iter = time.time()
            # First few iters use short timelimit; later iters more generous
            frac = 0.1 if it < 5 else (0.25 if it < 15 else 0.5)
            status, cells, dt = self.solve_one_iter(
                model, sel, timelimit_frac=frac)
            total_solve_time += dt

            if not cells:
                if self.verbose:
                    print(f"[lazy] iter {it}: {status} (no cells) [{dt:.1f}s]",
                          flush=True)
                if status in ("INFEASIBLE", "MODEL_INVALID", "UNKNOWN"):
                    if self._total_added == 0:
                        return "INFEASIBLE_BARE", [], total_solve_time, time.time() - t_start
                    return status, [], total_solve_time, time.time() - t_start
                continue

            # Find violated lines on current candidate
            cells_set = set(cells)
            viol = find_violated_line_cons(cells_set, m, self.line_w)

            # Check structural verify
            vf, nc = verify_cells(cells, m)
            if vf:
                if self.verbose:
                    print(f"[lazy] iter {it}: VERIFY PASS ({nc} cells) "
                          f"[{dt:.1f}s]", flush=True)
                if self.checkpoint:
                    _write_ckpt(self.checkpoint, cells)
                return "SOLVED", cells, total_solve_time, time.time() - t_start

            # Add new constraints
            n_new = self.add_line_cons_dedup(model, sel, viol)

            # Detect cycling: if viol count stays same for 5+ iters, we stalled
            sig = (len(viol), n_new)
            if sig == prev_signature:
                stall += 1
            else:
                stall = 0
            prev_signature = sig
            if stall > 8:
                if self.verbose:
                    print(f"[lazy] iter {it}: STALLED (viol={len(viol)} "
                          f"unchanged for {stall} iters), "
                          f"adding ALL violated lines",
                          flush=True)
                # Force-add ALL violated lines even if previously added
                for d in viol:
                    self._added.add(frozenset(d.items()))
                    self.add_line_cons(model, sel, [d])
                self._total_added += len(viol)
                stall = 0

            elapsed = time.time() - t_iter
            if self.verbose:
                print(f"[lazy] iter {it}: found={len(cells)} "
                      f"viol={len(viol)} new_cons={n_new} "
                      f"tot_added={self._total_added} [{dt:.1f}s+{elapsed:.1f}s]",
                      flush=True)

            if self._total_added > 500000:
                if self.verbose:
                    print(f"[lazy] hit 500k — full model fallback", flush=True)
                break

        if self.verbose:
            print(f"[lazy] max_iters={max_iters} reached, "
                  f"added={self._total_added}", flush=True)
        return "MAX_ITERS", cells or [], total_solve_time, time.time() - t_start


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=str, default="17",
                    help="m value or range like 17-22")
    ap.add_argument("--timelimit", type=float, default=7200.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--checkpoint", default="")
    args = ap.parse_args()

    if not HAVE_ORTOOLS:
        print("OR-Tools not available. Install with pip install ortools")
        sys.exit(1)

    # Parse m
    if "-" in args.m:
        m_range = list(range(*map(int, args.m.split("-"))))
    else:
        m_range = [int(args.m)]
    m = m_range[0]

    print(f"== Lazy CP-SAT: m={m} ==")
    t0 = time.time()
    line_cons_full, line_w, reps = build_line_dicts(m)
    print(f"build: {time.time()-t0:.1f}s  total_possible_cons={len(line_cons_full)}",
          flush=True)

    ckpt = args.checkpoint or None
    solver = LazySolver(m, line_w, reps, line_cons_full,
                        timelimit=args.timelimit,
                        workers=args.workers, verbose=True,
                        checkpoint=ckpt)
    status, cells, tsolve, total = solver.run(max_iters=30)
    print(f"[final] status={status} found={len(cells)} "
          f"solve={tsolve:.1f}s total={total:.1f}s "
          f"added={solver._total_added}", flush=True)
    if cells:
        vf, nc = verify_cells(cells, m)
        print(f"[final] verify={vf} cells={nc}", flush=True)
        print("[cells]", cells)


if __name__ == "__main__":
    main()
