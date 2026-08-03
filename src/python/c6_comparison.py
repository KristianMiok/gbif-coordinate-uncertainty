"""
C6 -- the comparison. Four estimators on the same virtual species.

The design point is HETEROGENEOUS radii. Real occurrence data is a mixture:
some records locate to metres, some to kilometres. Every method handles the
mixture differently, and that is what separates them.

  threshold : keep only precise records, discard the rest (standard practice;
              trades sample size for apparent precision)
  naive     : keep everything at face value, ignore uncertainty
  NEP       : Smith et al. 2023 GEB -- give each imprecise record the
              environment within its region CLOSEST TO THE MEAN of the precise
              records. Deliberately conservative, hence biased inward.
  COS       : Hefley et al. 2017 MEE -- integrate the intensity over each
              record's own support. Unbiased under overlap (established as C1).

beta0 is profiled out in closed form throughout, so every estimator is compared
on the same footing:

    ell(b1) = n log n - n log T - n + sum_i log g_i
    T   = sum_cells exp(b1 x)
    g_i = the estimator's contribution for record i

Reported: bias in b1, and RANGE INFLATION -- predicted over true suitable area,
the quantity that actually reaches conservation decisions.
"""
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize_scalar
import time

N = 200                    # raster side, cells
B1_TRUE = 1.0
N_EXPECT = 500
R_COARSE = 8.0             # imprecise records, in cells
SA_RANGES = [2.0, 8.0, 32.0]
FRAC_COARSE = [0.25, 0.50, 0.75]
N_REP = 30
SEED = 20260803

rng = np.random.default_rng(SEED)


def field(rng, scale):
    z = gaussian_filter(rng.normal(size=(N, N)), scale, mode="wrap")
    return (z - z.mean()) / z.std()


def offsets(r):
    R = int(np.ceil(r))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    m = (xx ** 2 + yy ** 2) <= r ** 2
    return yy[m], xx[m]


DY, DX = offsets(R_COARSE)
NDISC = len(DY)


def simulate(x, frac_coarse, rng):
    lam = np.exp(B1_TRUE * x)
    lam = lam / lam.sum() * N_EXPECT
    n = rng.poisson(lam.sum())
    idx = rng.choice(lam.size, size=n, replace=True, p=(lam / lam.sum()).ravel())
    tr, tc = np.unravel_index(idx, (N, N))
    coarse = rng.random(n) < frac_coarse
    rr, cc = tr.copy(), tc.copy()
    k = int(coarse.sum())
    if k:
        pick = rng.integers(0, NDISC, k)
        rr[coarse] = (tr[coarse] + DY[pick]) % N
        cc[coarse] = (tc[coarse] + DX[pick]) % N
    return rr, cc, coarse


def _fit(loglik_terms, x):
    """Maximise the beta0-profiled IPP likelihood over b1."""
    def nll(b1):
        e = np.exp(b1 * x)
        g = loglik_terms(b1, e)
        if np.any(g <= 0):
            return np.inf
        n = len(g)
        return -(n * np.log(n) - n * np.log(e.sum()) - n + np.log(g).sum())
    return minimize_scalar(nll, bounds=(-2.0, 4.0), method="bounded",
                           options=dict(xatol=1e-4)).x


def fit_points(x, xs):
    """Records treated as points with covariate values xs."""
    return _fit(lambda b1, e: np.exp(b1 * xs), x)


def fit_cos(x, rows, cols, coarse):
    """Integrate over each record's own support; precise records keep one cell."""
    if coarse.any():
        Rd = (rows[coarse][:, None] + DY[None, :]) % N
        Cd = (cols[coarse][:, None] + DX[None, :]) % N
    pr, pc = rows[~coarse], cols[~coarse]

    def terms(b1, e):
        out = []
        if (~coarse).any():
            out.append(e[pr, pc])
        if coarse.any():
            out.append(e[Rd, Cd].mean(axis=1))   # mean, not sum: unit-mass kernel
        return np.concatenate(out)
    return _fit(terms, x)


def nep_environments(x, rows, cols, coarse):
    """Smith et al.: an imprecise record takes the environment in its region
    closest to the MEAN environment of the precise records."""
    xs = np.empty(len(rows))
    xs[~coarse] = x[rows[~coarse], cols[~coarse]]
    target = xs[~coarse].mean() if (~coarse).any() else 0.0
    if coarse.any():
        R = (rows[coarse][:, None] + DY[None, :]) % N
        C = (cols[coarse][:, None] + DX[None, :]) % N
        cand = x[R, C]
        j = np.abs(cand - target).argmin(axis=1)
        xs[coarse] = cand[np.arange(len(j)), j]
    return xs


def range_inflation(x, b1_hat):
    """Predicted over true suitable area, both thresholded at half of max."""
    t = np.exp(B1_TRUE * x); t /= t.max()
    p = np.exp(b1_hat * x);  p /= p.max()
    return (p > 0.5).mean() / max((t > 0.5).mean(), 1e-9)


METHODS = ("threshold", "naive", "NEP (Smith)", "COS")

print(f"raster {N}x{N} | b1_true {B1_TRUE:.2f} | coarse radius {R_COARSE} cells")
print(f"{N_REP} replicates per cell | expected n = {N_EXPECT}\n")
print(f"{'SA':>4}{'coarse':>8} | " +
      " | ".join(f"{m:^20}" for m in METHODS))
print(f"{'':>4}{'':>8} | " +
      " | ".join(f"{'b1':>5} {'infl':>5} {'n':>6}" for _ in METHODS))
print("-" * 100)

t0 = time.perf_counter()
summary = {}

for sa in SA_RANGES:
    for fc in FRAC_COARSE:
        acc = {k: {"b": [], "infl": [], "n": []} for k in METHODS}
        for _ in range(N_REP):
            x = field(rng, sa)
            rows, cols, coarse = simulate(x, fc, rng)
            if len(rows) < 20 or (~coarse).sum() < 5:
                continue

            b = fit_points(x, x[rows[~coarse], cols[~coarse]])
            acc["threshold"]["b"].append(b)
            acc["threshold"]["n"].append(int((~coarse).sum()))
            acc["threshold"]["infl"].append(range_inflation(x, b))

            b = fit_points(x, x[rows, cols])
            acc["naive"]["b"].append(b); acc["naive"]["n"].append(len(rows))
            acc["naive"]["infl"].append(range_inflation(x, b))

            b = fit_points(x, nep_environments(x, rows, cols, coarse))
            acc["NEP (Smith)"]["b"].append(b)
            acc["NEP (Smith)"]["n"].append(len(rows))
            acc["NEP (Smith)"]["infl"].append(range_inflation(x, b))

            b = fit_cos(x, rows, cols, coarse)
            acc["COS"]["b"].append(b); acc["COS"]["n"].append(len(rows))
            acc["COS"]["infl"].append(range_inflation(x, b))

        summary[(sa, fc)] = acc
        line = f"{sa:>4.0f}{fc:>8.0%} | "
        line += " | ".join(
            f"{np.mean(acc[k]['b']):>5.2f} {np.mean(acc[k]['infl']):>5.2f} "
            f"{np.mean(acc[k]['n']):>6.0f}" for k in METHODS)
        print(line)

print(f"\nelapsed {time.perf_counter() - t0:.0f}s")
print("b1_true = 1.00; range inflation 1.00 = unbiased area.\n")

# ------------------------------------------------------------------ read-out
print("=" * 100)
print("WINNER BY |bias in b1|, per design cell")
for sa in SA_RANGES:
    row = f"  SA {sa:>4.0f} : "
    for fc in FRAC_COARSE:
        acc = summary[(sa, fc)]
        best = min(METHODS, key=lambda k: abs(np.mean(acc[k]["b"]) - B1_TRUE))
        row += f"{fc:.0%} coarse -> {best:<12}  "
    print(row)

print("\nSPREAD between best and worst |bias|, per cell")
for sa in SA_RANGES:
    row = f"  SA {sa:>4.0f} : "
    for fc in FRAC_COARSE:
        acc = summary[(sa, fc)]
        errs = [abs(np.mean(acc[k]["b"]) - B1_TRUE) for k in METHODS]
        row += f"{fc:.0%} -> {max(errs) - min(errs):.3f}   "
    print(row)
print("\nWhere the spread is small the choice of method does not matter, and")
print("that is the applicability boundary the paper is about. Read the row")
print("where the winner changes, not the winner itself.")
print("=" * 100)
