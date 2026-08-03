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
R_COARSE = 8.0             # imprecise records, in cells
SA = 2.0                   # radius >> autocorrelation: the regime where the
                           # choice of method demonstrably matters (C6)
N_EXPECTS = [20, 50, 150, 500]
FRAC_COARSE = [0.25, 0.50, 0.75]
N_REP = 200
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


def simulate(x, frac_coarse, n_expect, rng):
    lam = np.exp(B1_TRUE * x)
    lam = lam / lam.sum() * n_expect
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
print(f"autocorrelation range fixed at {SA:.0f} cells (r/SA = {R_COARSE/SA:.0f})")
print(f"{N_REP} replicates per cell\n")
print("Smith et al. (2023) report their largest gains from retaining imprecise")
print("records when PRECISE records number fewer than about 15-20. The C6 run")
print("used n = 500, far outside that regime, so it could not test their claim.")
print("This run sweeps sample size to reach it.\n")

hdr = f"{'n_exp':>6}{'coarse':>8}{'n_prec':>8} | " + " | ".join(f"{m:^19}" for m in METHODS)
print(hdr)
print(f"{'':>6}{'':>8}{'':>8} | " + " | ".join(f"{'b1':>6}{'SD':>6}{'infl':>6}" for _ in METHODS))
print("-" * len(hdr))

t0 = time.perf_counter()
summary = {}

for ne in N_EXPECTS:
    for fc in FRAC_COARSE:
        acc = {k: {"b": [], "infl": []} for k in METHODS}
        nprec = []
        for _ in range(N_REP):
            x = field(rng, SA)
            rows, cols, coarse = simulate(x, fc, ne, rng)
            if len(rows) < 8 or (~coarse).sum() < 3:
                continue
            nprec.append(int((~coarse).sum()))

            b = fit_points(x, x[rows[~coarse], cols[~coarse]])
            acc["threshold"]["b"].append(b); acc["threshold"]["infl"].append(range_inflation(x, b))

            b = fit_points(x, x[rows, cols])
            acc["naive"]["b"].append(b); acc["naive"]["infl"].append(range_inflation(x, b))

            b = fit_points(x, nep_environments(x, rows, cols, coarse))
            acc["NEP (Smith)"]["b"].append(b); acc["NEP (Smith)"]["infl"].append(range_inflation(x, b))

            b = fit_cos(x, rows, cols, coarse)
            acc["COS"]["b"].append(b); acc["COS"]["infl"].append(range_inflation(x, b))

        summary[(ne, fc)] = (acc, np.mean(nprec))
        line = f"{ne:>6}{fc:>8.0%}{np.mean(nprec):>8.0f} | "
        line += " | ".join(
            f"{np.mean(acc[k]['b']):>6.2f}{np.std(acc[k]['b']):>6.2f}"
            f"{np.median(acc[k]['infl']):>6.2f}" for k in METHODS)
        print(line)

print(f"\nelapsed {time.perf_counter() - t0:.0f}s")
print("b1_true = 1.00; inflation is the MEDIAN, robust to the heavy tail at")
print("small n where a single bad fit dominates the mean.\n")

print("=" * len(hdr))
print("RMSE of b1 -- the quantity that trades bias against the variance cost")
print("of discarding records. This is the comparison Smith et al. are making.")
print(f"\n{'n_prec':>8}{'coarse':>8} | " + " | ".join(f"{m:>13}" for m in METHODS))
for ne in N_EXPECTS:
    for fc in FRAC_COARSE:
        acc, npr = summary[(ne, fc)]
        rm = {k: np.sqrt(np.mean((np.array(acc[k]["b"]) - B1_TRUE) ** 2)) for k in METHODS}
        best = min(rm, key=rm.get)
        row = f"{npr:>8.0f}{fc:>8.0%} | " + " | ".join(
            f"{rm[k]:>13.3f}" for k in METHODS)
        print(row + f"   <- {best}")

print("\nIf `threshold` stops winning as n_prec falls below ~15-20, Smith et")
print("al.'s claim is reproduced and discarding is costly at small samples.")
print("If it keeps winning, their result does not transfer to this design and")
print("the difference must be explained, not ignored.")
print("=" * len(hdr))
