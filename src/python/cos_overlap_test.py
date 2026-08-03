"""
DECIDING TEST — does the per-observation COS estimator stay unbiased when the
supports OVERLAP?

Hefley (2017 MEE) and Walker (2020 Biometrics) both derive the change of support
from a PARTITION: S = union of disjoint A_j, counts per polygon ~ Poisson.
Both shipped implementations exploit that (raster `aggregate` with a factor, or
a data.table group-by on a rasterised section-ID layer). Neither can represent
GBIF point-radius discs, which overlap.

Web Appendix A also admits a second route: factorise the IPP density over
observations and marginalise each u_i over its own A_i. That route never uses
disjointness. This script tests whether that is actually true in practice.

  log lambda(s) = b0 + b1 * x(s)
  ell(b)        = sum_i log( integral_{A_i} lambda ) - integral_S lambda

Four estimators per replicate:
  oracle  — true locations (upper bound on what is achievable)
  naive   — reported locations, error ignored (should be attenuated)
  COS-part— supports are disjoint grid cells (the published, known-good case)
  COS-disc— supports are OVERLAPPING discs (the case that decides the project)

PASS  : COS-disc bias ~ 0 and 95% CI coverage ~ 0.95
FAIL  : COS-disc biased, or coverage far from nominal
        -> the estimator needs disjointness after all, project stops here.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
import time

# ----------------------------------------------------------------- settings
H = W = 200          # raster cells
SMOOTH = 6.0         # spatial autocorrelation of the covariate (cells)
B0_TRUE = -3.0       # controls n
B1_TRUE = 1.0        # the parameter under test
R_DISC = 6.0         # uncertainty radius, in cells
N_REP = 200
SEED = 20260803

rng = np.random.default_rng(SEED)


def make_covariate(rng):
    z = gaussian_filter(rng.normal(size=(H, W)), SMOOTH, mode="wrap")
    return (z - z.mean()) / z.std()


def disc_offsets(r):
    R = int(np.ceil(r))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    m = (xx ** 2 + yy ** 2) <= r ** 2
    return yy[m], xx[m]


DY, DX = disc_offsets(R_DISC)
NDISC = len(DY)


def simulate(x, rng):
    """Draw an IPP realisation, then displace each point within its disc."""
    lam = np.exp(B0_TRUE + B1_TRUE * x)
    n = rng.poisson(lam.sum())
    if n < 30:
        return None
    p = (lam / lam.sum()).ravel()
    idx = rng.choice(lam.size, size=n, replace=True, p=p)
    tr, tc = np.unravel_index(idx, (H, W))
    # reported = true + uniform displacement inside the disc
    # (so the true location always lies inside the reported disc)
    pick = rng.integers(0, NDISC, n)
    rr = (tr + DY[pick]) % H
    cc = (tc + DX[pick]) % W
    return tr, tc, rr, cc


def disc_cells(rows, cols):
    """(n, NDISC) index arrays for the cells inside each reported disc."""
    R = (rows[:, None] + DY[None, :]) % H
    C = (cols[:, None] + DX[None, :]) % W
    return R, C


BLK = int(round(np.sqrt(NDISC)))       # square block of the same area as a disc
BY, BX = np.mgrid[0:BLK, 0:BLK]
BY, BX = BY.ravel(), BX.ravel()


def block_cells(rows, cols):
    """Disjoint square blocks tiling the raster — the published partition case."""
    r0 = (rows // BLK) * BLK
    c0 = (cols // BLK) * BLK
    R = (r0[:, None] + BY[None, :]) % H
    C = (c0[:, None] + BX[None, :]) % W
    return R, C


def nll_point(par, x, rows, cols):
    b0, b1 = par
    lam = np.exp(b0 + b1 * x)
    return -(np.sum(b0 + b1 * x[rows, cols]) - lam.sum())


def nll_cos(par, x, R, C):
    b0, b1 = par
    lam = np.exp(b0 + b1 * x)
    inner = lam[R, C].sum(axis=1)          # integral over each support
    return -(np.sum(np.log(inner)) - lam.sum())


def fit(fun, args, start=(-3.0, 0.5)):
    res = minimize(fun, start, args=args, method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
    # finite-difference Hessian for a Wald interval on b1
    h = 1e-4
    p = res.x
    f = lambda a, b: fun(np.array([a, b]), *args)
    d2 = (f(p[0], p[1] + h) - 2 * f(p[0], p[1]) + f(p[0], p[1] - h)) / h ** 2
    dab = (f(p[0] + h, p[1] + h) - f(p[0] + h, p[1] - h)
           - f(p[0] - h, p[1] + h) + f(p[0] - h, p[1] - h)) / (4 * h ** 2)
    d2a = (f(p[0] + h, p[1]) - 2 * f(p[0], p[1]) + f(p[0] - h, p[1])) / h ** 2
    Hm = np.array([[d2a, dab], [dab, d2]])
    try:
        se = np.sqrt(np.diag(np.linalg.inv(Hm)))[1]
    except np.linalg.LinAlgError:
        se = np.nan
    return res.x[1], se


# ------------------------------------------------------------------ overlap
x0 = make_covariate(rng)
sim0 = simulate(x0, rng)
R0, C0 = disc_cells(sim0[2], sim0[3])
flat = (R0 * W + C0).ravel()
_, counts = np.unique(flat, return_counts=True)
print(f"raster {H}x{W} | disc radius {R_DISC} cells ({NDISC} cells per disc)")
print(f"n = {len(sim0[0])} records")
print(f"cells claimed by >1 disc : {100*np.mean(counts > 1):.1f}% "
      f"(max {counts.max()} discs share one cell)")
print(f"--> supports genuinely overlap; a partition does not exist\n")

# ------------------------------------------------------------- replications
res = {k: [] for k in ("oracle", "naive", "cos_part", "cos_disc")}
cov = {k: [] for k in res}
t0 = time.perf_counter()

for rep in range(N_REP):
    x = make_covariate(rng)
    sim = simulate(x, rng)
    if sim is None:
        continue
    tr, tc, rr, cc = sim

    b, se = fit(nll_point, (x, tr, tc));  res["oracle"].append(b);  cov["oracle"].append(se)
    b, se = fit(nll_point, (x, rr, cc));  res["naive"].append(b);   cov["naive"].append(se)

    # published case: supports are DISJOINT square blocks (Walker's raster
    # `aggregate`), sized to match the disc so the comparison is fair
    Rp, Cp = block_cells(rr, cc)
    b, se = fit(nll_cos, (x, Rp, Cp));    res["cos_part"].append(b); cov["cos_part"].append(se)

    # the case under test: overlapping discs
    Rd, Cd = disc_cells(rr, cc)
    b, se = fit(nll_cos, (x, Rd, Cd));    res["cos_disc"].append(b); cov["cos_disc"].append(se)

    if (rep + 1) % 25 == 0:
        print(f"  {rep+1}/{N_REP} replicates ... {time.perf_counter()-t0:.0f}s")

print(f"\n{'estimator':<12}{'mean b1':>10}{'bias':>10}{'SD':>9}{'95% cov':>10}")
print("-" * 51)
for k in ("oracle", "naive", "cos_part", "cos_disc"):
    b = np.array(res[k]); s = np.array(cov[k])
    ok = np.isfinite(b) & np.isfinite(s)
    covg = np.mean(np.abs(b[ok] - B1_TRUE) < 1.96 * s[ok])
    print(f"{k:<12}{b[ok].mean():>10.3f}{b[ok].mean()-B1_TRUE:>10.3f}"
          f"{b[ok].std():>9.3f}{covg:>10.3f}")

bd = np.array(res["cos_disc"]); sd = np.array(cov["cos_disc"])
ok = np.isfinite(bd) & np.isfinite(sd)
bias = bd[ok].mean() - B1_TRUE
covg = np.mean(np.abs(bd[ok] - B1_TRUE) < 1.96 * sd[ok])
mcse = bd[ok].std() / np.sqrt(ok.sum())

print(f"\ntrue b1 = {B1_TRUE}   |   Monte Carlo SE of the bias = {mcse:.3f}")
print("\n" + "=" * 51)
if abs(bias) < 2 * mcse and 0.90 <= covg <= 0.98:
    print("PASS — COS is unbiased with OVERLAPPING supports.")
    print("       The partition in Hefley/Walker is a property of their data,")
    print("       not a requirement of the estimator. GBIF discs are usable.")
else:
    print("FAIL — COS degrades under overlap.")
    print(f"       bias {bias:+.3f} (MCSE {mcse:.3f}), coverage {covg:.3f}")
    print("       Disjointness is load-bearing. Stop before writing the package.")
print("=" * 51)
