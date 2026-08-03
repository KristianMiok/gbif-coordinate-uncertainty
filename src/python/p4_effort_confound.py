"""
P4 — is the uncertainty radius separable from unobserved sampling effort?

Pre-registered in prereg_support_identifiability.md. Do not change the design
after reading the output.

The observed reported locations are a realisation of

    ( lambda(s) * b(s) )  convolved with  K_{r_true}

lambda = exp(b0 + b1 x)  is the ecological signal
b(s)                     is sampling effort -- UNOBSERVED in practice
K_r                      is the uniform disc kernel of the true radius

The fitted model has no effort term (you cannot map collector effort), so it is

    lambda_beta  convolved with  K_{r_assumed}

Both effort structure and positional blurring make the observed pattern smoother
than lambda. The question is whether the likelihood can still tell them apart.

beta0 profiles out in closed form:
    ell(b1, r) = n log n - n log T - n + sum_i log g_i
    T   = sum_cells exp(b1 x)
    g_i = ( exp(b1 x) conv K_r ) evaluated at record i

so each evaluation costs one circular FFT convolution.

READ-OUT
  P1/P2 : under gamma = 0, does the profile in r peak at r_true, and does the
          peak sharpen as the autocorrelation range of x falls below r_true?
  P4    : as effort structure gamma rises, does the fitted radius drift away
          from r_true? Drift = the two are confounded = STRONG branch.
"""
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize_scalar
import time

# ------------------------------------------------------------------ settings
N = 256                       # raster side
R_TRUE = 8.0                  # true uncertainty radius, cells
B1_TRUE = 1.0
N_EXPECT = 800                # expected records
SA_RANGES = [2.0, 8.0, 32.0]  # autocorrelation range of x, vs R_TRUE = 8
GAMMAS = [0.0, 0.5, 1.0]      # strength of unobserved effort structure
SA_EFFORT = 16.0              # autocorrelation range of the effort field
R_GRID = np.geomspace(1.5, 32.0, 15)
N_REP = 15
SEED = 20260803

rng = np.random.default_rng(SEED)
FREQ = np.fft.rfft2(np.zeros((N, N))).shape


def smooth_field(rng, scale):
    z = gaussian_filter(rng.normal(size=(N, N)), scale, mode="wrap")
    return (z - z.mean()) / z.std()


def disc_kernel_fft(r):
    """FFT of a disc kernel NORMALISED TO UNIT MASS.

    This normalisation is the whole point. The raw COS integral
    integral_{A_i} lambda grows like r^2, so sum_i log(.) increases without
    bound in r while the offset integral_S lambda does not depend on r at all.
    The published COS likelihood is therefore MONOTONE in r: it is conditional
    on the reporting mechanism and cannot compare different assumed radii.

    To make r a likelihood parameter you must model the REPORTING process --
    the density of the reported location w given the true location u, uniform
    on disc(u, r) and integrating to one. The observed points are then an IPP
    with intensity lambda conv Kbar_r, which preserves total mass, so the
    offset stays r-free and the profile in r has an interior maximum.
    """
    yy, xx = np.mgrid[0:N, 0:N]
    yy = np.minimum(yy, N - yy)
    xx = np.minimum(xx, N - xx)
    k = ((xx ** 2 + yy ** 2) <= r ** 2).astype(np.float64)
    k /= k.sum()
    return np.fft.rfft2(k)


KFFT = {float(r): disc_kernel_fft(r) for r in R_GRID}
KFFT[R_TRUE] = disc_kernel_fft(R_TRUE)


def circ_conv(field, kfft):
    return np.fft.irfft2(np.fft.rfft2(field) * kfft, s=(N, N))


def simulate(x, gamma, z, rng):
    """Points from lambda*b, then displaced uniformly inside disc(R_TRUE)."""
    lam = np.exp(B1_TRUE * x)
    b = np.exp(gamma * z)
    inten = lam * b
    inten = inten / inten.sum() * N_EXPECT
    n = rng.poisson(inten.sum())
    idx = rng.choice(inten.size, size=n, replace=True,
                     p=(inten / inten.sum()).ravel())
    tr, tc = np.unravel_index(idx, (N, N))
    # uniform displacement inside the true disc
    R = int(np.ceil(R_TRUE))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    m = (xx ** 2 + yy ** 2) <= R_TRUE ** 2
    dy, dx = yy[m], xx[m]
    pick = rng.integers(0, len(dy), n)
    return (tr + dy[pick]) % N, (tc + dx[pick]) % N


def profile_ll(b1, rows, cols, x, kfft):
    """Profile log-likelihood with beta0 concentrated out."""
    e = np.exp(b1 * x)
    T = e.sum()
    g = circ_conv(e, kfft)[rows, cols]
    if np.any(g <= 0):
        return -np.inf
    n = len(rows)
    return n * np.log(n) - n * np.log(T) - n + np.log(g).sum()


def profile_over_r(rows, cols, x):
    """For each assumed radius, maximise over b1. Returns (ll, b1hat) arrays."""
    lls, b1s = [], []
    for r in R_GRID:
        res = minimize_scalar(
            lambda b: -profile_ll(b, rows, cols, x, KFFT[float(r)]),
            bounds=(0.0, 4.0), method="bounded",
            options=dict(xatol=1e-3))
        lls.append(-res.fun)
        b1s.append(res.x)
    return np.array(lls), np.array(b1s)


def support_width(lls):
    """Width of the 2-log-unit region in r, as a fraction of the grid span."""
    inside = lls >= lls.max() - 2.0
    lo, hi = R_GRID[inside][0], R_GRID[inside][-1]
    return lo, hi, np.log(hi / lo) / np.log(R_GRID[-1] / R_GRID[0])


print(f"raster {N}x{N} | r_true = {R_TRUE} cells | b1_true = {B1_TRUE}")
print(f"expected n = {N_EXPECT} | {N_REP} replicates per cell\n")
print(f"{'SA range':>9}{'gamma':>7}{'r_hat':>9}{'r_hat SD':>10}"
      f"{'b1_hat':>9}{'2-unit interval':>20}{'width':>8}")
print("-" * 72)

t0 = time.perf_counter()
results = {}

for sa in SA_RANGES:
    for gamma in GAMMAS:
        rh, bh, widths, los, his = [], [], [], [], []
        for rep in range(N_REP):
            x = smooth_field(rng, sa)
            z = smooth_field(rng, SA_EFFORT)
            rows, cols = simulate(x, gamma, z, rng)
            lls, b1s = profile_over_r(rows, cols, x)
            k = int(np.argmax(lls))
            rh.append(R_GRID[k])
            bh.append(b1s[k])
            lo, hi, w = support_width(lls)
            los.append(lo); his.append(hi); widths.append(w)
        rh = np.array(rh); bh = np.array(bh)
        results[(sa, gamma)] = (rh, bh, np.array(widths))
        print(f"{sa:>9.0f}{gamma:>7.1f}{rh.mean():>9.2f}{rh.std():>10.2f}"
              f"{bh.mean():>9.3f}"
              f"{f'[{np.mean(los):.1f}, {np.mean(his):.1f}]':>20}"
              f"{np.mean(widths):>8.2f}")

print(f"\nelapsed {time.perf_counter()-t0:.0f}s")

# ------------------------------------------------------------------ verdict
print("\n" + "=" * 72)
print("P1/P2 — identifiability vs autocorrelation range (gamma = 0)")
for sa in SA_RANGES:
    rh, _, w = results[(sa, 0.0)]
    tag = "SHARP" if w.mean() < 0.25 else ("BROAD" if w.mean() < 0.6 else "RIDGE")
    print(f"  SA range {sa:>4.0f} (r_true/SA = {R_TRUE/sa:>4.1f}) : "
          f"r_hat {rh.mean():5.2f}  2-unit width {w.mean():.2f}  -> {tag}")

print("\nP4 — drift of r_hat with unobserved effort structure")
for sa in SA_RANGES:
    base = results[(sa, 0.0)][0].mean()
    line = f"  SA range {sa:>4.0f}: "
    for gamma in GAMMAS:
        rh = results[(sa, gamma)][0].mean()
        line += f"gamma={gamma:.1f} -> r_hat {rh:5.2f} ({rh/base:4.2f}x)   "
    print(line)

drift = max(abs(np.log(results[(sa, GAMMAS[-1])][0].mean()
                       / results[(sa, 0.0)][0].mean())) for sa in SA_RANGES)
ridge = np.mean([results[(sa, 0.0)][2].mean() for sa in SA_RANGES]) > 0.6

print("\n" + "=" * 72)
if drift > np.log(1.5):
    print("STRONG branch — r_hat drifts with unobserved effort.")
    print("  Radius misspecification and sampling bias are confounded.")
    print("  The reported radius cannot be validated from the data under the")
    print("  normal GBIF condition. This is an impossibility result -> MEE.")
elif ridge:
    print("STRONG/WEAK — profile is a ridge even without effort structure.")
    print("  r is not estimable; check whether that is decision-relevant.")
else:
    print("WEAK/KILL branch — r_hat is stable under effort structure.")
    print("  The two are separable; the deliverable is a diagnostic, not an")
    print("  impossibility result. Ecography, and do not oversell it.")
print("=" * 72)
