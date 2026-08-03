"""
Cost of one COS log-likelihood evaluation at GBIF scale.

Walker/Hefley shipped code recomputes lambda over the WHOLE raster each
evaluation, then does a group-by aggregation keyed on a section-ID layer.
That trick needs a partition (each cell -> exactly one polygon).

GBIF point-radius discs OVERLAP, so the group-by is unavailable.
Two replacements are benchmarked here:

  (A) naive: for each record, sum lambda over the cells in its own disc
  (B) FFT: integral over disc = convolution of lambda with a disc kernel,
      evaluated at the record centre. One convolution per DISTINCT RADIUS,
      independent of the number of records.

Only the per-evaluation cost is measured; multiply by optimiser iterations.
"""
import numpy as np
from scipy.signal import fftconvolve
import time

rng = np.random.default_rng(0)


def disc_kernel(r_cells):
    """Uniform disc kernel, normalised to integrate to the disc area (in cells)."""
    R = int(np.ceil(r_cells))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    k = ((xx ** 2 + yy ** 2) <= r_cells ** 2).astype(np.float64)
    return k


def naive_disc_sums(lam, rows, cols, r_cells):
    """Per-record loop: slice the disc bounding box, mask, sum."""
    R = int(np.ceil(r_cells))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    mask = (xx ** 2 + yy ** 2) <= r_cells ** 2
    H, W = lam.shape
    out = np.empty(len(rows))
    padded = np.pad(lam, R, mode="constant")
    for i in range(len(rows)):
        r0 = rows[i]
        c0 = cols[i]
        block = padded[r0:r0 + 2 * R + 1, c0:c0 + 2 * R + 1]
        out[i] = block[mask].sum()
    return out


def fft_disc_sums(lam, rows, cols, r_cells):
    """One convolution for the whole radius class, then point lookup."""
    k = disc_kernel(r_cells)
    conv = fftconvolve(lam, k, mode="same")
    return conv[rows, cols]


def run(grid, n_records, radius_classes, label):
    H = W = grid
    # a smooth-ish intensity surface; content is irrelevant to timing
    beta_x = rng.normal(size=(H, W)).astype(np.float64)
    lam = np.exp(0.5 * beta_x)  # this exp() is unavoidable per evaluation

    rows = rng.integers(0, H, n_records)
    cols = rng.integers(0, W, n_records)
    # assign each record to a radius class (GBIF radii pile up at few values)
    cls = rng.integers(0, len(radius_classes), n_records)

    print(f"\n{label}")
    print(f"  raster {H}x{W} = {H*W/1e6:.2f}M cells | n = {n_records:,} "
          f"| radius classes (cells) = {radius_classes}")

    t0 = time.perf_counter()
    lam_local = np.exp(0.5 * beta_x)
    t_exp = time.perf_counter() - t0
    print(f"  exp() over full raster            : {t_exp:8.3f} s")

    # ---- naive
    t0 = time.perf_counter()
    for j, r in enumerate(radius_classes):
        sel = cls == j
        if sel.sum():
            naive_disc_sums(lam, rows[sel], cols[sel], r)
    t_naive = time.perf_counter() - t0
    print(f"  (A) naive per-record disc sums    : {t_naive:8.3f} s")

    # ---- fft
    t0 = time.perf_counter()
    for j, r in enumerate(radius_classes):
        sel = cls == j
        if sel.sum():
            fft_disc_sums(lam, rows[sel], cols[sel], r)
    t_fft = time.perf_counter() - t0
    print(f"  (B) FFT, one conv per radius class: {t_fft:8.3f} s")

    total_naive = t_exp + t_naive
    total_fft = t_exp + t_fft
    print(f"  --> per likelihood eval: naive {total_naive:.3f} s | FFT {total_fft:.3f} s"
          f"  (speedup {total_naive/total_fft:.1f}x)")
    for iters in (2000,):
        print(f"  --> {iters} optimiser iterations: naive "
              f"{total_naive*iters/3600:8.2f} h | FFT {total_fft*iters/3600:8.2f} h")
    return total_naive, total_fft


# Walker CWD scale, but presence-only (1 layer, not their 24-layer stack)
run(2930, 2500, [8], "WALKER SCALE (2974x2890 raster, n=2497) -- presence-only equivalent")

# GBIF-realistic: CHELSA 1km over a large region, radii 1km/5km/10km
run(2000, 100_000, [1, 5, 10], "GBIF SCALE (1km raster, 4M cells, n=100k)")

run(2000, 1_000_000, [1, 5, 10], "GBIF SCALE (1km raster, 4M cells, n=1M)")
