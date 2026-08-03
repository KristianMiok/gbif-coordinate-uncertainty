"""
The bottleneck is NOT the disc integration -- it is recomputing exp(beta0 + x*beta)
over the entire raster at every likelihood evaluation.

Two fixes, both exact-enough:
  1. sparse: only cells inside the UNION of discs enter the per-record integrals.
     Precompute that index set ONCE; per evaluation touch only those cells.
  2. coarse offset: the global term integral_S lambda(s) ds is a single scalar.
     It does not need native resolution. Evaluate it on a coarsened raster.
"""
import numpy as np
import time

rng = np.random.default_rng(0)
H = W = 2000                 # 4M cells, ~1km CHELSA over a large region
n = 100_000
radius_classes = [1, 5, 10]  # cells

x = rng.normal(size=(H, W)).astype(np.float64)
rows = rng.integers(0, H, n)
cols = rng.integers(0, W, n)
cls = rng.integers(0, len(radius_classes), n)

# ---------- one-off precomputation (done ONCE, not per evaluation) ----------
t0 = time.perf_counter()
flat_idx_per_record = []
R_max = max(radius_classes)
xpad = np.pad(x, R_max, mode="constant")
Hp = H + 2 * R_max

masks = {}
for r in radius_classes:
    R = int(np.ceil(r))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    m = (xx ** 2 + yy ** 2) <= r ** 2
    dy, dx = np.nonzero(m)
    masks[r] = (dy - R, dx - R)

all_idx = []
owner = []
for i in range(n):
    r = radius_classes[cls[i]]
    dy, dx = masks[r]
    rr = rows[i] + dy + R_max
    cc = cols[i] + dx + R_max
    all_idx.append(rr * Hp + cc)
    owner.append(np.full(len(dy), i))
all_idx = np.concatenate(all_idx)
owner = np.concatenate(owner)

uniq, inverse = np.unique(all_idx, return_inverse=True)
x_sparse = xpad.ravel()[uniq]
t_pre = time.perf_counter() - t0

# coarsened raster for the global offset
fac = 10
x_coarse = x[:H // fac * fac, :W // fac * fac].reshape(
    H // fac, fac, W // fac, fac).mean(axis=(1, 3))

print(f"raster {H}x{W} = {H*W/1e6:.1f}M cells | n = {n:,}")
print(f"union of discs touches {len(uniq)/1e6:.2f}M cells "
      f"({100*len(uniq)/(H*W):.1f}% of raster)")
print(f"one-off precomputation: {t_pre:.1f} s  (done once, reused every evaluation)")
print()

beta0, beta1 = 0.3, 0.5

# ---------- per-evaluation cost, dense (what the shipped code does) ----------
t0 = time.perf_counter()
lam_dense = np.exp(beta0 + beta1 * x)
offset_dense = lam_dense.mean()
t_dense = time.perf_counter() - t0

# ---------- per-evaluation cost, sparse + coarse offset ----------
t0 = time.perf_counter()
lam_sparse = np.exp(beta0 + beta1 * x_sparse)
per_cell = lam_sparse[inverse]
disc_int = np.bincount(owner, weights=per_cell, minlength=n)
offset_coarse = np.exp(beta0 + beta1 * x_coarse).mean()
ll = np.log(disc_int + 1e-300).sum() - offset_coarse
t_sparse = time.perf_counter() - t0

print(f"per likelihood evaluation:")
print(f"  dense  (exp over full raster)      : {t_dense*1000:7.1f} ms")
print(f"  sparse (union of discs + coarse Z) : {t_sparse*1000:7.1f} ms"
      f"   ({t_dense/t_sparse:.1f}x faster)")
print()
for iters in (2000,):
    print(f"  {iters} optimiser iterations -> dense {t_dense*iters/60:.1f} min"
          f" | sparse {t_sparse*iters/60:.1f} min")

print()
print(f"offset check: dense mean {offset_dense:.6f} vs coarse {offset_coarse:.6f} "
      f"(rel. diff {abs(offset_dense-offset_coarse)/offset_dense:.2%})")
print()
n_fits = 4 * 3 * 3 * 4 * 100   # methods x r-classes x directedness x delta x reps
print(f"factorial sweep of ~{n_fits:,} fits:")
print(f"  dense  : {t_dense*2000*n_fits/3600:,.0f} core-hours")
print(f"  sparse : {t_sparse*2000*n_fits/3600:,.0f} core-hours")
