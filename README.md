# Coordinate uncertainty in GBIF: for which records does it matter?

Working repository. **Not a finished study** — one keystone result (C6, the
comparison against existing alternatives) is missing. See
`docs/results_log.md` for exactly what is established, what has been killed,
and what remains.

Downstream of the "cleaning becomes biasing" line: the same argument applied to
a continuous metadata field rather than a binary quality flag.

## Layout

```
docs/    prereg_support_identifiability.md   pre-registration + dated amendment
         results_log.md                      ledger: established / killed / open
         paper_skeleton.md                   framing, contributions, gaps

src/R/   01_download_and_radius_profile.R    GBIF download + radius distribution
         02_grid_centroids.R                 grid-centroid candidates, partition
         03_precision_support.R              coordinate decimals as a bound

src/python/
         cos_overlap_test.py                 DECIDING TEST: COS under overlap
         p4_effort_confound.py               identifiability of the radius
         benchmarks/                         scaling; both reported a negative

results/ script output (gitignored except summaries)
data/    see data/README.md — the data is not in the repo
```

## Reproducing

R side, in order. `01` needs GBIF credentials in `~/.Renviron`
(`GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL`) and blocks while GBIF prepares the
download; `02` and `03` read the zip left on disk.

```bash
Rscript src/R/01_download_and_radius_profile.R
Rscript src/R/02_grid_centroids.R
Rscript src/R/03_precision_support.R
```

Python side, independent of the data:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python src/python/cos_overlap_test.py      # ~90 s
python src/python/p4_effort_confound.py    # ~10 s
```

## Data

Not committed. GBIF download DOI **10.15468/dl.99ezk2** — 507,980 crayfish
occurrence records (Astacidae, Cambaridae, Cambaroididae, Parastacidae) with
coordinates and no geospatial issue. The zip is ~49 MB; regenerate it with
`01_...R` or fetch it by key.

## How to read a negative result here

Several scripts print a verdict. Two of those verdicts have already been wrong:

- `p4_effort_confound.py` printed `STRONG` on the strength of a design cell in
  which the estimator carries no information. The pre-registered rule lacked an
  identifiability precondition and an MCSE comparison. Amendment 1 corrects it;
  the printed banner does **not** yet implement the corrected rule.
- `03_precision_support.R` originally asserted a direction that its own table
  contradicted. Fixed — it now computes retention rates instead of asserting.

Read the numbers, not the banner.

## Conventions

- Prior-art scan before code. Directions have been lost in this project by
  reversing that order.
- Kill criteria written down before the run, in `docs/`.
- Amendments are appended and dated; predictions are never edited in place.
