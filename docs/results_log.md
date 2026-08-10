# Results ledger

One line per claim. **Nothing enters a manuscript from here without being
re-run.** Killed items stay in the file — they are the reason not to repeat the
work.

Last updated 2026-08-03.

---

## Established

| id | claim | evidence | script |
|---|---|---|---|
| C1 | Change of support is unbiased with **overlapping** per-record supports. Disjointness in Hefley (2017) / Walker (2020) is a property of their data (PLSS sections, administrative units), not a requirement of the estimator. | bias −0.001, SD 0.018, 95% coverage 0.965, n = 200 replicates; 94.1% of raster cells claimed by >1 disc, max 68 discs per cell | `src/python/cos_overlap_test.py` |
| C2 | The published COS likelihood **cannot estimate the radius**. ∫_{A_i}λ grows like r² while the offset ∫_S λ is r-free, so the profile is monotone and pins at the search bound. r becomes a parameter only under a normalised (reporting-model) kernel, where total mass is preserved. | first run returned r̂ = 32.00 (grid maximum) in all 9 cells; interior maximum after normalising | `src/python/p4_effort_confound.py` (see the docstring of `disc_kernel_fft`) |
| C3 | Identifiability of the radius tracks the ratio of radius to predictor autocorrelation range — the Naimi (2011) condition restated as an identifiability statement. | 2-log-unit width 0.02 → 0.17 → 0.81 as r/SA falls 4 → 1 → 0.25; predicted in advance in the pre-registration | `src/python/p4_effort_confound.py` |
| C4 | Empirical partition of GBIF crayfish at 1 km: 42.7% inert (radius smaller than half a cell, correction impossible and unnecessary), 28.7% no radius reported, 20.9% usable for COS, 7.7% candidate grid centroids. At 90 m the COS-relevant share rises to 25.6%. | 507,980 records, GBIF DOI **10.15468/dl.99ezk2** | `src/R/01_...`, `src/R/02_...` |
| C5a | Documented geocoder defaults (301, 999, 3036, 9999) survive a `r < 1000` filter at a **higher** rate than real measurements (79.4% vs 73.7%). 999 in particular is below the threshold and passes untouched. A threshold cannot distinguish a software artefact from a measurement. | 1,098 records carry a documented default; 3036 sits at the 75th percentile of the reported-radius distribution | `src/R/03_precision_support.R` |
| C5b | 17,137 records (4.7% of reporting) claim a radius smaller than their own coordinate rounding can carry. Median understatement 1.3×, 90th percentile 11.2×. Dominated by records claiming r = 5000 while given to one decimal. | — | `src/R/03_precision_support.R` |

## Killed

| what | why | date |
|---|---|---|
| "we avoid assuming an error distribution" as a novelty claim | uniform-on-disc **is** a distribution; the real distinction is geographic vs covariate space, and Hefley (2017) already specifies error geographically | 2026-08-03 |
| per-record heterogeneous radius as the differentiator vs `refitME` | `refitME`'s `sigma.sq.u` is per **variable**, not per observation, so refitME is homoscedastic across records — but Hefley's framework takes a per-record polygon natively, so this differentiates against the wrong opponent | 2026-08-03 |
| FFT convolution as a general speed-up | loses to a naive per-record loop below ~150k records; wins 5× at n = 1M. Ship naive as default, FFT above a measured threshold | 2026-08-03 |
| sparse evaluation on the union of discs | at n = 100k the union covers 97.3% of the raster; sparse indexing is 5–13× **slower** | 2026-08-03 |
| coarsening the global offset ∫_S λ | Jensen: averaging x then exponentiating underestimates by 11.6%, biasing β₀ | 2026-08-03 |
| **P4** — radius and sampling effort are inseparable | not confirmed. Where the profile is sharp (SA = 2, 8) effort moves r̂ within noise. The apparent drift was confined to the cell where r̂ has SD 5–12 on a 1.5–32 grid, i.e. where the estimator has no information. The pre-registered verdict rule was itself wrong (no identifiability precondition, no MCSE comparison) — corrected in Amendment 1 | 2026-08-03 |
| √2 grid-centroid signature | independent check failed: 72.3% coordinate duplication in the candidate subset vs 68.9% elsewhere. Baseline is saturated because crayfish are monitored at fixed localities, so the test had no power. **Do not build a lattice detector — GridDER (Ecological Informatics 2023) already does this.** | 2026-08-03 |
| "NA as an informative third state" as a **method** | The second available signal — coordinate decimals — is informative only in the few-decimals direction. For NA records the bound is a *lower* bound with median 0.7 m, i.e. no constraint. Using it as a support would manufacture false precision (δ_size < 1). Survives as a **finding** (28.7% is irreducibly indeterminate) and as an identifiability limit in the Discussion, not as a method | 2026-08-03 |

## Open

| id | question | next action |
|---|---|---|
| C6 | How does COS compare with the existing alternatives — threshold filter, naive point, and Smith et al. (2023) NGP/NEP — on the same virtual species? | **the keystone; nothing else should start before it** |
| — | Does the reporting-model kernel (C2) actually recover β₁ better than plain COS when the radius is misspecified? | narrow δ_size test; only if the draft needs C2 to carry weight |
| — | Second taxon. m20 was already criticised for one system; crayfish report a radius 71.3% of the time against Marcer et al.'s 18% for preserved specimens, so this dataset is **not** representative | Odonata is already in hand from m20 |
| — | Are the 17,137 internally inconsistent records (C5b) real false precision, or publishers generalising coordinates for protected species? | break down by `datasetKey` before any claim |

## Method notes worth keeping

- The correct term is **Berkson error**, not "location error": the truth lies
  inside the reported region, so the error is independent of the observed value.
  Berkson error is benign for linear models but not for nonlinear ones, and SDMs
  have a log link.
- Berkson *misspecification* already has a literature (arXiv 2306.01468 and the
  robust-Bayesian line). What is absent is the **spatial** case, where the
  support has geometry.
- Both shipped implementations aggregate by a partition key — `aggregate(raster,
  fact=)` in Walker's simulation, a `data.table` group-by on a rasterised
  section-ID layer in the CWD analysis. Neither can represent overlapping
  supports, so C1 is real at the level of what can be computed.
- Walker's reported 5.5 h runtime is not the cost of COS: it is 24 raster layers
  (4 ages × 2 sexes × marked bivariate structure) over 8.59M cells, in R. The
  presence-only equivalent runs in ~2 minutes in numpy on a laptop.

## Added 2026-08-03 (evening) — C6 closed

| id | claim | evidence |
|---|---|---|
| C6a | Ignoring reported uncertainty is the only error that matters for **bias**. Threshold, NEP and COS all recover b1 = 1.00 ± 0.01 across all nine design cells; naive attenuates to 0.43–0.81. | `src/python/c6_comparison.py`, 30 reps |
| C6b | The applicability boundary is real: spread between best and worst method is 0.185–0.566 when radius >> autocorrelation range, and 0.012–0.025 when radius << range. Factor of ~20. | same |
| C6c | On **RMSE**, COS wins all 12 cells of the sample-size sweep, but the gain tracks the FRACTION OF IMPRECISE RECORDS, not sample size: ~5% at 25% coarse, ~10% at 50%, ~25% at 75%, essentially flat in n. | `src/python/c6b_samplesize.py`, 200 reps |
| C6d | NEP (Smith et al. 2023) tracks the threshold filter almost exactly (RMSE 0.437 vs 0.447 at the hardest cell). It retains the records but not their information: assigning each imprecise record the environment nearest the precise-record mean adds near-mean covariate values, which carry almost nothing about the slope. | same |
| C6e | Range inflation is the decision-relevant metric, not RMSE. Naive inflates predicted suitable area by 64–72% at only 25% imprecise records, and 14–17x at 75%, independent of sample size. All three corrections hold inflation at 1.00. | both scripts |

**Calibration to the crayfish data — state this in the paper.** At 1 km, ~25% of
reporting records fall in the COS-relevant class, which puts the expected RMSE
gain of COS over a threshold filter at 4–7%. At 90 m the class is ~36% and the
gain perhaps 10%. The COS advantage on this dataset is modest. The strong claim
is about range inflation under the naive treatment, not about estimator choice
among the corrections.

**Framing correction.** The C6 design was built to test Smith et al.'s small-n
claim and used the wrong axis. Their regime (fewer than ~15–20 precise records)
turned out not to govern the comparison; the imprecise fraction does. Report the
imprecise fraction, not the sample size, as the applicability variable.

## Added 2026-08-03 (late) — the partition is a publisher property

| id | claim | evidence |
|---|---|---|
| C7 | The partition varies 11-fold across 13 taxonomic groups at 1 km (4.5% salmonids to 47.4% bats) and does not track realm or record type; within-freshwater variation (4.5–37.8%) exceeds between-category variation. | `src/R/04_multitaxon_partition.R`, accessed 2026-08-03 |
| C8 | **Within a single taxon, the between-dataset spread is larger than the between-taxon spread.** Usable share across 81 crayfish datasets (>=500 records, 93.1% of all records) runs 0%–100%, 10th–90th percentile 0%–84.4%. Dataset identity explains 52.2% of the variance in whether a record is usable. | `src/R/05_per_dataset.R` |
| C9 | Publishers operate in regimes, not gradients: five of the fifteen largest datasets report no radius at all, one is 86.6% inert with median r = 15 m, one is 79.2% usable with median r = 5000 m. | same |
| C5b | **REVISED — withdraw the original wording.** The 17,137 "internally inconsistent" records are concentrated: the top three datasets hold 79.4%, and one dataset (8a863029) holds 11,865 at a 44.2% internal rate. That same dataset has median r = 5000 m. This is a publisher georeferencing convention (locality-name georeferencing at a fixed 5 km radius, coordinates truncated on publication), NOT record-level false precision. Present it as a publisher property. | same |

**Consequence for the package.** The diagnostic must report PER DATASET, not
aggregated over a download. The aggregate 17.4% usable for crayfish conceals
that half the datasets sit at 0% and part of the corpus sits near 80%. An
aggregate number is actively misleading here.

**Consequence for the framing.** Two ecologists working on the same species from
different sources get opposite partitions, and neither can anticipate it. That
is the paper.

## Added 2026-08-03 — C8 replicates across taxa

| id | claim | evidence |
|---|---|---|
| C8r | **C8 replicates in all twelve groups.** Dataset identity explains 47.3%–82.8% of the variance in whether a record is usable (median 61.1%). Crayfish, at 51.9%, sit near the LOW end — they are a conservative case, not an outlier. | `src/R/06_c8_replication.R`, accessed 2026-08-03 |
| C9r | The bimodality is general, not a crayfish quirk: 30.4%–82.6% of datasets have essentially no usable records, while 3.2%–41.4% are essentially all usable. Every group has p10 = 0.0%. | same |
| — | **Method validation:** facet-based per-dataset counts reproduce the download-based measurement for crayfish (R2 51.9% vs 52.2%; usable 17.4% in both). The facet route measures the same quantity without downloading. | same |

**Caveats to carry into the manuscript.**
- Four groups (dragonflies, orchids, mosses, bats) hit the 500-dataset facet
  cap, so only the largest publishers are seen. Small publishers are more
  heterogeneous, so this biases R2 DOWNWARD — conservative, but state it, and
  re-run those four at a higher cap to show stability.
- Bats at 82.8% with 41.4% of datasets fully usable is worth explaining rather
  than reporting bare; likely acoustic monitoring with a fixed detection radius.
- **The R2 is optimistic by construction:** a one-way ANOVA on a binary outcome
  with up to 500 groups spends 500 degrees of freedom. Replace with an adjusted
  R2 or the ICC from a mixed model with dataset as a random effect before
  submission. The bimodality is too strong to be an artefact, but the headline
  number will move.

## Added 2026-08-03 — the ICC attempt FAILED, do not use it

`src/R/07_icc.R` returns ICC 92.8%–99.5% (median 98.6%) with sigma2_u from 42.5
to 642.8. **These are degenerate fits, not results.**

The logistic ICC is sigma2_u / (sigma2_u + pi^2/3), and pi^2/3 = 3.29. Once
sigma2_u reaches the hundreds the ICC is driven to 1 arithmetically, whatever
the data say. The cause is separation: a large share of datasets sit at exactly
k = 0 or k = n, the logit is +-Inf for those, and glmer pushes the random-effect
variance to the numerical boundary. The script did not check convergence
warnings, only errors, so the failure printed as a clean table.

**Neither variance model is usable on this outcome.** The one-way ANOVA R2 is
optimistic with hundreds of groups; the binomial GLMM degenerates under the same
bimodality that is the finding. Do not report either.

**Report description instead.** The bimodality is visible without a model:
- share of datasets below 1% usable: 30.4%–82.6% across taxa
- share above 99% usable: 3.2%–41.4%
- every group has p10 = 0.0%

If a single index is wanted, use a Gini coefficient or the entropy of the
per-dataset usable share, both of which are assumption-free and do not degenerate
under separation.

**Process note.** The ICC was proposed as the more defensible statistic without
first checking whether the model applies to data with this much separation. Check
model applicability before proposing a replacement statistic, not after.

## Added 2026-08-03 — the dataset effect, measured without a model

| id | claim | evidence |
|---|---|---|
| C8f | **Knowing the publisher removes a median 59% of the uncertainty about whether a record is usable.** Uncertainty coefficient U = 1 - H(usable|dataset)/H(usable) runs 0.459-0.803 across twelve groups (median 0.593). | `src/R/08_dataset_effect.R` |
| C8g | **None of it is a many-groups artefact.** A permutation null holding dataset sizes and the total usable count fixed gives U = 0.000 in all twelve groups; the excess equals the statistic. Validated on synthetic data: 600 datasets with no real effect give U = 0.002 against a null of 0.002. | same |
| C9f | 57.3%-88.0% of datasets (median 69.3%) sit in a pure regime -- either essentially nothing usable or essentially everything. Gini of the per-dataset usable share runs 0.516-0.931. | same |

**Report these, not the failed variance models.** The uncertainty coefficient
works precisely because separation is not a problem for it: a dataset at 0% or
100% contributes zero conditional entropy, which is the finding rather than an
obstacle. Worth a methods note -- anyone measuring a "publisher effect" on data
this bimodal will hit the same wall with ANOVA and with a GLMM.

**Figure 1 candidate:** `results/dataset_regimes.png`, twelve histograms of the
per-dataset usable share. The bimodality is visible without any model.

**Follow-up worth one query:** bats are the outlier (U = 0.803, 36.7% of datasets
fully usable). Likely acoustic monitoring published with a fixed detection
radius. Confirm before using it as an illustration.

## Added 2026-08-03 — Figure 1, and a claim that had to be narrowed

**The strongest single number in the project.** Share of datasets falling within
±10 points of their own group's aggregate usable share:

  cetaceans 2.6% | dragonflies 3.6% | bats 5.1% | ground beetles 5.5%
  orchids 6.0% | crayfish 6.3% | amphibians 7.5% | mosses 8.0%
  swallowtails 9.2% | freshwater mussels 15.9% | salmonids 83.8% | corals 88.6%

In ten of twelve groups the aggregate describes fewer than one dataset in ten.
For cetaceans, 97 of every 100 datasets are nowhere near the group mean.

*Caveat:* salmonids and corals look like exceptions but are not. Their aggregates
are 4.4% and 4.9%, i.e. already at the floor, so a ±10-point window catches the
zero spike by construction. Use a RELATIVE window before publication, or a
reviewer will find this.

**The zero spike splits into two opposite behaviours.** Share of the spike that
is "no radius reported" rather than "radius reported and small": 18.1%
(cetaceans) to 63.3% (swallowtails). For cetaceans the spike is genuinely
precise data needing no correction; for swallowtails it is publishers declining
to state uncertainty, where correction is needed but impossible. The exploratory
figure in 08_ conflated these.

**The regime claim had to be narrowed.** Weighting by DATASETS gives clean
bimodality. Weighting by RECORDS fills the middle: roughly 35-55% of records in
several groups sit at intermediate shares, because large datasets are more often
intermediate while fully-usable datasets are numerous but small. "Publishers
occupy regimes" is true of publishers, not of records. The claim that survives
both weightings is the narrower one: **the aggregate describes no dataset.**

Figures: `results/fig1_regimes_by_record.png` (main; what a user will actually
find in a download) and `results/fig1_regimes_by_dataset.png` (supplement; how
publishers behave).

## Added 2026-08-08 -- vocabulary and boundary harmonisation (found while building the ingrain package)

Four definitional mismatches surfaced by porting the classifier; none
changes a published number, all must be fixed or footnoted before review.

1. **04_ comment vs implementation.** Line 48 says `inert if r < R/2`;
   the cumulative-query implementation counts `r <= R/2` as inert
   (GBIF ranges inclusive). The implementation is the reference; fix the
   comment to `<=`.
2. **05_ differs from 04_ twice.** It uses strict `r < RES/2`, and it
   extracts a fifth state (geocoder defaults 301/999/3036/9999) BEFORE
   partitioning, so its per-dataset inert/marginal shares exclude
   records that 04_ includes (301 falls in 04_'s inert band at 1 km).
   Per-dataset and multi-taxon percentages are therefore not on the same
   definition. 05_ also nulls r <= 0 and r > 1e7.
3. **"Usable" means two opposite things.** 06_/07_/08_: usable
   = reported radius > 3R (the actionable class) over ALL records,
   datasets >= 500 -- this is the variable behind U = 0.593. 12_'s
   header ("at or below the threshold, as in the earlier scripts")
   mischaracterises the earlier scripts; 12_'s swept quantity r <= t is
   a different variable. Any manuscript sentence glossing U = 0.593
   must say: whether a record's reported radius exceeds three times the
   1 km grain.
4. **Package vocabulary (adopted in ingrain).** States: inert /
   marginal / actionable / unreported; "usable" retired from state
   names. The published binary U is exactly the collapse of the
   four-state U (actionable vs all other states pooled). "in-grain
   share" is reserved for the r <= t quantity from 12_, pending a
   decision on whether it enters the manuscript prose.

## Added 2026-08-10 -- pre-submission scoping sweep (kill criteria fixed in advance)

**K1 (pre-emption of the resolution-relative audit): NOT triggered.**
Tool and preprint space checked (CoordinateCleaner 3.0.1, GeoThinneR
2025, HOGS MEE 2025, arXiv/bioRxiv sweeps): nothing classifies
occurrence records by correction-relevance relative to the analysis
grain. Claim wording adopted: "first tool to operationalise a
resolution-relative partition of occurrence records by
correction-relevance", Hefley et al. (2017) credited as the seed.

**Gábor et al. 2022 (MEE, 10.1111/2041-210X.13956) must be cited and
pre-empted.** They show coarser analysis grains do not overcome
positional error -- superficially in tension with records migrating to
inert. Defence (now in the vignette): inert is a statement about
correctability, not harm; coarsening shrinks the share of records a
treatment can act on, not the damage.

**Citation chimera caught and fixed.** The vignette cited "Quality
issues in georeferencing" (Marcer et al. 2021, Diversity and
Distributions 27:564-567) under the venue/ID of "Uncertainty matters"
(Marcer et al. 2022, Ecography e06025). Same error class as the m20
Zizka citation. ACTION STANDING: verify every reference against
CrossRef before any submission. Claim-4 boundary confirmed: Marcer
2022 is descriptive, preserved-specimens, institution/country level;
ours is information-theoretic (U + permutation null), dataset-level,
resolution-relative, all record types.

**MEE Applications format (checked 2026-08-10):** 3000-4000 words
INCLUDING references; IMRAD not required; application description plus
worked examples recommended; Code Checklist at submission; software
tests expected. The vignette maps onto this near one-to-one.
