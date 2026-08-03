# Paper skeleton — coordinate uncertainty in GBIF: for which records does it matter?

Draft framing, 2026-08-03. Written before any further simulation, to test
whether the results already in hand carry a paper and what is still missing.

Working reference: this is the GBIF/reliability paper. Sits downstream of m20
(quality filtering as covariate shift) and applies the same argument to a
continuous metadata field rather than a binary flag.

---

## 1. The argument in one paragraph

GBIF records carry a per-record positional uncertainty radius. Standard practice
is to threshold on it and discard the rest. Two published alternatives exist and
do not cite each other properly: a biogeographic line that assigns imprecise
records to a single conservative point (Smith et al. 2023 GEB, `enmSdmX`), and a
statistical line that integrates the intensity over the region containing the
truth (Hefley et al. 2017 MEE; Walker et al. 2020 Biometrics). We connect them,
extend the statistical correction to the per-record overlapping supports that
GBIF actually supplies, and then show that the correction is inert for most
records, unavailable for a large minority, and that the metadata deciding which
is which cannot be validated from the data. The recommendation is not a better
threshold but a reported partition.

## 2. Candidate titles

- *For which occurrence records does coordinate uncertainty matter?*
- *Coordinate uncertainty in GBIF is usable, inert, or unknowable: a partition*
- *Positional uncertainty as Berkson error: what can and cannot be corrected in
  aggregated occurrence data*

Prefer a title that names the partition. The correction is the method; the
partition is the finding.

## 3. Contributions, with status

| # | Claim | Status | Evidence |
|---|---|---|---|
| C1 | Change of support extends to per-record **overlapping** supports, enabling GBIF point-radius data for the first time | **DONE** | `cos_overlap_test.py`: bias −0.001, coverage 0.965, n=200; 94% of cells shared by >1 disc |
| C2 | The published COS likelihood is conditional on the reporting mechanism and cannot estimate the radius; a normalised (reporting-model) kernel is required | **DONE, theory + demo** | monotone profile, r̂ pinned at grid bound; interior maximum after normalisation |
| C3 | Identifiability of the radius tracks the ratio of radius to predictor autocorrelation range | **DONE** | 2-unit width 0.02 → 0.17 → 0.81 as r/SA falls 4 → 1 → 0.25 |
| C4 | Empirical partition: usable / inert / unreported | **DONE** | 507,980 crayfish records, DOI 10.15468/dl.99ezk2 |
| C5 | The threshold filter selects on an unverifiable label | **PARTLY** | filter drops 103,019 honest reporters; 145,710 non-reporters pass or fail on `na.rm` alone |
| C6 | Comparison against the existing alternatives | **NOT DONE** | needs threshold vs naive vs NGP/NEP vs COS on virtual species |

## 4. Section structure

**Introduction.** Positional uncertainty is Berkson error, not classical: the
truth lies inside the reported region, independent of the reported value. Berkson
error is benign in linear models but not in nonlinear ones, and SDMs have a log
link — so the common intuition that positional error "just adds noise" does not
apply. Two literatures, no contact between them. Smith et al. 2023 cites Hefley
2017 but classifies it among simulation studies that add artificial error, which
it is not.

**Methods 1 — the estimator.** IPP; COS; the per-observation form; proof that
disjointness is a property of the published datasets, not a requirement.
Reporting-model kernel and why the radius only becomes a parameter under it.

**Methods 2 — applicability.** The inert regime: when the support sits inside a
raster cell the intensity is constant and no correction is possible. This is
Hefley's own caveat, promoted here to a design criterion.

**Results 1 — simulation.** C1, C2, C3. Plus C6 once run.

**Results 2 — the partition.** C4, C5 on the crayfish download, at three raster
resolutions.

**Discussion.** What to report instead of a threshold. The unreported third
state. Freshwater is the favourable case: at 90 m the COS-relevant share more
than doubles relative to 1 km.

## 5. What is missing, in priority order

1. **C6, the comparison.** Without it the paper is four negative-ish results and
   a table; a reviewer will ask what they should actually do. This is the
   benchmark previously scoped as "Paper A" and it must be inside this paper,
   not a separate one.
2. **A second taxon.** m20 was already criticised for resting on one system.
   Repeating that here is avoidable: pick a group with a different uncertainty
   regime (Odonata is already in hand from m20; a preserved-specimen group would
   contrast better, since crayfish are 71% reporting against Marcer's 18% for
   specimens).
3. **GridDER** on the download, to settle the grid-centroid subset without
   writing a lattice detector.
4. **Does the reporting-model kernel actually help?** One narrow δ_size test.
   Only worth running if the draft needs C2 to carry weight rather than sit as a
   remark.

## 6. Vulnerabilities to state before a reviewer does

- C1 may be judged obvious. Pre-empt: the published implementations both use a
  partition-keyed aggregation and neither can represent overlapping supports, so
  the extension is real at the level of what can be computed, whatever its
  difficulty as mathematics.
- The paper is largely about when a method does *not* apply. That is defensible
  only if the applicability criterion is quantitative and actionable. Keep the
  partition table as the central object.
- Crayfish reporting rate (71.3%) is far above the GBIF norm for preserved
  specimens (18%, Marcer et al. 2022). State this; do not let the partition be
  read as a global figure.
- The 42.7% "inert" share is resolution-dependent by construction. Report all
  three resolutions, never one.
- Nothing here shows COS improves a real-data SDM. Either add that or say
  plainly that the paper is about applicability, not performance.

## 7. Target

With C6 and a second taxon: **Methods in Ecology and Evolution**, on the strength
of C1+C2 as method and the partition as the applied case.

Without C6: **Ecography** or **Diversity and Distributions**, as an applicability
and data-quality paper. Do not submit to a methods journal on C4 alone.

## 8. Decision to take now

Write the introduction and the partition section first, from results already in
hand. If the draft reads as complete without C6, run only the narrow δ_size
test. If it reads as unfinished — which is the likely outcome — schedule C6 as
the next block of work and do not start the package until it is done.
