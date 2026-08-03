# Pre-registration — is a misspecified uncertainty radius identifiable?

**Written before any code is run.** Predictions and kill criteria are fixed here
so that whatever comes out is not a post-hoc fit.

Date: 2026-08-03
Depends on: `cos_overlap_test.py` (PASS — COS is unbiased with overlapping
supports, bias −0.001, coverage 0.965, n = 200 replicates)

---

## 1. Question

Change of support (Hefley 2017 MEE; Walker 2020 Biometrics) integrates the
intensity over a region assumed to contain the true location. GBIF supplies that
region as `coordinateUncertaintyInMeters`, a point-radius (Wieczorek 2004).

**The assumption is that the reported radius is correct. It frequently is not.**
GridDER finds ~13.5% of GBIF records carry grid-centroid coordinates whose
reported radius does not describe the real support. Records reporting `r = 1000`
that were georeferenced from a locality name are the same failure.

> Q1. Does COS remain useful when the reported support is wrong?
> Q2. Can the error be **detected from the data alone**, without ground truth?

Q1 is a magnitude question. Q2 is the one that decides the paper's ceiling.

## 2. Why it is decision-relevant

The standard practice is a threshold filter (`r < 1000`). That filter *trusts*
the reported radius: it keeps records claiming precision and drops records
honest enough to report a large radius. If reported radii cannot be validated,
the filter is selecting on an unverifiable label — which is the m20 thesis
applied to a continuous field.

## 3. Formal statement

Observed reported locations are a realisation of an IPP with intensity

    λ ⊛ K_{r_true}

where K_r is the uniform disc kernel. COS with an assumed radius fits

    λ_β ⊛ K_{r_assumed}

Identifiability of r is therefore a deconvolution question: does

    λ_β ⊛ K_{r_a}  =  λ_true ⊛ K_{r_t}

admit solutions with (β, r_a) ≠ (β_true, r_t)?

For small β₁, exp(β₁x) ≈ 1 + β₁x and convolution is linear, giving model
intensity ∝ 1 + β₁·x_r where x_r is x smoothed at scale r. When r is small
relative to the autocorrelation range of x, x_r ≈ x for all such r, so different
radii are near-indistinguishable and β₁ absorbs the difference.

## 4. Predictions, fixed in advance

**P1.** r is weakly or non-identifiable when r ≪ autocorrelation range of the
predictor; the profile likelihood in (r, β₁) is a ridge, not a peak.

**P2.** r becomes identifiable as r approaches or exceeds the autocorrelation
range. The identifiable regime coincides with the regime where COS matters.

**P3 (asymmetry).** Under-reported radius (false precision, r_a ≪ r_t) is more
detectable than over-reported radius, because the observed pattern is smoother
than any model with small r_a can reproduce, whereas an over-large r_a only
dilutes.

**P4 (the confound, and the real risk).** Excess smoothness in the observed
pattern is also produced by spatially structured sampling effort b(s). Radius
misspecification and effort bias enter the observed intensity in the same place.
If they are not separable, r is unidentifiable *for reasons that have nothing to
do with autocorrelation*, and P1–P3 are beside the point.

**P4 is the prediction most likely to make this a strong paper and least likely
to be true in a convenient way. It is tested first.**

## 5. Design

Factors, crossed:

- `r_a / r_t` ∈ {0.25, 0.5, 1, 2, 4}  — size misspecification
- support offset ∈ {0, 0.5·r, 1·r}    — shift misspecification (snapping)
  *(kept separate from size; conflating them is the flaw in the exploratory
  `cos_part` result and it must not be repeated)*
- autocorrelation range of x ∈ {r/4, r, 4r}
- sampling effort b(s): {uniform, smooth gradient, patchy}

Estimand: β₁ bias, 95% CI coverage, and the shape of the profile log-likelihood
surface over (r_a, β₁).

Identifiability is read off the profile surface: curvature of the profile in r_a
at the maximum, and the width of the 2-unit log-likelihood region.

## 6. Kill criteria

**KILL** — if under a uniform effort field the profile likelihood in r_a has a
clear peak across all autocorrelation regimes, r is straightforwardly estimable,
there is no identifiability story, and the contribution collapses to "estimate r
instead of trusting it". That is a short methods note, not a paper. Stop.

**WEAK (→ Ecography)** — if identifiability follows P1/P2, i.e. r is recoverable
exactly where COS matters and not recoverable where it does not. The deliverable
is a diagnostic and a decision rule, not an impossibility result. Acceptable but
do not oversell.

**STRONG (→ MEE)** — if P4 holds and effort bias is not separable from radius
misspecification, then the reported radius cannot be validated from the data at
all in the presence of sampling bias, which is the normal condition for GBIF.
That is an impossibility statement of the same form as m20's bracket, and it is
the only branch that justifies a methods-journal submission.

**ABORT** — if COS bias under realistic misspecification is smaller than the
Monte Carlo SE, the question is not decision-relevant regardless of
identifiability. Report as a robustness section of the benchmark paper.

## 7. What must not be claimed

- Not "misspecification causes bias". That is a tautology; a statistician
  reviewer will say so.
- Not "we discovered COS is sensitive to its assumptions". The contribution has
  to be the identifiability verdict, with a number and a direction.
- Not "location error". The correct term is **Berkson error** — the truth lies
  inside the reported region, so the error is independent of the observed value.
  Berkson misspecification already has a literature (arXiv 2306.01468 and the
  robust-Bayesian line); what is absent is the *spatial* case where the support
  has geometry. Say that precisely or get caught.
- Berkson error is benign for linear models but not for nonlinear ones. SDMs
  have a log link. Cite that, do not assume it.

## 8. Order of work

1. P4 first, on uniform vs structured effort. It is cheap and it decides the
   branch.
2. Only if P4 survives, run the full crossed design.
3. No package code until the branch is fixed.

---

# AMENDMENT 1 — 2026-08-03, after first execution

Predictions in §4 and criteria in §6 are left exactly as written. This section
records what happened and corrects one rule that was wrong on its face.

## A1.1 Structural finding, not anticipated

The published COS likelihood **cannot estimate r at all**. The term
∫_{A_i} λ grows like r², so Σᵢ log(·) increases without bound, while the offset
∫_S λ does not depend on r. The profile is monotone and pins at whatever upper
bound the search grid imposes. The first run returned r̂ = 32.00 (the grid
maximum) in all nine design cells for this reason.

This is not an error in Hefley (2017) or Walker (2020): in both, r is fixed by
metadata and never varied. It does mean the published likelihood is *conditional
on the reporting mechanism* and cannot compare assumed radii.

r becomes a parameter only when the **reporting process** is modelled: the
density of the reported location given the true one, uniform on the disc and
integrating to one. The observed points are then an IPP with intensity
λ ⊛ K̄_r, total mass is preserved, the offset is r-free, and the profile has an
interior maximum. Operationally this is one line — normalise the kernel — but
conceptually it is a different model, and it is a claimable contribution.

## A1.2 P1 / P2 — CONFIRMED

Width of the 2-log-unit interval in r, at γ = 0, n = 15 replicates,
r_true = 8 cells:

| SA range | r_true / SA | r̂ | 2-unit width |
|---|---|---|---|
| 2 | 4.0 | 7.94 | 0.02 |
| 8 | 1.0 | 7.92 | 0.17 |
| 32 | 0.25 | 6.30 | 0.81 |

Identifiability of the radius tracks the Naimi ratio, as predicted in advance.

## A1.3 P4 — NOT CONFIRMED

Where the profile is sharp, effort structure moves nothing:

- SA = 2 : 7.94 → 8.40 → 7.61 across γ = 0, 0.5, 1.0
- SA = 8 : 7.92 → 7.08 → 7.10

The apparent "drift" was confined to SA = 32, where r̂ has SD 5.1–12.4 on a grid
spanning 1.5–32. In that cell the estimator carries almost no information; the
cell mean of a near-uniform variable moved, not the estimate. A 4-replicate
pilot gave 1.50 in the same cell where 15 replicates gave 14.11.

**Effort does not confound the radius where the radius is identifiable.**

## A1.4 The verdict rule in §6 was wrong — corrected

The shipped rule tested drift in r̂ without requiring the cell to be identifiable
and without comparing the shift to its Monte Carlo standard error. It therefore
printed STRONG on the strength of a cell with no information. Corrected rule:

> Drift is evidence of confounding **only** in cells whose 2-log-unit width at
> γ = 0 is below 0.3, and only when |r̂(γ) − r̂(0)| exceeds 3 × MCSE, where
> MCSE = SD(r̂) / √n_rep. Cells with width ≥ 0.3 are reported as
> *uninformative*, never as drift.

Under the corrected rule the current result is **WEAK** (§6): a diagnostic, not
an impossibility. Target Ecography, not MEE, unless something else changes it.

## A1.5 Grid-centroid subset — signature NOT established

Records whose radius equals a round cell side divided by √2 were proposed as
grid centroids (39,033 records; 7.68% of the crayfish download,
DOI 10.15468/dl.99ezk2). The independent check failed: coordinate duplication
was 72.3% in the candidate subset versus 68.9% elsewhere.

The test had no power. Baseline duplication is already ~69% because crayfish are
monitored at fixed localities, so the metric is saturated. The √2 pattern
remains arithmetic coincidence until a **lattice test** — common divisor of
projected coordinate differences, per dataset — passes.

**Do not build that test.** GridDER (Ecological Informatics, 2023) is a
published, cited tool that performs exactly this detection. Use it or cite it;
reimplementing it is scope creep with no publication value.

## A1.6 What the empirical partition actually shows

At 1 km resolution, 507,980 crayfish records:

| state | n | % |
|---|---|---|
| radius inert (< 0.5 cell) — correction not needed | 216,847 | 42.7 |
| no radius reported — state unknown | 145,710 | 28.7 |
| usable for COS | 106,390 | 20.9 |
| candidate grid centroid (unconfirmed) | 39,033 | 7.7 |

At 90 m (Hydrography90m), the COS-relevant share rises to 25.6% of all records
plus 10.9% marginal. Freshwater is a more favourable case than terrestrial.

This partition, not the correction itself, is the strongest empirical result so
far. The question the data supports is *for which records does per-record
uncertainty matter, and why is the answer unknowable for nearly a third of
them* — not *how to propagate uncertainty in GBIF*.

## A1.7 Revised next steps

1. Do not run the full crossed design yet. Under the corrected rule the branch
   is WEAK and the design was sized for a question that has been answered.
2. The open question worth a deciding test is now different: does the
   **reporting-model** correction (normalised kernel, §A1.1) recover β₁ better
   than plain COS when the radius is misspecified? That is the δ_size axis
   alone, and it is cheap.
3. GridDER on the crayfish download, to settle §A1.5 without new code.
