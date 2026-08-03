# Data

Nothing here is committed.

## GBIF occurrence download

**DOI 10.15468/dl.99ezk2** — 507,980 records, four crayfish families
(Astacidae, Cambaridae, Cambaroididae, Parastacidae), filtered to
`hasCoordinate = TRUE` and `hasGeospatialIssue = FALSE`. SIMPLE_CSV, ~49 MB
compressed. Download key `0022943-260721160103020`, created 2026-08-03.

Regenerate with `src/R/01_download_and_radius_profile.R`, or fetch by key:

```r
d <- rgbif::occ_download_import(rgbif::occ_download_get("0022943-260721160103020"))
```

GBIF retains prepared downloads for six months. After that, re-request — the
predicate is in the script and the record count may have moved.

Cite the DOI in any manuscript. It is the only reproducibility guarantee here.

## Not representative

Crayfish report `coordinateUncertaintyInMeters` for 71.3% of records. Marcer et
al. (2022, Ecography) report 18% across GBIF preserved specimens. This dataset
is dominated by human observations and is a favourable case. Any partition
computed here must be labelled as such.
