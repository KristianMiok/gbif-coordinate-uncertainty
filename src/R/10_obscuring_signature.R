# =============================================================================
# Is there a detectable obscuring signature in European Orchidaceae?
#
# WHY THIS RUNS BEFORE ANYTHING ELSE
#
# The earlier multi-taxon work (06_, 07_, 08_) measured the three-state
# partition at thresholds of 45-3000 m. Those magnitudes are far below the
# grain of any current biodiversity foundation model (BioAnalyst: 0.25 deg,
# ~25 km), so ordinary positional uncertainty is invisible at that scale and
# cannot drive a calibration audit of such a model.
#
# Deliberate generalization is different. GBIF's sensitive-species guidance
# documents that iNaturalist obscures listed taxa -- orchids named explicitly --
# by replacing the coordinate with a RANDOM POINT inside a 0.2 x 0.2 degree
# cell and setting coordinateUncertaintyInMeters to that cell's DIAGONAL.
#
# Computed from that rule, the expected European signature is:
#   lat 0.2 deg ~ 22.2 km; lon 0.2 deg ~ 22.2 * cos(lat) km
#     40N -> ~27.6 km   50N -> ~26.4 km   60N -> ~24.8 km   70N -> ~23.5 km
# So obscured European records should pile up in roughly 23-28 km, with the
# 20-35 km band capturing the mass. That is ~one BioAnalyst cell: obscuring
# CAN move a record to a neighbouring cell, ordinary GPS error cannot.
#
# THIS SCRIPT DECIDES ONLY ONE THING: does that pile-up exist, and is it big
# enough to condition on. It does not model anything.
#
# KILL CRITERION, fixed before running:
#   if the 20-35 km band holds < 2% of georeferenced European Orchidaceae
#   records that report an uncertainty, OR shows no density excess over the
#   adjacent bands, the obscuring-based design is dead and we fall back to the
#   TabPFN / tabular route. Write the outcome down either way.
#
# NO DOWNLOADS. occ_search(limit = 0) returns counts only.
#
# v2: GBIF returns intermittent 503s on large count queries. All calls now
# retry with exponential backoff, and band counts are cached to disk after
# each query, so an interrupted run resumes instead of restarting.
# =============================================================================

library(rgbif)
options(scipen = 999)          # GBIF rejects "1e+05" in range predicates

stopifnot(packageVersion("rgbif") >= "3.7.0")

dir.create("results", showWarnings = FALSE)
CACHE <- "results/c10_band_cache.rds"

# ---- retry wrapper ---------------------------------------------------------

retry <- function(fn, tries = 7, label = "") {
  for (k in seq_len(tries)) {
    res <- try(fn(), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    wait <- min(2^k, 90)
    message("  GBIF unavailable", if (nzchar(label)) paste0(" [", label, "]") else "",
            " attempt ", k, "/", tries, " -- waiting ", wait, "s")
    Sys.sleep(wait)
  }
  stop("GBIF unreachable after ", tries, " attempts: ", label)
}

# ---- 0. resolve the family key rather than hardcoding it --------------------

orchid <- retry(function() name_backbone(name = "Orchidaceae", rank = "family"),
                label = "name_backbone")
FAMILY_KEY <- orchid$usageKey
message("Orchidaceae usageKey: ", FAMILY_KEY)

gbif_n <- function(..., label = "") {
  args <- list(...)
  res <- retry(function() do.call(occ_search, c(list(
            familyKey = FAMILY_KEY, continent = "EUROPE",
            hasCoordinate = TRUE, limit = 0), args)), label = label)
  Sys.sleep(1)                 # be polite; large counts are expensive server-side
  res$meta$count
}

# ---- 1. denominators -------------------------------------------------------

n_total <- gbif_n(label = "total")
n_unc   <- gbif_n(coordinateUncertaintyInMeters = "0,10000000", label = "any radius")
n_na    <- n_total - n_unc

cat(sprintf(
  "\nEuropean Orchidaceae, georeferenced\n  total            %10d\n  reports radius   %10d  (%.1f%%)\n  no radius (NA)   %10d  (%.1f%%)\n",
  n_total, n_unc, 100 * n_unc / n_total, n_na, 100 * n_na / n_total))

# NOTE: the NA share is itself a result. Crayfish reported a radius for 71.3%;
# Marcer et al. (2022) find 18% across GBIF preserved specimens. Whatever
# orchids do, it is a third data point on that spread, not a nuisance.
#
# CAUTION on reading a HIGH reporting rate as diligence: obscured records are
# ASSIGNED a radius as part of the generalization procedure. Some of the
# reported share may be generalization rather than measurement. Step 2 is what
# separates these, which is why the band scan and not this percentage is the
# decision point.

# ---- 2. the band scan ------------------------------------------------------
# Deliberately asymmetric: fine inside the predicted obscuring window, coarse
# elsewhere. Adjacent bands are the control -- a real signature is a LOCAL
# density excess, not merely a large count in a wide band.

bands <- data.frame(
  lo = c(    0,   100,  1000,  5000, 10000, 15000, 20000, 23000, 26000, 29000, 35000,  50000, 100000),
  hi = c(  100,  1000,  5000, 10000, 15000, 20000, 23000, 26000, 29000, 35000, 50000, 100000, 10000000)
)
bands$n <- NA_real_

if (file.exists(CACHE)) {
  cached <- readRDS(CACHE)
  if (identical(cached$lo, bands$lo) && identical(cached$hi, bands$hi)) {
    bands$n <- cached$n
    message("Resuming from cache: ", sum(!is.na(bands$n)), "/", nrow(bands), " bands done")
  }
}

for (i in seq_len(nrow(bands))) {
  if (!is.na(bands$n[i])) next
  lab <- paste0(bands$lo[i], "-", bands$hi[i], " m")
  message("  querying band ", i, "/", nrow(bands), ": ", lab)
  bands$n[i] <- gbif_n(coordinateUncertaintyInMeters =
                         paste0(bands$lo[i], ",", bands$hi[i]), label = lab)
  saveRDS(bands, CACHE)        # checkpoint after every query
}

bands$pct      <- 100 * bands$n / n_unc
bands$width_km <- (bands$hi - bands$lo) / 1000
bands$density  <- bands$n / bands$width_km          # records per km of band width

cat("\nBand scan (denominator = records reporting a radius)\n")
out <- data.frame(
  band    = sprintf("%7.0f - %9.0f m", bands$lo, bands$hi),
  n       = bands$n,
  pct     = round(bands$pct, 2),
  density = round(bands$density, 1)
)
print(out, row.names = FALSE)

# the decision number
in_win <- bands$lo >= 20000 & bands$hi <= 35000
obscure_window <- sum(bands$n[in_win])
win_pct <- 100 * obscure_window / n_unc

cat(sprintf("\n>>> 20-35 km window: %d records, %.2f%% of radius-reporting records\n",
            obscure_window, win_pct))

# density must stand ABOVE the neighbours, not merely be large
nb_lo  <- bands$density[bands$lo == 15000]
nb_hi  <- bands$density[bands$lo == 35000]
win_d  <- mean(bands$density[in_win])
ratio  <- win_d / mean(c(nb_lo, nb_hi))
cat(sprintf(">>> density in window %.1f vs neighbours %.1f / %.1f  (ratio %.2f)\n",
            win_d, nb_lo, nb_hi, ratio))

passed <- (win_pct >= 2) && (ratio > 1.5)
cat(ifelse(passed,
  ">>> BOTH conditions met. Proceed to step 3.\n",
  ">>> KILL CRITERION TRIGGERED (need >=2% AND density ratio >1.5). Record this and stop.\n"))

# ---- 3. is it publisher-structured? ----------------------------------------
# Only meaningful if step 2 passed. Reuses the logic of 08_: the effect that
# matters is dataset-level, because that is what makes it a covariate shift
# with geographic structure rather than random noise.

if (passed) {

  facet_ds <- function(..., label = "") {
    args <- list(...)
    r <- retry(function() do.call(occ_search, c(list(
           familyKey = FAMILY_KEY, continent = "EUROPE", hasCoordinate = TRUE,
           limit = 0, facet = "datasetKey", facetLimit = 1000), args)),
           label = label)
    Sys.sleep(1)
    d <- r$facets$datasetKey
    names(d) <- c("datasetKey", "n")
    d$n <- as.numeric(d$n)
    d
  }

  ds_all <- facet_ds(label = "facet all")
  ds_win <- facet_ds(coordinateUncertaintyInMeters = "20000,35000", label = "facet window")

  m <- merge(ds_all, ds_win, by = "datasetKey", all.x = TRUE,
             suffixes = c("_all", "_win"))
  m$n_win[is.na(m$n_win)] <- 0
  m$share <- m$n_win / m$n_all

  cat(sprintf("\nDatasets seen: %d (facet cap 1000 -- if hit, raise and re-run)\n",
              nrow(ds_all)))
  cat(sprintf("Datasets with zero window records:   %.1f%%\n",
              100 * mean(m$share == 0)))
  cat(sprintf("Datasets with >50%% window records:  %.1f%%\n",
              100 * mean(m$share > 0.5)))
  cat("\nTop 10 datasets by window share (min 200 records):\n")
  top <- m[m$n_all >= 200, ]
  print(head(top[order(-top$share), c("datasetKey", "n_all", "n_win", "share")], 10),
        row.names = FALSE)

  saveRDS(list(bands = bands, datasets = m, n_total = n_total, n_unc = n_unc,
               win_pct = win_pct, density_ratio = ratio),
          "results/c10_obscuring_signature.rds")
  cat("\nSaved results/c10_obscuring_signature.rds\n")
}

# =============================================================================
# READING THE OUTPUT
#
# The design lives if BOTH hold:
#   (a) the 20-35 km window carries >= 2% of radius-reporting records, AND
#   (b) its density stands clearly above the adjacent bands (ratio > 1.5),
#       because a merely large count in a wide band is not a signature.
#
# Step 3 then says whether the window is concentrated in a few publishers.
# Bimodality there is the point: it means obscuring is a DATASET property with
# geographic footprint, i.e. a spatially structured covariate shift a model
# cannot average away -- which is the thing worth auditing a foundation model
# against.
#
# If the window is diffuse across publishers, it is probably not obscuring but
# genuine coarse georeferencing, and the story weakens to ordinary imprecision.
# Say so and stop; do not rescue it.
#
# CAVEAT TO CARRY FORWARD: the diagonal rule is iNaturalist's. National schemes
# (BSBI, NDFF, national atlases) generalize on their own grids -- 1 km, 10 km,
# hectad -- and will sit in DIFFERENT bands. If step 3 shows a large European
# publisher with a hard uncertainty mode outside 20-35 km, that is a second
# obscuring convention, not a failure. Widen the scan rather than discard it.
# =============================================================================
