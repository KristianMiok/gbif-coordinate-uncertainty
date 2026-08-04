# =============================================================================
# Do the 28 species BioAnalyst was trained on differ in occurrence-data quality?
#
# WHAT THIS IS
#
# BioAnalyst (bfm-model, Trantas et al. 2025) lists its target species in
# bfm_model/bfm/configs/train_config.yaml as `species_vars`. Those strings are
# GBIF taxonKeys. The model's training targets are therefore directly
# auditable in GBIF, which removes the need to pick a taxon by hand.
#
# The earlier work (08_) established that occurrence quality is a DATASET
# property, not a biological one: knowing the publisher removes a median 59%
# of the uncertainty about whether a record is usable (U = 0.593, permutation
# null 0.000). If that heterogeneity also separates the 28 model species, it
# gives a grouping variable for a group-conditional calibration audit of the
# model -- the same move as language in the BAN hate-speech work, where BERT's
# own scores were informative for English (p = 1.4e-11) and worthless for
# Croatian (p = 1).
#
# WHAT THIS IS NOT
#
# This is a SCREEN, not the stratum itself. The stratification an audit
# actually needs is spatial -- grid cells whose training records came from
# good vs. poor publishers -- and that requires a real occ_download with
# coordinates. These count queries only decide whether that download is worth
# doing. Passing here does not mean the audit is ready to run.
#
# KILL CRITERION, fixed before running:
#   >= 8 species below 60% usable AND >= 8 species above 85% usable,
#   each with >= 5,000 European records.
#   -> two comparable groups exist; proceed to the download
#   If all 28 fall inside a 15-point band, or either group has < 8 qualifying
#   species, there is nothing to condition on and the design dies here.
#
# "Usable" is defined as in the earlier scripts: a record reports
# coordinateUncertaintyInMeters at or below the working threshold. The
# threshold is swept rather than fixed, because the partition should not be an
# artefact of one cut.
#
# NOTE ON RANGE PREDICATES: GBIF ranges are inclusive at BOTH ends, so adjacent
# disjoint bands double-count boundary values (this bit 10_, which summed to
# 132% of its denominator). Everything here uses CUMULATIVE queries from zero,
# which are unaffected, matching what 04_ and 06_ did.
#
# NO DOWNLOADS.
# =============================================================================

library(rgbif)
options(scipen = 999)

dir.create("results", showWarnings = FALSE)
CACHE <- "results/c12_bfm_species_cache.rds"

# ---- the 28 keys, verbatim from train_config.yaml --------------------------

BFM_KEYS <- c(
  1340361, 1340503, 1536449, 1898286, 1920506, 2430567, 2431885, 2433433,
  2434779, 2435240, 2435261, 2437394, 2441454, 2473958, 2491534, 2891770,
  3034825, 4408498, 5218786, 5219073, 5219173, 5219219, 5844449, 8002952,
  8077224, 8894817, 8909809, 9809229
)
stopifnot(length(BFM_KEYS) == 28)

THRESHOLDS <- c(1000, 5000, 10000)   # swept; 1000 m is the conventional cut

retry <- function(fn, tries = 10, label = "") {
  for (k in seq_len(tries)) {
    res <- try(fn(), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    wait <- min(2^k, 300)
    message("    GBIF unavailable [", label, "] ", k, "/", tries, " -- ", wait, "s")
    Sys.sleep(wait)
  }
  stop("GBIF unreachable: ", label)
}

n_for <- function(key, ..., label = "") {
  args <- list(...)
  res <- retry(function() do.call(occ_search, c(list(
           taxonKey = key, continent = "EUROPE", hasCoordinate = TRUE,
           limit = 0), args)), label = label)
  Sys.sleep(2.0)
  res$meta$count
}

# ---- 1. resolve keys to names ----------------------------------------------
# Printed so the species list can be sanity-checked against the model paper.

cat("Resolving 28 taxon keys...\n")
info <- do.call(rbind, lapply(BFM_KEYS, function(k) {
  u <- retry(function() name_usage(key = k)$data, label = paste("name", k))
  data.frame(key = k,
             name = if (!is.null(u$canonicalName)) u$canonicalName else u$scientificName,
             rank = u$rank,
             kingdom = if (!is.null(u$kingdom)) u$kingdom else NA,
             class = if (!is.null(u$class)) u$class else NA,
             stringsAsFactors = FALSE)
}))
print(info[, c("key", "name", "rank", "class")], row.names = FALSE)

# ---- 2. quality profile per species ----------------------------------------

cached  <- if (file.exists(CACHE)) readRDS(CACHE) else list()
if (!is.null(cached$done)) { res <- cached$done; partial <- cached$partial } else {
  res <- cached; partial <- list() }   # migrate old per-species cache
if (is.null(partial)) partial <- list()

for (i in seq_len(nrow(info))) {
  k <- info$key[i]
  kk <- as.character(k)
  if (!is.null(res[[kk]])) next

  message(sprintf("  [%2d/28] %s (%d)", i, info$name[i], k))
  row <- if (!is.null(partial[[kk]])) partial[[kk]] else
           list(key = k, name = info$name[i])

  # checkpoint after EVERY query, not every species: GBIF throttles after a few
  # hundred calls, and losing a whole species row to one refusal wastes the
  # queries that already succeeded.
  put <- function(field, ...) {
    if (!is.null(row[[field]])) return(row)
    row[[field]] <<- n_for(...)
    partial[[kk]] <<- row
    saveRDS(list(done = res, partial = partial), CACHE)
    row
  }

  put("n_total", k, label = paste("total", k))
  if (row$n_total == 0) {
    row$n_radius <- 0
    for (t in THRESHOLDS) row[[paste0("n_le_", t)]] <- 0
  } else {
    put("n_radius", k, coordinateUncertaintyInMeters = "0,10000000",
        label = paste("radius", k))
    for (t in THRESHOLDS) {
      put(paste0("n_le_", t), k,
          coordinateUncertaintyInMeters = paste0("0,", t),
          label = paste0("le", t, " ", k))
    }
  }
  res[[kk]] <- row
  partial[[kk]] <- NULL
  saveRDS(list(done = res, partial = partial), CACHE)
}

d <- do.call(rbind, lapply(res, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
d <- d[order(match(d$key, BFM_KEYS)), ]

# ---- 3. the numbers --------------------------------------------------------
# Two distinct quantities, and they must not be conflated:
#   reported  = share of records that state a radius at all
#   usable    = share of ALL records with a radius at or below the threshold
# A record with no radius is NOT usable for a correction but is also not known
# to be bad -- that is the "state unknown" class from the earlier work.

d$pct_reported <- 100 * d$n_radius / pmax(d$n_total, 1)
for (t in THRESHOLDS) {
  d[[paste0("usable_", t)]] <- 100 * d[[paste0("n_le_", t)]] / pmax(d$n_total, 1)
}

cat("\n=== Quality profile of BioAnalyst's 28 training species (Europe) ===\n")
show <- data.frame(
  name      = substr(d$name, 1, 26),
  n         = d$n_total,
  reported  = round(d$pct_reported, 1),
  u1km      = round(d$usable_1000, 1),
  u5km      = round(d$usable_5000, 1),
  u10km     = round(d$usable_10000, 1)
)
print(show[order(-show$u1km), ], row.names = FALSE)

# ---- 4. the decision -------------------------------------------------------

eligible <- d[d$n_total >= 5000, ]
cat(sprintf("\nSpecies with >= 5,000 European records: %d / 28\n", nrow(eligible)))

verdict_at <- function(t) {
  col <- paste0("usable_", t)
  lo  <- sum(eligible[[col]] <  60)
  hi  <- sum(eligible[[col]] >  85)
  rng <- if (nrow(eligible)) diff(range(eligible[[col]])) else 0
  cat(sprintf("  threshold %5d m:  low(<60%%) = %2d   high(>85%%) = %2d   spread = %.1f pts\n",
              t, lo, hi, rng))
  c(lo = lo, hi = hi, spread = rng)
}

cat("\nGrouping at each threshold (eligible species only):\n")
v <- sapply(THRESHOLDS, verdict_at)
colnames(v) <- THRESHOLDS

passes <- apply(v, 2, function(x) x["lo"] >= 8 && x["hi"] >= 8 && x["spread"] >= 15)

cat("\n")
if (any(passes)) {
  ok <- THRESHOLDS[passes]
  cat(sprintf(">>> PASS at threshold(s): %s\n", paste(ok, "m", collapse = ", ")))
  cat(">>> Two comparable groups exist. The spatial download is worth doing.\n")
  cat(">>> REMINDER: this is a screen. The audit stratum is spatial and still\n")
  cat(">>> requires occ_download with coordinates for these species.\n")
  t_use <- ok[1]
  col <- paste0("usable_", t_use)
  cat(sprintf("\nAt %d m -- LOW group (<60%% usable):\n", t_use))
  print(eligible[eligible[[col]] < 60, c("name", "n_total", col)], row.names = FALSE)
  cat(sprintf("\nAt %d m -- HIGH group (>85%% usable):\n", t_use))
  print(eligible[eligible[[col]] > 85, c("name", "n_total", col)], row.names = FALSE)
} else {
  cat(">>> KILL CRITERION TRIGGERED at every threshold.\n")
  cat(">>> The 28 species do not separate on data quality. There is nothing to\n")
  cat(">>> condition on. Record this and stop; do not loosen the thresholds.\n")
}

saveRDS(list(info = info, profile = d, verdict = v, passes = passes),
        "results/c12_bfm_species_quality.rds")
cat("\nSaved results/c12_bfm_species_quality.rds\n")

# =============================================================================
# READING THE OUTPUT
#
# PASS: proceed to an occ_download for these 28 keys, Europe, keeping
# datasetKey and coordinateUncertaintyInMeters. Then build the real stratum:
# for each 0.25 deg cell on BioAnalyst's 160 x 280 grid, the share of its
# training records that are usable. That per-cell share -- not the per-species
# share screened here -- is what a group-conditional coverage test conditions on.
#
# FAIL: honour it. The model's species were chosen for data availability, so
# uniformly good coverage is a plausible outcome and would be an honest
# negative result worth a paragraph. Do not lower the bar to rescue it.
#
# EITHER WAY, one number is already publishable: the spread of `reported` and
# `usable` across the species a foundation model was trained on. Whatever it
# shows, no one has reported it, and it is the empirical hook for the claim
# that training-data quality is uneven in a way the model cannot see.
#
# CAVEAT: `continent = "EUROPE"` is GBIF's own continent assignment and is not
# identical to BioAnalyst's grid (H=160, W=280 at 0.25 deg). The screen is
# approximate by design; the download step must use an explicit bounding box
# matched to the model grid.
# =============================================================================
