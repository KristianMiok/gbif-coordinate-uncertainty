# =============================================================================
# Is the 26-29 km spike a POINT MASS or a spread?
#
# 10_ found a sharp excess in European Orchidaceae at 26-29 km:
#   20-23 km:    1,573
#   23-26 km:    6,877
#   26-29 km:   98,405   <-- density 32,802/km vs ~500-2,300 either side
#   29-35 km:    1,499
# Density ratio to the outer neighbours: 10.65. This is bimodality, not a tail.
#
# The pre-registered 2% prevalence gate did NOT pass (1.75%). That is recorded
# in logs/10_. This script tests the quantity that actually decides the design,
# with its own criterion fixed BEFORE running:
#
#   >= 50% of the 26-29 km mass in a single 10 m bin
#        -> machine-assigned value; the partition is exact, no threshold needed;
#           the obscuring design lives
#   largest 10 m bin < 20%
#        -> coarse georeferencing on some grid; grey zone; design dies,
#           fall back to the TabPFN / tabular route
#   between 20% and 50%
#        -> several conventions at once; record and decide separately
#
# Rationale for looking for a point mass: iNaturalist assigns the diagonal of a
# 0.2 x 0.2 deg cell. If computed per record it varies with latitude and spreads;
# if a fixed nominal value is used it concentrates. A concentrated value is a
# procedural fingerprint and gives a partition with no arbitrary cut.
#
# KNOWN ISSUE INHERITED FROM 10_: GBIF range predicates are inclusive at BOTH
# ends, so adjacent bands double-count records sitting exactly on a boundary
# (10_ summed to 132% of its denominator). Here the bins are narrow and the
# boundaries are round numbers unlikely to be modes, but the same caveat holds:
# treat bin counts as approximate near boundaries. The winner-take-all reading
# below is robust to this; the exact percentage is not.
#
# NO DOWNLOADS.
# =============================================================================

library(rgbif)
options(scipen = 999)

dir.create("results", showWarnings = FALSE)
CACHE <- "results/c11_mode_cache.rds"

retry <- function(fn, tries = 7, label = "") {
  for (k in seq_len(tries)) {
    res <- try(fn(), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    wait <- min(2^k, 90)
    message("  GBIF unavailable [", label, "] ", k, "/", tries, " -- waiting ", wait, "s")
    Sys.sleep(wait)
  }
  stop("GBIF unreachable: ", label)
}

FAMILY_KEY <- 7689   # Orchidaceae, confirmed in 10_

gbif_n <- function(..., label = "") {
  args <- list(...)
  res <- retry(function() do.call(occ_search, c(list(
           familyKey = FAMILY_KEY, continent = "EUROPE",
           hasCoordinate = TRUE, limit = 0), args)), label = label)
  Sys.sleep(1)
  res$meta$count
}

scan_bins <- function(lo, hi, step, cache_key) {
  edges <- seq(lo, hi, by = step)
  bins  <- data.frame(lo = head(edges, -1), hi = tail(edges, -1), n = NA_real_)

  cache <- if (file.exists(CACHE)) readRDS(CACHE) else list()
  if (!is.null(cache[[cache_key]]) &&
      identical(cache[[cache_key]]$lo, bins$lo)) {
    bins$n <- cache[[cache_key]]$n
    message("Resuming ", cache_key, ": ", sum(!is.na(bins$n)), "/", nrow(bins), " bins")
  }

  for (i in seq_len(nrow(bins))) {
    if (!is.na(bins$n[i])) next
    lab <- paste0(bins$lo[i], "-", bins$hi[i])
    message("  ", cache_key, " bin ", i, "/", nrow(bins), ": ", lab, " m")
    bins$n[i] <- gbif_n(coordinateUncertaintyInMeters =
                          paste0(bins$lo[i], ",", bins$hi[i]), label = lab)
    cache[[cache_key]] <- bins
    saveRDS(cache, CACHE)
  }
  bins
}

# ---- round 1: 100 m bins across the spike ----------------------------------

cat("\n=== Round 1: 100 m bins, 26,000-29,000 m ===\n")
r1 <- scan_bins(26000, 29000, 100, "r1_100m")
r1$pct <- 100 * r1$n / sum(r1$n)

top1 <- r1[order(-r1$n), ][1:8, ]
cat("\nTop 8 bins of 100 m:\n")
print(data.frame(bin = sprintf("%d - %d m", top1$lo, top1$hi),
                 n = top1$n, pct = round(top1$pct, 1)), row.names = FALSE)
cat(sprintf("\nTotal in 26-29 km scan: %d\n", sum(r1$n)))
cat(sprintf("Largest 100 m bin holds %.1f%% of the spike\n", max(r1$pct)))

# ---- round 2: 10 m bins inside the winning 100 m bin ------------------------

best <- r1[which.max(r1$n), ]
cat(sprintf("\n=== Round 2: 10 m bins, %d-%d m ===\n", best$lo, best$hi))
r2 <- scan_bins(best$lo, best$hi, 10, "r2_10m")
r2$pct_of_spike <- 100 * r2$n / sum(r1$n)

top2 <- r2[order(-r2$n), ][1:5, ]
cat("\nTop 5 bins of 10 m (pct of the whole 26-29 km spike):\n")
print(data.frame(bin = sprintf("%d - %d m", top2$lo, top2$hi),
                 n = top2$n, pct = round(top2$pct_of_spike, 1)), row.names = FALSE)

share <- max(r2$pct_of_spike)
cat(sprintf("\n>>> Largest 10 m bin: %d-%d m, %d records, %.1f%% of the spike\n",
            top2$lo[1], top2$hi[1], top2$n[1], share))

verdict <- if (share >= 50) "POINT MASS -- design lives, partition is exact" else
           if (share <  20) "SPREAD -- coarse georeferencing; design dies, go to TabPFN route" else
                            "MIXED (20-50%) -- several conventions; record and decide separately"
cat(">>> ", verdict, "\n", sep = "")

# ---- publisher structure of the mode ---------------------------------------
# Run whenever a mode is identifiable at all (>= 20%), because "which publisher"
# is what turns this from a curiosity into a spatially structured covariate
# shift. If one publisher owns the mode, the degradation has a geographic
# footprint a model cannot average away.

if (share >= 20) {
  cat("\n=== Publisher structure of the mode ===\n")
  r <- retry(function() occ_search(
        familyKey = FAMILY_KEY, continent = "EUROPE", hasCoordinate = TRUE,
        coordinateUncertaintyInMeters = paste0(top2$lo[1], ",", top2$hi[1]),
        limit = 0, facet = "datasetKey", facetLimit = 50), label = "facet mode")
  d <- r$facets$datasetKey
  names(d) <- c("datasetKey", "n")
  d$n <- as.numeric(d$n)
  d$pct <- 100 * d$n / sum(d$n)
  cat("\nTop datasets holding the mode:\n")
  print(head(d[order(-d$n), ], 10), row.names = FALSE)
  cat(sprintf("\nTop dataset holds %.1f%% of mode records\n", max(d$pct)))

  saveRDS(list(r1 = r1, r2 = r2, mode_share = share, mode_bin = top2[1, ],
               datasets = d, verdict = verdict),
          "results/c11_mode.rds")
  cat("Saved results/c11_mode.rds\n")
}

# =============================================================================
# WHAT TO DO WITH EACH VERDICT
#
# POINT MASS: the obscured group is defined by an exact value, not a threshold.
# That is the cleanest possible partition for a group-conditional calibration
# audit -- no grey zone, no cut to defend in review. Next step is a real
# occ_download restricted to European Orchidaceae, keeping the mode value as a
# label, and the audit proceeds on that partition.
#
# SPREAD: honour it. Log the negative result, note that the 26-29 km excess is
# grid-based coarse georeferencing rather than machine obscuring, and move to
# the TabPFN route. Do not widen bins to rescue a mode.
#
# MIXED: likely iNaturalist plus one or more national schemes on their own
# grids. This is not a failure but it does mean the partition is multi-valued;
# decide whether to audit against the largest single convention or against
# "any machine-assigned value", and say which in the manuscript.
# =============================================================================
