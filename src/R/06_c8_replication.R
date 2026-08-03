# =============================================================================
# Does C8 replicate? Is the partition a dataset-level property in general,
# or only in crayfish?
#
# C8 was measured on one download: across 81 crayfish datasets the usable share
# ran 0%-100%, the 10th-90th percentile ran 0%-84.4%, and dataset identity
# explained 52.2% of the variance in whether a record is usable. If that is
# specific to crayfish the central claim of the paper collapses to an anecdote.
#
# NO DOWNLOADS. The GBIF search API facets by datasetKey and accepts
# coordinateUncertaintyInMeters as a range, so per-dataset partition counts come
# from three facet queries per taxon.
#
# The ANOVA R^2 on a binary outcome is computable from counts alone:
#   grand mean  p = sum(usable) / sum(n)
#   SSB = sum_j n_j (p_j - p)^2
#   SST = N p (1 - p)              [binary outcome]
#   R^2 = SSB / SST
# =============================================================================

library(rgbif)

TAXA <- list(
  "crayfish"           = c("Astacidae", "Cambaridae", "Cambaroididae", "Parastacidae"),
  "dragonflies"        = "Odonata",
  "salmonid fishes"    = "Salmonidae",
  "freshwater mussels" = "Unionidae",
  "orchids"            = "Orchidaceae",
  "ground beetles"     = "Carabidae",
  "mosses"             = "Bryophyta",
  "swallowtails"       = "Papilionidae",
  "bats"               = "Chiroptera",
  "amphibians"         = "Ranidae",
  "hard corals"        = "Scleractinia",
  "cetaceans"          = "Cetacea"
)
KINGDOM <- c(orchids = "Plantae", mosses = "Plantae")

RES        <- 1000        # partition resolution, m
FACET_LIM  <- 500         # datasets returned per query
MIN_N      <- 500         # datasets smaller than this are unstable
RMAX       <- 10000000

resolve <- function(nms, kd = "Animalia") {
  vapply(nms, function(n) {
    b <- name_backbone(name = n, kingdom = kd, strict = TRUE)
    if (is.null(b$usageKey) || is.na(b$usageKey))
      stop(sprintf("could not resolve '%s' in %s", n, kd))
    if (!identical(tolower(b$canonicalName), tolower(n)))
      stop(sprintf("'%s' resolved to '%s' -- wrong taxon, refusing", n, b$canonicalName))
    as.character(b$usageKey)
  }, character(1), USE.NAMES = FALSE)
}

facet_counts <- function(keys, ...) {
  s <- occ_search(taxonKey = keys, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                  facet = "datasetKey", facetLimit = FACET_LIM, facetMincount = 1,
                  limit = 0, ...)
  f <- if (!is.null(s$facets)) s$facets$datasetKey else
    do.call(rbind, lapply(s, function(z) z$facets$datasetKey))
  if (is.null(f) || !nrow(f)) return(data.frame(ds = character(), n = numeric()))
  f <- aggregate(list(n = as.numeric(f$count)), by = list(ds = f$name), FUN = sum)
  f[order(-f$n), ]
}

join_on <- function(a, b, col) {
  a[[col]] <- b$n[match(a$ds, b$ds)]
  a[[col]][is.na(a[[col]])] <- 0
  a
}

res <- list()
for (nm in names(TAXA)) {
  kd <- if (nm %in% names(KINGDOM)) KINGDOM[[nm]] else "Animalia"
  keys <- resolve(TAXA[[nm]], kd)

  tot <- facet_counts(keys)
  if (!nrow(tot)) { cat(sprintf("%-20s no facets returned, skipped\n", nm)); next }
  Sys.sleep(1)
  rep_ <- facet_counts(keys, coordinateUncertaintyInMeters = paste0("0,", format(RMAX, scientific = FALSE)))
  Sys.sleep(1)
  hi   <- facet_counts(keys, coordinateUncertaintyInMeters = paste0("0,", RES * 3))
  Sys.sleep(1)

  d <- data.frame(ds = tot$ds, n = tot$n)
  d <- join_on(d, rep_, "n_rep")
  d <- join_on(d, hi,   "n_below_hi")
  d$usable <- (d$n_rep - d$n_below_hi) / d$n
  d <- d[d$n >= MIN_N & is.finite(d$usable), ]
  d$usable <- pmin(pmax(d$usable, 0), 1)
  if (nrow(d) < 5) { cat(sprintf("%-20s only %d usable datasets, skipped\n", nm, nrow(d))); next }

  N  <- sum(d$n)
  p  <- sum(d$usable * d$n) / N
  r2 <- sum(d$n * (d$usable - p)^2) / max(N * p * (1 - p), 1e-12)
  q  <- quantile(d$usable, c(.1, .5, .9))

  res[[nm]] <- list(k = nrow(d), N = N, p = p, r2 = r2, q = q,
                    zero = mean(d$usable < 0.01), full = mean(d$usable > 0.99))
  cat(sprintf("%-20s %3d datasets | overall usable %5.1f%% | R2 %5.1f%%\n",
              nm, nrow(d), 100 * p, 100 * r2))
}

cat("\n=== C8 REPLICATION ACROSS TAXA (partition at 1 km) ===\n")
cat(sprintf("%-20s %5s %11s %8s %8s %8s %8s %8s\n",
            "taxon", "n_ds", "records", "usable", "p10", "p90", "R2", "%ds=0"))
for (nm in names(res)) {
  o <- res[[nm]]
  cat(sprintf("%-20s %5d %11s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n",
              nm, o$k, format(o$N, big.mark = ","), 100 * o$p,
              100 * o$q[[1]], 100 * o$q[[2 + 1]], 100 * o$r2, 100 * o$zero))
}

r2s <- vapply(res, function(o) o$r2, numeric(1))
cat(sprintf("\nvariance explained by dataset identity: %.1f%% to %.1f%% (median %.1f%%)\n",
            100 * min(r2s), 100 * max(r2s), 100 * median(r2s)))
cat(sprintf("crayfish, measured on the full download: 52.2%%\n"))

cat("\n=== REGIMES, NOT GRADIENTS ===\n")
cat(sprintf("%-20s %14s %14s\n", "taxon", "% ds at ~0%", "% ds at ~100%"))
for (nm in names(res)) {
  o <- res[[nm]]
  cat(sprintf("%-20s %13.1f%% %13.1f%%\n", nm, 100 * o$zero, 100 * o$full))
}

saveRDS(res, "results/c8_replication.rds")
cat("\nwrote results/c8_replication.rds\n")

cat("\n--- HOW TO READ THIS ---\n")
cat("C8 REPLICATES if the between-dataset spread (p10 to p90) is wide in most\n")
cat("groups and R2 is substantial throughout. The claim then is general: the\n")
cat("partition is a publisher property and cannot be predicted from the taxon.\n")
cat("C8 FAILS if crayfish is an outlier -- if most groups are homogeneous\n")
cat("across datasets, the central claim reverts to an anecdote and the paper\n")
cat("goes back to Ecography without a package.\n")
cat("\nCAVEAT: facets are capped at FACET_LIM datasets per query, so groups with\n")
cat("many small publishers are truncated. Records below MIN_N are excluded.\n")
cat(sprintf("accessed: %s\n", Sys.Date()))
