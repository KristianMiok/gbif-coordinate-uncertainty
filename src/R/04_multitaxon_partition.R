# =============================================================================
# Does the partition vary across taxa?
#
# The central empirical claim so far -- that occurrence records split into
# "correction needed", "correction inert", and "state unknown" -- was measured
# on one group. Crayfish report a radius for 71.3% of records; Marcer et al.
# (2022) find 18% across GBIF preserved specimens. The crayfish figure is
# therefore not transferable, and the partition may be a property of the taxon
# rather than of GBIF.
#
# This decides whether a diagnostic package has a reason to exist:
#   - if the partition is broadly similar across groups, a table in the paper
#     is enough and no software is warranted
#   - if it varies widely and cannot be guessed from the taxon, users need to
#     compute it for their own data, which is what a package is for
#
# NO DOWNLOADS. The GBIF search API accepts coordinateUncertaintyInMeters as a
# range, so cumulative counts at the thresholds that define the partition can be
# queried directly. ~13 calls per taxon.
# =============================================================================

library(rgbif)

# groups chosen to span the axes that plausibly drive the partition:
# realm (freshwater / terrestrial / marine), record type (specimen-dominated vs
# observation-dominated), and georeferencing tradition.
TAXA <- list(
  # freshwater
  "crayfish"          = c("Astacidae", "Cambaridae", "Cambaroididae", "Parastacidae"),
  "dragonflies"       = "Odonata",
  "salmonid fishes"   = "Salmonidae",
  "freshwater mussels"= "Unionidae",
  # terrestrial, specimen-heavy
  "orchids"           = "Orchidaceae",
  "ground beetles"    = "Carabidae",
  "mosses"            = "Bryophyta",
  # terrestrial, observation-heavy
  "swallowtails"      = "Papilionidae",
  "bats"              = "Chiroptera",
  "amphibians"        = "Ranidae",
  # marine
  "hard corals"       = "Scleractinia",
  "cetaceans"         = "Cetacea",
  # accidental but valuable: no taxonKey = the whole of GBIF, as a baseline
  "ALL GBIF"          = NA
)

# partition thresholds: for resolution R, inert if r < R/2, usable if r > 3R
RES <- c(90, 250, 1000)
CUTS <- sort(unique(c(RES / 2, RES * 3, 15, 100, 1000, 10000, 100000)))
RMAX <- 10000000

cnt <- function(keys, ...) {
  s <- occ_search(taxonKey = keys, hasCoordinate = TRUE,
                  hasGeospatialIssue = FALSE, limit = 0, ...)
  if (is.null(s$meta)) sum(vapply(s, function(z) z$meta$count, numeric(1)))
  else s$meta$count
}

KINGDOM <- c(orchids = "Plantae", mosses = "Plantae")   # default Animalia

resolve <- function(nms, kd = "Animalia") {
  if (length(nms) == 1 && is.na(nms)) return(NULL)   # deliberate: all of GBIF
  k <- vapply(nms, function(n) {
    b <- name_backbone(name = n, kingdom = kd, strict = TRUE)
    if (is.null(b$usageKey) || is.na(b$usageKey))
      stop(sprintf("could not resolve '%s' in %s", n, kd))
    if (!identical(tolower(b$canonicalName), tolower(n)))
      stop(sprintf("'%s' resolved to '%s' (rank %s) -- wrong taxon, refusing",
                   n, b$canonicalName, b$rank))
    if (!b$rank %in% c("FAMILY", "ORDER", "CLASS", "PHYLUM", "SUPERFAMILY"))
      stop(sprintf("'%s' resolved at rank %s -- too high or too low", n, b$rank))
    as.character(b$usageKey)
  }, character(1))
  unname(k)
}

out <- list()
for (nm in names(TAXA)) {
  keys <- resolve(TAXA[[nm]], if (nm %in% names(KINGDOM)) KINGDOM[[nm]] else "Animalia")
  n_geo <- cnt(keys)
  n_rep <- cnt(keys, coordinateUncertaintyInMeters = paste0("0,", format(RMAX, scientific = FALSE)))
  cum <- vapply(CUTS, function(u)
    cnt(keys, coordinateUncertaintyInMeters = paste0("0,", format(u, scientific = FALSE))), numeric(1))
  names(cum) <- CUTS
  out[[nm]] <- list(n_geo = n_geo, n_rep = n_rep, cum = cum)
  cat(sprintf("%-20s georeferenced %10s | radius reported %6.1f%%\n",
              nm, format(n_geo, big.mark = ","), 100 * n_rep / n_geo))
  Sys.sleep(1.0)
}

# --- reporting rate ----------------------------------------------------------
cat("\n=== REPORTING RATE ===\n")
cat(sprintf("%-20s %12s %12s %8s\n", "taxon", "georef", "with radius", "%"))
for (nm in names(out)) {
  o <- out[[nm]]
  cat(sprintf("%-20s %12s %12s %7.1f%%\n", nm,
              format(o$n_geo, big.mark = ","),
              format(o$n_rep, big.mark = ","), 100 * o$n_rep / o$n_geo))
}
rr <- vapply(out, function(o) 100 * o$n_rep / o$n_geo, numeric(1))
cat(sprintf("\nrange across taxa: %.1f%% to %.1f%%  (spread %.1f points)\n",
            min(rr), max(rr), max(rr) - min(rr)))

# --- partition, per resolution ----------------------------------------------
for (res in RES) {
  lo <- as.character(res / 2); hi <- as.character(res * 3)
  cat(sprintf("\n=== PARTITION AT %d m (share of ALL georeferenced records) ===\n", res))
  cat(sprintf("%-20s %8s %10s %8s %10s\n",
              "taxon", "inert", "marginal", "usable", "no radius"))
  for (nm in names(out)) {
    o <- out[[nm]]
    inert <- o$cum[[lo]]
    marg  <- o$cum[[hi]] - o$cum[[lo]]
    usable <- o$n_rep - o$cum[[hi]]
    nona  <- o$n_geo - o$n_rep
    cat(sprintf("%-20s %7.1f%% %9.1f%% %7.1f%% %9.1f%%\n", nm,
                100 * inert / o$n_geo, 100 * marg / o$n_geo,
                100 * usable / o$n_geo, 100 * nona / o$n_geo))
  }
  us <- vapply(out, function(o)
    100 * (o$n_rep - o$cum[[hi]]) / o$n_geo, numeric(1))
  cat(sprintf("\nusable share ranges %.1f%% to %.1f%% (%.0fx)  min: %s  max: %s\n",
              min(us), max(us), max(us) / max(min(us), 0.01),
              names(which.min(us)), names(which.max(us))))
}

# --- radius distribution shape ----------------------------------------------
cat("\n=== CUMULATIVE %% OF REPORTING RECORDS BELOW EACH RADIUS ===\n")
show <- as.character(c(15, 45, 100, 500, 1000, 3000, 10000, 100000))
cat(sprintf("%-20s", "taxon")); cat(sprintf("%9s", paste0(show, "m")), sep = ""); cat("\n")
for (nm in names(out)) {
  o <- out[[nm]]
  cat(sprintf("%-20s", nm))
  cat(sprintf("%8.1f%%", 100 * o$cum[show] / o$n_rep), sep = "")
  cat("\n")
}

saveRDS(out, "results/multitaxon_counts.rds")
cat("\nwrote results/multitaxon_counts.rds\n")

cat("\n--- HOW TO READ THIS ---\n")
cat("If the usable share is similar everywhere, a table in the paper suffices\n")
cat("and a diagnostic package has no reason to exist. If it varies by an order\n")
cat("of magnitude and does not track realm or record type in a guessable way,\n")
cat("users cannot predict it for their own data, and that is the package.\n")
cat("\nCounts are live GBIF queries and will drift. Record the access date.\n")
cat(sprintf("accessed: %s\n", Sys.Date()))
