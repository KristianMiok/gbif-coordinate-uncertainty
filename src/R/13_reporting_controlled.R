# =============================================================================
# Does the reporting gap survive controlling for record type and era?
#
# WHERE THIS COMES FROM
#
# C12 tested the `usable` axis (share of records with a radius at or below a
# threshold) and its criterion fired. But the OTHER column separated sharply,
# and along a different mechanism than predicted:
#
#   reported = share of records stating ANY coordinateUncertaintyInMeters
#     Gulo gulo         4.8%      Pieris brassicae   96.9%
#     Ursus arctos     30.9%      Vanessa atalanta   94.6%
#     Aquila fasciata  50.4%      Triturus cristatus 94.1%
#     Lynx lynx        55.4%      Impatiens gland.   88.0%
#
# The prediction was coordinate GENERALIZATION (large radii on sensitive
# species). The data show NON-REPORTING instead: for Gulo gulo, 95% of records
# carry no precision statement at all. That is the "state unknown" class, and
# it matters for a model because BioAnalyst cannot distinguish a precise record
# from one with no stated precision -- it sees neither.
#
# THE CONFOUND THIS SCRIPT EXISTS TO KILL OR CONFIRM
#
# Low reporting may be sensitivity-driven withholding, OR simply a composition
# artefact: large carnivores are documented by older museum and atlas records
# that predate the coordinateUncertaintyInMeters convention, while butterflies
# and bumblebees are documented by modern GPS-carrying citizen scientists.
# Those two explanations make the same marginal prediction and must be
# separated before any claim is made.
#
# Control: restrict to basisOfRecord = HUMAN_OBSERVATION and year >= 2015.
# Same record type, same era, everyone carrying a phone. If the gap persists
# inside that stratum it is not composition.
#
# CRITERION, fixed before running. Deliberately NOT a threshold on prevalence:
# C10 and C12 both failed because they cut an unseen distribution at invented
# points. This tests SURVIVAL of a difference instead.
#
#   BOTH must hold for the reporting axis to be real:
#     (a) controlled spread (max - min reported share, species with >= 1000
#         records in the stratum) >= 40 points
#     (b) Spearman rho between raw and controlled reported share >= 0.5,
#         i.e. the ordering is not rearranged by the control
#
#   controlled spread < 20 points -> composition explains it, axis dies
#   spread 20-40, or rho < 0.5    -> attenuated; report as partially confounded
#                                     and do not build a design on it
#
# NO DOWNLOADS. Cumulative queries only.
# =============================================================================

library(rgbif)
options(scipen = 999)

dir.create("results", showWarnings = FALSE)
CACHE <- "results/c13_reporting_cache.rds"

prev <- readRDS("results/c12_bfm_species_quality.rds")
base <- prev$profile[, c("key", "name", "n_total", "n_radius", "pct_reported")]

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

q <- function(..., label = "") {
  args <- list(...)
  res <- retry(function() do.call(occ_search, c(list(
           continent = "EUROPE", hasCoordinate = TRUE, limit = 0), args)),
           label = label)
  Sys.sleep(2)
  res
}

cached  <- if (file.exists(CACHE)) readRDS(CACHE) else list(done = list(), partial = list())
res     <- cached$done
partial <- cached$partial
if (is.null(res)) res <- list()
if (is.null(partial)) partial <- list()

# ---- per-species: controlled counts + record-type composition --------------

for (i in seq_len(nrow(base))) {
  k  <- base$key[i]
  kk <- as.character(k)
  if (!is.null(res[[kk]])) next

  message(sprintf("  [%2d/%d] %s", i, nrow(base), base$name[i]))
  row <- if (!is.null(partial[[kk]])) partial[[kk]] else
           list(key = k, name = base$name[i])

  put <- function(field, value_fn) {
    if (!is.null(row[[field]])) return(invisible(NULL))
    row[[field]] <<- value_fn()
    partial[[kk]] <<- row
    saveRDS(list(done = res, partial = partial), CACHE)
  }

  # (1) controlled stratum: same record type, same era
  put("n_ctrl", function()
    q(taxonKey = k, basisOfRecord = "HUMAN_OBSERVATION", year = "2015,2026",
      label = paste("ctrl n", k))$meta$count)

  put("n_ctrl_radius", function() {
    if (row$n_ctrl == 0) return(0)
    q(taxonKey = k, basisOfRecord = "HUMAN_OBSERVATION", year = "2015,2026",
      coordinateUncertaintyInMeters = "0,10000000",
      label = paste("ctrl radius", k))$meta$count
  })

  # (2) record-type composition, to document that the confound is real
  put("bor", function() {
    f <- q(taxonKey = k, facet = "basisOfRecord", facetLimit = 20,
           label = paste("bor", k))$facets$basisOfRecord
    if (is.null(f)) return(list())
    setNames(as.numeric(f[[2]]), f[[1]])
  })

  # (3) modern share, to show era composition separately from type
  put("n_recent", function()
    q(taxonKey = k, year = "2015,2026", label = paste("recent", k))$meta$count)

  res[[kk]] <- row
  partial[[kk]] <- NULL
  saveRDS(list(done = res, partial = partial), CACHE)
}

# ---- assemble --------------------------------------------------------------

d <- do.call(rbind, lapply(res, function(r) data.frame(
  key = r$key, name = r$name,
  n_ctrl = r$n_ctrl, n_ctrl_radius = r$n_ctrl_radius,
  n_recent = r$n_recent,
  pct_human = {
    b <- r$bor
    if (length(b) == 0) NA else 100 * sum(b[names(b) == "HUMAN_OBSERVATION"]) / sum(b)
  },
  stringsAsFactors = FALSE)))

d <- merge(base, d, by = c("key", "name"))
d$pct_recent   <- 100 * d$n_recent / pmax(d$n_total, 1)
d$ctrl_reported <- 100 * d$n_ctrl_radius / pmax(d$n_ctrl, 1)

# ---- 1. is the confound real? ----------------------------------------------
# It has to be shown to exist before it is shown not to explain the gap.

cat("\n=== Is the composition confound real? ===\n")
cat(sprintf("HUMAN_OBSERVATION share: %.1f%% to %.1f%% across species\n",
            min(d$pct_human, na.rm = TRUE), max(d$pct_human, na.rm = TRUE)))
cat(sprintf("Post-2015 share:         %.1f%% to %.1f%%\n",
            min(d$pct_recent), max(d$pct_recent)))
cc1 <- cor(d$pct_reported, d$pct_human, method = "spearman", use = "complete.obs")
cc2 <- cor(d$pct_reported, d$pct_recent, method = "spearman")
cat(sprintf("Spearman(reported, %%human)  = %+.2f\n", cc1))
cat(sprintf("Spearman(reported, %%recent) = %+.2f\n", cc2))
cat("A strong positive correlation here means the confound is real and the\ncontrol below is necessary.\n")

# ---- 2. does the gap survive the control? ----------------------------------

elig <- d[d$n_ctrl >= 1000, ]
cat(sprintf("\n=== Controlled stratum (HUMAN_OBSERVATION, year >= 2015) ===\n"))
cat(sprintf("Species with >= 1,000 records in stratum: %d / %d\n", nrow(elig), nrow(d)))

show <- data.frame(
  name        = substr(elig$name, 1, 26),
  n_all       = elig$n_total,
  raw_rep     = round(elig$pct_reported, 1),
  n_ctrl      = elig$n_ctrl,
  ctrl_rep    = round(elig$ctrl_reported, 1),
  delta       = round(elig$ctrl_reported - elig$pct_reported, 1)
)
print(show[order(show$ctrl_rep), ], row.names = FALSE)

spread_raw  <- diff(range(elig$pct_reported))
spread_ctrl <- diff(range(elig$ctrl_reported))
rho <- cor(elig$pct_reported, elig$ctrl_reported, method = "spearman")

cat(sprintf("\nSpread raw:        %.1f points\n", spread_raw))
cat(sprintf("Spread controlled: %.1f points\n", spread_ctrl))
cat(sprintf("Spearman(raw, controlled) = %+.2f\n", rho))

verdict <- if (spread_ctrl >= 40 && rho >= 0.5) {
  "SURVIVES -- the reporting gap is not a composition artefact"
} else if (spread_ctrl < 20) {
  "DIES -- composition explains the gap; the reporting axis is an artefact"
} else {
  "ATTENUATED -- partially confounded; report, do not build a design on it"
}
cat("\n>>> ", verdict, "\n", sep = "")

saveRDS(list(data = d, eligible = elig, spread_raw = spread_raw,
             spread_ctrl = spread_ctrl, rho = rho, verdict = verdict,
             cor_human = cc1, cor_recent = cc2),
        "results/c13_reporting_controlled.rds")
cat("Saved results/c13_reporting_controlled.rds\n")

# =============================================================================
# READING THE OUTPUT
#
# SURVIVES: the claim is that occurrence records for sensitive species carry no
# precision statement at a rate unexplained by when or how they were collected,
# and that a model trained on them cannot see this. Next step is spatial: which
# 0.25 deg cells of BioAnalyst's 160 x 280 grid are dominated by non-reporting
# records. That is the audit stratum, and it needs a real occ_download.
#
# DIES: honour it. The reporting axis is composition, the whole BioAnalyst
# stratification idea is exhausted, and the honest write-up is that data
# quality varies across the model's training species but not in a way that is
# separable from how those species happen to be surveyed. That is still a
# publishable paragraph and it closes the question.
#
# ATTENUATED: the most likely outcome and the least comfortable. Report the
# attenuation, do not proceed to the download, and do not go looking for a
# fourth axis. Three designs have now been tested against pre-registered
# criteria; a fourth attempt after three failures is fishing, not screening.
#
# CAVEAT: year >= 2015 and HUMAN_OBSERVATION are themselves choices. If the
# verdict lands near a boundary, the honest move is to state that the result
# depends on the control window, not to search for a window that passes.
# =============================================================================
