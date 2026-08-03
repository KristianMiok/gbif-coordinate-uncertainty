# =============================================================================
# How much of the partition is really a dataset property?
#
# 06_ reported a one-way ANOVA R^2 on a binary outcome with up to 500 groups.
# That spends 500 degrees of freedom and is optimistic by construction. This
# replaces it with the intraclass correlation from a binomial mixed model with
# dataset as a random intercept, which is the quantity actually being claimed:
#
#   cbind(usable, n - usable) ~ 1 + (1 | dataset)
#   ICC = sigma2_u / (sigma2_u + pi^2 / 3)        [latent-scale, logistic]
#
# Also raises the facet cap for the four groups that hit it at 500, to show
# whether the estimate is stable when smaller publishers enter.
#
# NO DOWNLOADS.
# =============================================================================

for (p in c("rgbif", "lme4")) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
library(rgbif); library(lme4)

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
KINGDOM   <- c(orchids = "Plantae", mosses = "Plantae")
RES       <- 1000
FACET_LIM <- 1000        # raised from 500
MIN_N     <- 500
RMAX      <- 10000000

resolve <- function(nms, kd = "Animalia") {
  vapply(nms, function(n) {
    b <- name_backbone(name = n, kingdom = kd, strict = TRUE)
    if (is.null(b$usageKey) || is.na(b$usageKey)) stop(sprintf("unresolved: %s", n))
    if (!identical(tolower(b$canonicalName), tolower(n)))
      stop(sprintf("'%s' resolved to '%s' -- refusing", n, b$canonicalName))
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
  aggregate(list(n = as.numeric(f$count)), by = list(ds = f$name), FUN = sum)
}

out <- list(); tabs <- list()
for (nm in names(TAXA)) {
  kd   <- if (nm %in% names(KINGDOM)) KINGDOM[[nm]] else "Animalia"
  keys <- resolve(TAXA[[nm]], kd)

  tot  <- facet_counts(keys); Sys.sleep(1)
  rep_ <- facet_counts(keys, coordinateUncertaintyInMeters =
                         paste0("0,", format(RMAX, scientific = FALSE))); Sys.sleep(1)
  hi   <- facet_counts(keys, coordinateUncertaintyInMeters = paste0("0,", RES * 3))
  Sys.sleep(1)
  if (!nrow(tot)) next

  d <- data.frame(ds = tot$ds, n = tot$n)
  d$n_rep <- rep_$n[match(d$ds, rep_$ds)];  d$n_rep[is.na(d$n_rep)] <- 0
  d$n_hi  <- hi$n[match(d$ds, hi$ds)];      d$n_hi[is.na(d$n_hi)]   <- 0
  d$k <- pmax(pmin(d$n_rep - d$n_hi, d$n), 0)          # usable count
  d <- d[d$n >= MIN_N, ]
  if (nrow(d) < 10) next
  tabs[[nm]] <- d

  N <- sum(d$n); p <- sum(d$k) / N
  r2_naive <- sum(d$n * (d$k / d$n - p)^2) / max(N * p * (1 - p), 1e-12)

  icc <- NA; sig <- NA
  fit <- try(glmer(cbind(k, n - k) ~ 1 + (1 | ds), data = d, family = binomial,
                   control = glmerControl(optimizer = "bobyqa")), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    sig <- as.numeric(VarCorr(fit)$ds[1])
    icc <- sig / (sig + pi^2 / 3)
  }

  out[[nm]] <- list(k = nrow(d), N = N, p = p, r2 = r2_naive, icc = icc, sig = sig,
                    capped = nrow(tot) >= FACET_LIM)
  cat(sprintf("%-20s %4d ds | usable %5.1f%% | naive R2 %5.1f%% | ICC %5.1f%%%s\n",
              nm, nrow(d), 100 * p, 100 * r2_naive, 100 * icc,
              if (nrow(tot) >= FACET_LIM) "  [CAPPED]" else ""))
}

cat("\n=== NAIVE R2 vs INTRACLASS CORRELATION ===\n")
cat(sprintf("%-20s %6s %12s %10s %10s %10s\n",
            "taxon", "n_ds", "records", "naive R2", "ICC", "sigma2_u"))
for (nm in names(out)) {
  o <- out[[nm]]
  cat(sprintf("%-20s %6d %12s %9.1f%% %9.1f%% %10.2f%s\n",
              nm, o$k, format(o$N, big.mark = ","), 100 * o$r2, 100 * o$icc,
              o$sig, if (o$capped) "  *" else ""))
}
cat("* still at the facet cap of ", FACET_LIM, "\n", sep = "")

ic <- vapply(out, function(o) o$icc, numeric(1))
r2 <- vapply(out, function(o) o$r2,  numeric(1))
cat(sprintf("\nICC       : %.1f%% to %.1f%% (median %.1f%%)\n",
            100 * min(ic, na.rm = TRUE), 100 * max(ic, na.rm = TRUE),
            100 * median(ic, na.rm = TRUE)))
cat(sprintf("naive R2  : %.1f%% to %.1f%% (median %.1f%%)\n",
            100 * min(r2), 100 * max(r2), 100 * median(r2)))

cat("\n--- WHICH NUMBER GOES IN THE ABSTRACT ---\n")
cat("The ICC is the defensible one: it is the share of latent-scale variance\n")
cat("attributable to the publisher, penalised for the number of groups.\n")
cat("Report the ICC. Keep the naive R2 out of the manuscript entirely, or put\n")
cat("it in the supplement as the crude version with its caveat attached.\n")
cat("If the ICC is much lower than the naive R2 but still large in absolute\n")
cat("terms, the claim survives with a smaller headline number, which is fine.\n")

saveRDS(list(summary = out, tables = tabs), "results/c8_icc.rds")
cat(sprintf("\nwrote results/c8_icc.rds\naccessed: %s\n", Sys.Date()))
