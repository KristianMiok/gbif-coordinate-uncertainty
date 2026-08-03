# =============================================================================
# Is the partition set by the publisher rather than the taxon?
#
# The multi-taxon run showed the usable share varying 11-fold across thirteen
# groups without tracking realm or record type, and the cumulative curves
# pointed at publisher defaults as the driver -- salmonids jumping from 60.7%
# below 500 m to 94.0% below 1000 m, hard corals at 69.8% below 100 m but 4.9%
# below 15 m. That is an inference from shapes, not a measurement.
#
# This measures it. If the partition varies as much BETWEEN DATASETS WITHIN one
# taxon as it does between taxa, the driver is georeferencing convention and a
# user genuinely cannot predict their own partition. If datasets within a taxon
# are homogeneous, the taxon is the unit after all and the claim must be weakened.
#
# Second question in the same pass: are the internally inconsistent records real
# false precision, or publishers deliberately generalising coordinates for
# protected species? Concentration in a few datasets means the latter.
#
# Reads the download already on disk. No network access.
# =============================================================================

library(rgbif)
library(dplyr)

DL   <- "0022943-260721160103020"   # DOI 10.15468/dl.99ezk2
FAKE <- c(301, 3036, 999, 9999)     # GBIF-documented geocoder defaults
RES  <- 1000                        # partition resolution, m
MIN_N <- 500                        # datasets smaller than this are unstable

d <- occ_download_import(occ_download_get(DL))
cat(sprintf("imported %s records\n", format(nrow(d), big.mark = ",")))

r <- d$coordinateUncertaintyInMeters
r[!is.na(r) & (r <= 0 | r > 1e7)] <- NA

# coordinate-implied lower bound (few-decimals direction only; see 03_)
sig_dec <- function(v) {
  s <- format(v, scientific = FALSE, trim = TRUE)
  s <- sub("0+$", "", s); s <- sub("\\.$", "", s)
  frac <- sub("^[^.]*$", "", s); frac <- sub("^[^.]*\\.", "", frac)
  nchar(frac)
}
M <- 111320
r_min <- sqrt((10^(-sig_dec(d$decimalLatitude)) * M)^2 +
              (10^(-sig_dec(d$decimalLongitude)) * M *
                 cos(d$decimalLatitude * pi / 180))^2) / 2

state <- ifelse(is.na(r), "unreported",
         ifelse(r %in% FAKE, "geocoder_default",
         ifelse(r < RES / 2, "inert",
         ifelse(r > RES * 3, "usable", "marginal"))))
inconsistent <- !is.na(r) & r < r_min

# --- per dataset -------------------------------------------------------------
by_ds <- data.frame(ds = d$datasetKey, state = state,
                    inc = inconsistent, r = r) %>%
  group_by(ds) %>%
  summarise(n = n(),
            usable      = mean(state == "usable"),
            inert       = mean(state == "inert"),
            unreported  = mean(state == "unreported"),
            inconsistent = mean(inc),
            med_r       = suppressWarnings(median(r, na.rm = TRUE)),
            .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  arrange(desc(n))

cat(sprintf("\n%s datasets with >= %s records (%.1f%% of all records)\n",
            nrow(by_ds), MIN_N,
            100 * sum(by_ds$n) / nrow(d)))

cat("\n=== THE 15 LARGEST DATASETS ===\n")
cat(sprintf("%-38s %8s %8s %8s %10s %8s\n",
            "datasetKey", "n", "usable", "inert", "unreport", "med r"))
for (i in seq_len(min(15, nrow(by_ds)))) {
  b <- by_ds[i, ]
  cat(sprintf("%-38s %8s %7.1f%% %7.1f%% %9.1f%% %8s\n",
              b$ds, format(b$n, big.mark = ","),
              100 * b$usable, 100 * b$inert, 100 * b$unreported,
              ifelse(is.na(b$med_r), "-", format(round(b$med_r)))))
}

# --- the comparison that decides the claim ----------------------------------
cat("\n=== SPREAD: BETWEEN DATASETS (one taxon) vs BETWEEN TAXA ===\n")
u <- by_ds$usable
qs <- quantile(u, c(0, .1, .25, .5, .75, .9, 1))
cat("usable share across crayfish datasets:\n")
print(round(100 * qs, 1))
cat(sprintf("\n  full range      : %.1f%% to %.1f%%  (%.0f-fold)\n",
            100 * min(u), 100 * max(u), max(u) / max(min(u), 0.001)))
cat(sprintf("  10th-90th pct   : %.1f%% to %.1f%%  (%.1f-fold)\n",
            100 * qs[2], 100 * qs[6], qs[6] / max(qs[2], 0.001)))
cat("\n  for reference, ACROSS THIRTEEN TAXA at 1 km: 4.5%% to 47.4%% (11-fold)\n")
cat("\nIf the within-taxon spread is comparable, the taxon is not the unit and\n")
cat("a user cannot predict their partition from the group they work on.\n")

# --- how concentrated is the variation --------------------------------------
cat("\n=== HOW MUCH OF THE VARIATION IS BETWEEN DATASETS? ===\n")
rec <- data.frame(ds = d$datasetKey, usable = state == "usable") %>%
  filter(ds %in% by_ds$ds)
grand <- mean(rec$usable)
ssb <- sum(by_ds$n * (by_ds$usable - grand)^2)
sst <- sum((rec$usable - grand)^2)
cat(sprintf("overall usable share      : %.1f%%\n", 100 * grand))
cat(sprintf("variance explained by dataset identity : %.1f%%\n", 100 * ssb / sst))
cat("(this is a one-way ANOVA R^2 on a binary outcome -- a crude but honest\n")
cat(" summary of how much of the partition is a publisher property)\n")

# --- internally inconsistent records: real, or generalisation policy? --------
cat("\n=== ARE THE INCONSISTENT RECORDS CONCENTRATED? ===\n")
inc_tot <- sum(inconsistent)
cat(sprintf("total internally inconsistent: %s (%.2f%% of all records)\n",
            format(inc_tot, big.mark = ","), 100 * inc_tot / nrow(d)))
inc_ds <- data.frame(ds = d$datasetKey, inc = inconsistent) %>%
  group_by(ds) %>% summarise(n = n(), k = sum(inc), .groups = "drop") %>%
  filter(k > 0) %>% arrange(desc(k))
cat(sprintf("spread over %s datasets\n", nrow(inc_ds)))
top <- head(inc_ds, 8)
cat(sprintf("\n%-38s %9s %9s %8s\n", "datasetKey", "records", "inconsist", "rate"))
for (i in seq_len(nrow(top))) {
  b <- top[i, ]
  cat(sprintf("%-38s %9s %9s %7.1f%%\n", b$ds,
              format(b$n, big.mark = ","), format(b$k, big.mark = ","),
              100 * b$k / b$n))
}
cat(sprintf("\ntop 3 datasets hold %.1f%% of all inconsistent records\n",
            100 * sum(head(inc_ds$k, 3)) / inc_tot))
cat("A high concentration means deliberate coordinate generalisation by a few\n")
cat("publishers, NOT false precision. In that case the flag is a publisher\n")
cat("property and must not be presented as a record-level data-quality defect.\n")

saveRDS(by_ds, "results/crayfish_by_dataset.rds")
cat("\nwrote results/crayfish_by_dataset.rds\n")
