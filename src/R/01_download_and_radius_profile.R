# =============================================================================
# GBIF coordinate-uncertainty profile for freshwater crayfish
#
# Produces the three numbers that decide the COS design:
#   1. what fraction of modelable records report a radius at all
#   2. the shape of the radius distribution (Marcer et al. 2022 found spikes,
#      not a smooth distribution -- spikes are what make FFT classes cheap)
#   3. radius expressed in RASTER CELLS, which is what the integrator costs
#
# Run PART 1, wait for the email / poll, then run PART 2.
# =============================================================================

library(rgbif)
library(dplyr)

# Credentials must be in ~/.Renviron as GBIF_USER, GBIF_PWD, GBIF_EMAIL
# usethis::edit_r_environ()   # then restart R

# ---------------------------------------------------------------- PART 1 ----
# Resolve taxon keys. Astacoidea (N hemisphere) + Parastacoidea (S hemisphere).
fams <- c("Astacidae", "Cambaridae", "Cambaroididae", "Parastacidae")
bb <- name_backbone_checklist(fams)
print(bb[, c("verbatim_name", "usageKey", "rank", "status", "matchType")])
keys <- bb$usageKey[!is.na(bb$usageKey)]
stopifnot(length(keys) == length(fams))

# Denominators BEFORE any filtering -- needed for the "NA as third state" claim
tot_all <- sum(sapply(keys, function(k) occ_count(taxonKey = k)))
tot_geo <- sum(sapply(keys, function(k) occ_count(taxonKey = k, hasCoordinate = TRUE)))
cat(sprintf("\nTotal crayfish records        : %s\n", format(tot_all, big.mark = ",")))
cat(sprintf("With coordinates              : %s (%.1f%%)\n",
            format(tot_geo, big.mark = ","), 100 * tot_geo / tot_all))

dl <- occ_download(
  pred_in("taxonKey", keys),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV"
)
print(dl)
cat("\nDownload key:", dl[1], "\n")
cat("Poll with:  occ_download_wait('", dl[1], "')\n", sep = "")

# occ_download_wait(dl)          # blocks until ready (minutes to hours)

# ---------------------------------------------------------------- PART 2 ----
# dl <- "0012345-260101000000000"   # paste the key here if resuming a new session

d <- occ_download_get(dl, overwrite = TRUE) %>% occ_download_import()
cat(sprintf("\nImported %s records\n", format(nrow(d), big.mark = ",")))

r <- d$coordinateUncertaintyInMeters

# --- 1. reporting rate -------------------------------------------------------
n_tot <- length(r)
n_rep <- sum(!is.na(r))
cat(sprintf("\n--- REPORTING ---\n"))
cat(sprintf("radius reported     : %s / %s  (%.1f%%)\n",
            format(n_rep, big.mark = ","), format(n_tot, big.mark = ","),
            100 * n_rep / n_tot))
cat(sprintf("radius NA           : %s  (%.1f%%)  <- the third state\n",
            format(n_tot - n_rep, big.mark = ","), 100 * (n_tot - n_rep) / n_tot))

# how the standard r < 1000 filter treats each state
kept_reported <- sum(!is.na(r) & r < 1000)
cat(sprintf("\nunder a naive `r < 1000` filter:\n"))
cat(sprintf("  honest reporters kept : %s (%.1f%% of reporters)\n",
            format(kept_reported, big.mark = ","), 100 * kept_reported / n_rep))
cat(sprintf("  honest reporters DROPPED: %s\n",
            format(n_rep - kept_reported, big.mark = ",")))
cat(sprintf("  non-reporters silently kept or dropped depending on na.rm: %s\n",
            format(n_tot - n_rep, big.mark = ",")))

# --- 2. distribution shape ---------------------------------------------------
rr <- r[!is.na(r) & r > 0]
cat(sprintf("\n--- DISTRIBUTION (m) ---\n"))
print(round(quantile(rr, c(0, .05, .25, .5, .75, .9, .95, .99, 1))))
cat(sprintf("mean %.0f  |  n distinct values %s\n",
            mean(rr), format(length(unique(rr)), big.mark = ",")))

cat("\n--- TOP 20 EXACT VALUES (spikes) ---\n")
spikes <- sort(table(rr), decreasing = TRUE)[1:20]
print(data.frame(radius_m = as.numeric(names(spikes)),
                 n = as.integer(spikes),
                 pct = round(100 * as.integer(spikes) / length(rr), 2)))
cat(sprintf("\ntop 20 values cover %.1f%% of reporting records\n",
            100 * sum(spikes) / length(rr)))

# --- 3. cost in raster cells -------------------------------------------------
cat("\n--- RADIUS IN CELLS, AND INTEGRATOR COST ---\n")
# constants fitted to the measured benchmark on your machine
c_naive <- 9.2e-9   # s per (record x cell) for the per-record loop
c_fft   <- 1.5e-8   # s per (class x raster cell) for one convolution

for (res in c(90, 1000)) {          # Hydrography90m, CHELSA ~1km
  cells <- rr / res
  cpd <- mean(pi * cells^2)         # mean cells per disc
  raster_cells <- 4e6               # adjust to your study extent
  n_class <- min(20, length(unique(round(cells))))
  n_cross <- n_class * raster_cells * c_fft / (cpd * c_naive)
  cat(sprintf("\nresolution %4d m:\n", res))
  cat(sprintf("  median radius %.1f cells | mean cells per disc %.0f\n",
              median(cells), cpd))
  cat(sprintf("  distinct radius classes (rounded) %d\n", n_class))
  cat(sprintf("  naive/FFT crossover at n = %s records\n",
              format(round(n_cross), big.mark = ",")))
  cat(sprintf("  -> at your n (%s): use %s\n", format(n_rep, big.mark = ","),
              ifelse(n_rep < n_cross, "NAIVE per-record loop", "FFT convolution")))
}

# --- plot --------------------------------------------------------------------
png("gbif_radius_hist.png", width = 900, height = 500)
hist(log10(rr), breaks = 80, col = "grey70", border = NA,
     main = "GBIF coordinateUncertaintyInMeters, crayfish",
     xlab = "log10(radius, m)")
abline(v = log10(c(100, 1000, 10000)), col = "firebrick", lty = 2)
dev.off()
cat("\nwrote gbif_radius_hist.png\n")

# keep the citable DOI for the manuscript
cat("\nGBIF download DOI (cite this):", occ_download_meta(dl)$doi, "\n")
