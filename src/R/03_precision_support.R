# =============================================================================
# Coordinate precision as a support, not a flag
#
# The reported radius is one signal of positional uncertainty. The COORDINATES
# THEMSELVES are a second, independent one: a value rounded to k significant
# decimals cannot resolve position below the size of the rounding cell, and the
# rounding is irreversible.
#
# ASYMMETRY -- the only claim that is safe:
#   few decimals  -> a genuine LOWER BOUND on uncertainty (information is gone)
#   many decimals -> NO information (format transformation manufactures digits;
#                    GBIF Georeferencing Best Practices, section on false
#                    precision, warns against reading precision off the value)
#
# Everything below uses only the safe direction.
#
# Three uses, none of which is a flag:
#   1. records with r = NA get a support from their coordinates
#   2. records claiming r below their own coordinate resolution are internally
#      inconsistent -- false precision, detectable without ground truth
#   3. documented geocoder defaults (301, 3036, 999, 9999) are not measurements
#
# Reads the download already on disk. No network access.
# =============================================================================

library(rgbif)
library(dplyr)

DL <- "0022943-260721160103020"     # DOI 10.15468/dl.99ezk2
FAKE <- c(301, 3036, 999, 9999)     # GBIF-documented geocoder defaults

d <- occ_download_import(occ_download_get(DL))
cat(sprintf("imported %s records\n\n", format(nrow(d), big.mark = ",")))

r   <- d$coordinateUncertaintyInMeters
lat <- d$decimalLatitude
lon <- d$decimalLongitude

# --- significant decimals, trailing zeros stripped ---------------------------
# 45.50000 carries the same information as 45.5; counting stored digits would
# manufacture precision, which is exactly the failure mode being avoided.
sig_dec <- function(v) {
  s <- format(v, scientific = FALSE, trim = TRUE)
  s <- sub("0+$", "", s)                 # strip trailing zeros
  s <- sub("\\.$", "", s)                # strip a bare trailing point
  frac <- sub("^[^.]*$", "", s)          # "" when there is no decimal point
  frac <- sub("^[^.]*\\.", "", frac)
  nchar(frac)
}

dlat <- sig_dec(lat)
dlon <- sig_dec(lon)

# --- lower bound on uncertainty implied by the rounding cell -----------------
# rounding to k decimals puts the truth inside a cell of side 10^-k degrees;
# the smallest circle containing that cell has radius = half its diagonal,
# which is the Darwin Core definition of coordinateUncertaintyInMeters.
M_PER_DEG <- 111320
cell_lat_m <- 10^(-dlat) * M_PER_DEG
cell_lon_m <- 10^(-dlon) * M_PER_DEG * cos(lat * pi / 180)
r_min <- sqrt(cell_lat_m^2 + cell_lon_m^2) / 2

cat("--- COORDINATE RESOLUTION ---\n")
cat("significant decimals (lat):\n")
print(table(pmin(dlat, 8)))
cat(sprintf("\nimplied lower bound r_min (m): median %.1f | 25%% %.1f | 75%% %.1f\n",
            median(r_min), quantile(r_min, .25), quantile(r_min, .75)))

# --- 1. what the coordinates give the NA records ----------------------------
na_r <- is.na(r)
cat(sprintf("\n--- 1. RECORDS WITH NO REPORTED RADIUS (n = %s, %.1f%%) ---\n",
            format(sum(na_r), big.mark = ","), 100 * mean(na_r)))
cat("these are currently either discarded or silently kept; the coordinates\n")
cat("give every one of them a defensible conservative support:\n\n")
qn <- quantile(r_min[na_r], c(.05, .25, .5, .75, .95))
print(round(qn, 1))
for (res in c(90, 1000)) {
  cells <- r_min[na_r] / res
  cat(sprintf("  at %4d m resolution: %.1f%% inert (<0.5 cell), %.1f%% usable (>3 cells)\n",
              res, 100 * mean(cells < 0.5), 100 * mean(cells > 3)))
}

# --- 2. internal inconsistency (false precision) ----------------------------
rep_r <- !is.na(r)
inconsistent <- rep_r & (r < r_min)
cat(sprintf("\n--- 2. INTERNALLY INCONSISTENT RECORDS ---\n"))
cat(sprintf("reported radius smaller than the coordinates can support:\n"))
cat(sprintf("  %s records  (%.2f%% of reporting, %.2f%% of all)\n",
            format(sum(inconsistent), big.mark = ","),
            100 * sum(inconsistent) / sum(rep_r),
            100 * sum(inconsistent) / nrow(d)))
if (sum(inconsistent) > 0) {
  ratio <- r_min[inconsistent] / r[inconsistent]
  cat(sprintf("  understatement factor r_min/r: median %.1fx | 90%% %.1fx | max %.0fx\n",
              median(ratio), quantile(ratio, .9), max(ratio)))
  cat("\n  most common claimed radii among inconsistent records:\n")
  print(head(sort(table(r[inconsistent]), decreasing = TRUE), 8))
}
cat("\n  NOTE the direction: these records UNDERSTATE their uncertainty.\n")
cat("  That is delta_size < 1 -- the support is too small, the truth can fall\n")
cat("  outside it, and a change-of-support correction cannot recover it.\n")
cat("  A threshold filter keeps exactly these records.\n")

# --- 3. documented geocoder defaults ----------------------------------------
is_fake <- rep_r & r %in% FAKE
cat(sprintf("\n--- 3. DOCUMENTED GEOCODER DEFAULTS (301, 3036, 999, 9999) ---\n"))
cat(sprintf("  %s records  (%.2f%% of reporting)\n",
            format(sum(is_fake), big.mark = ","),
            100 * sum(is_fake) / sum(rep_r)))
if (sum(is_fake) > 0) print(table(r[is_fake]))
cat(sprintf("  percentile of 3036 in the reported-radius distribution: %.1f%%\n",
            100 * mean(r[rep_r] <= 3036)))
cat("  these are software artefacts, not measurements; 301 in particular\n")
cat("  usually marks a country centroid, where true uncertainty is far larger.\n")

# --- 4. what the standard filter does to each state -------------------------
cat("\n--- 4. WHAT `r < 1000` DOES TO EACH STATE ---\n")
state <- ifelse(na_r, "no radius reported",
         ifelse(is_fake, "geocoder default (not a measurement)",
         ifelse(inconsistent, "false precision (understated)",
                              "internally consistent")))
kept <- !na_r & r < 1000
tb <- table(state, ifelse(kept, "kept", "dropped"))
print(tb)
# compute the retention rate per state instead of asserting a direction
keep_rate <- tapply(kept[!na_r], state[!na_r], mean)
cat("\nretention rate under the filter, by state:\n")
print(round(100 * sort(keep_rate, decreasing = TRUE), 1))
cat("\nRead this off the numbers. Do not assume a direction: in this dataset\n")
cat("false-precision records are dominated by claims of r = 5000, which the\n")
cat("threshold removes, so they are retained LESS often than consistent ones.\n")
cat("The defensible claim is narrower -- the filter cannot distinguish a\n")
cat("software default from a measurement, and it discards honest reporters.\n")

# --- 5. what the lower bound can and cannot do ------------------------------
cat("\n--- 5. THE LOWER BOUND IS NOT A SUPPORT FOR NA RECORDS ---\n")
cat("r_min bounds uncertainty from BELOW. For a record with no reported\n")
cat("radius it says only 'at least this much', which is uninformative:\n")
cat("a locality-derived coordinate can be stored with six decimals and still\n")
cat("be kilometres wrong. Treating r_min as the support for NA records would\n")
cat("manufacture false precision -- delta_size < 1, the dangerous direction.\n")
cat(sprintf("\n  median r_min among NA records: %.1f m  (i.e. no constraint)\n",
            median(r_min[na_r])))

cat("\nWhere the bound DOES bind is on records that report a radius smaller\n")
cat("than their own coordinates can carry. There the correction is upward,\n")
cat("and enlarging a support is always safe.\n")
r_use <- r
r_use[rep_r] <- pmax(r[rep_r], r_min[rep_r])
r_use[is_fake] <- NA          # a software default is not a measurement
cat(sprintf("\n  supports enlarged : %s (%.1f%% of reporting)\n",
            format(sum(rep_r & r_use > r, na.rm = TRUE), big.mark = ","),
            100 * sum(rep_r & r_use > r, na.rm = TRUE) / sum(rep_r)))
cat(sprintf("  records still without any support: %s (%.1f%%)\n",
            format(sum(is.na(r_use)), big.mark = ","),
            100 * mean(is.na(r_use))))

for (res in c(90, 1000)) {
  cells <- r_use[!is.na(r_use)] / res
  cat(sprintf("\namong records WITH a support, at %4d m: inert %.1f%% | marginal %.1f%% | COS matters %.1f%%\n",
              res, 100 * mean(cells < 0.5),
              100 * mean(cells >= 0.5 & cells <= 3), 100 * mean(cells > 3)))
}

cat("\n--- CAVEATS THAT MUST TRAVEL WITH THIS ---\n")
cat("* only the FEW-decimals direction is used; many decimals are treated as\n")
cat("  uninformative, per GBIF Georeferencing Best Practices.\n")
cat("* publishers who generalise coordinates for conservation or privacy\n")
cat("  reasons will appear here as false precision when they are not.\n")
cat("  Check the inconsistent subset by dataset before interpreting.\n")
cat("* r_min assumes rounding, not truncation, and a spherical Earth.\n")
