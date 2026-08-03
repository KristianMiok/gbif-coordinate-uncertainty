# =============================================================================
# Grid-centroid detection and COS applicability, GBIF crayfish
#
# Two questions, both answerable from metadata alone:
#
#  A. Which records are GRID CENTROIDS rather than observations?
#     Signature: reported radius = s / sqrt(2), the half-diagonal of a square
#     cell of side s (GridDER's definition). For these the support is a SQUARE
#     and the coordinate is a deterministic centroid, not a displaced point --
#     so BOTH the shape and the reporting model assumed by COS are wrong.
#
#  B. For what fraction of the data can COS do anything at all?
#     Hefley (2017) notes the correction cannot work when the support sits
#     inside one raster cell, because the intensity is constant there.
#
# Reads the download already on disk. No network access.
# =============================================================================

library(rgbif)
library(dplyr)

DL <- "0022943-260721160103020"     # DOI 10.15468/dl.99ezk2

d <- occ_download_import(occ_download_get(DL))
r <- d$coordinateUncertaintyInMeters
cat(sprintf("imported %s records\n\n", format(nrow(d), big.mark = ",")))

# --- 0. sanity: strip impossible radii ---------------------------------------
EARTH_QUARTER <- 10018754          # quarter of Earth's circumference, m
bad <- !is.na(r) & (r <= 0 | r > EARTH_QUARTER)
cat(sprintf("--- IMPLAUSIBLE RADII ---\nr <= 0 or r > %s m : %s records (%.3f%%)\n",
            format(EARTH_QUARTER, big.mark = ","),
            format(sum(bad), big.mark = ","), 100 * mean(bad)))
cat(sprintf("largest retained: %s m\n\n",
            format(round(max(r[!is.na(r) & !bad])), big.mark = ",")))
r[bad] <- NA

# --- A1. sqrt(2) signature ---------------------------------------------------
# a value is a half-diagonal if r * sqrt(2) lands on a "round" cell side
roundish <- function(v, tol = 0.01) {
  cand <- c(1, 2, 5) * rep(10^(0:6), each = 3)          # 1,2,5,10,20,50,...
  m <- outer(v, cand, function(a, b) abs(a - b) / b)
  hit <- apply(m, 1, min) < tol
  side <- cand[apply(m, 1, which.min)]
  list(hit = hit, side = side)
}

tab <- as.data.frame(table(r[!is.na(r)]), stringsAsFactors = FALSE)
names(tab) <- c("radius_m", "n")
tab$radius_m <- as.numeric(tab$radius_m)
tab <- tab[order(-tab$n), ]

hd <- roundish(tab$radius_m * sqrt(2))
tab$implied_cell_m <- ifelse(hd$hit, hd$side, NA)

grid_vals <- tab$radius_m[!is.na(tab$implied_cell_m)]
is_grid <- !is.na(r) & r %in% grid_vals

cat("--- A. GRID-CENTROID SIGNATURE (r = cell_side / sqrt(2)) ---\n")
print(head(tab[!is.na(tab$implied_cell_m), c("radius_m", "n", "implied_cell_m")], 15),
      row.names = FALSE)
cat(sprintf("\ngrid-centroid records : %s  (%.2f%% of reporting, %.2f%% of all)\n",
            format(sum(is_grid), big.mark = ","),
            100 * sum(is_grid) / sum(!is.na(r)),
            100 * sum(is_grid) / nrow(d)))

# --- A2. independent check: centroids repeat, observations do not ------------
coord_key <- paste(d$decimalLatitude, d$decimalLongitude)
dup_rate <- function(sel) {
  k <- coord_key[sel]
  1 - length(unique(k)) / length(k)
}
other <- !is.na(r) & !is_grid
cat("\n--- A2. COORDINATE DUPLICATION (independent of the radius field) ---\n")
cat(sprintf("grid-centroid subset : %.1f%% of records share coordinates\n",
            100 * dup_rate(is_grid)))
cat(sprintf("other reporting      : %.1f%%\n", 100 * dup_rate(other)))
cat(sprintf("radius NA            : %.1f%%\n", 100 * dup_rate(is.na(r))))
cat("(centroids should repeat far more; if they do not, the signature is\n",
    " arithmetic coincidence and must not be used)\n")

# --- A3. where do they come from --------------------------------------------
cat("\n--- A3. TOP PUBLISHERS OF GRID-CENTROID RECORDS ---\n")
gc_src <- d[is_grid, ] %>%
  count(datasetKey, countryCode, sort = TRUE) %>%
  head(10)
print(as.data.frame(gc_src), row.names = FALSE)

# --- A4. geometric consequence ----------------------------------------------
cat("\n--- A4. DISC vs SQUARE ---\n")
cat("A disc of radius r circumscribes the square of side r*sqrt(2):\n")
cat(sprintf("  square area = 2r^2, disc area = pi*r^2  ->  disc is %.2fx larger\n",
            pi / 2))
cat("  the disc CONTAINS the truth, so COS stays valid, but the support is\n")
cat("  inflated by 57%% and the correction is diluted accordingly.\n")
cat("  The reporting model is the real breakage: the coordinate is a fixed\n")
cat("  centroid, not a uniform draw, so records in one cell are not\n")
cat("  independent displacements of distinct true locations.\n")

# --- B. can COS do anything at all? -----------------------------------------
cat("\n--- B. COS APPLICABILITY BY RASTER RESOLUTION ---\n")
classify <- function(cells) {
  cut(cells, breaks = c(-Inf, 0.5, 3, Inf),
      labels = c("inert (< 0.5 cell)", "marginal (0.5-3)", "COS matters (>3)"))
}
for (res in c(90, 250, 1000)) {
  cells <- r / res
  cl <- classify(cells[!is.na(cells)])
  tb <- table(cl)
  cat(sprintf("\nresolution %5d m   (median radius %.2f cells)\n",
              res, median(cells, na.rm = TRUE)))
  for (i in seq_along(tb)) {
    cat(sprintf("  %-20s %8s  (%.1f%% of reporting, %.1f%% of ALL records)\n",
                names(tb)[i], format(tb[i], big.mark = ","),
                100 * tb[i] / sum(tb), 100 * tb[i] / nrow(d)))
  }
}

# --- C. the honest partition of the dataset ----------------------------------
cat("\n--- C. PARTITION OF THE FULL DATASET (at 1 km) ---\n")
cells1k <- r / 1000
state <- ifelse(is.na(r), "no radius reported",
         ifelse(is_grid, "grid centroid (support misspecified)",
         ifelse(cells1k < 0.5, "radius inert at this resolution",
                              "usable for COS")))
p <- sort(table(state), decreasing = TRUE)
for (i in seq_along(p)) {
  cat(sprintf("  %-38s %8s  (%.1f%%)\n", names(p)[i],
              format(p[i], big.mark = ","), 100 * p[i] / nrow(d)))
}

# --- plot --------------------------------------------------------------------
png("gbif_radius_grid.png", width = 1000, height = 500)
rr <- r[!is.na(r)]
hist(log10(rr), breaks = 100, col = "grey80", border = NA,
     main = "GBIF crayfish coordinate uncertainty; red = grid-centroid signature",
     xlab = "log10(radius, m)")
abline(v = log10(grid_vals), col = "firebrick", lwd = 1)
abline(v = log10(c(90, 1000)), col = "steelblue", lty = 2, lwd = 2)
legend("topright", c("grid-centroid values", "raster resolutions"),
       col = c("firebrick", "steelblue"), lty = c(1, 2), bty = "n")
dev.off()
cat("\nwrote gbif_radius_grid.png\n")
