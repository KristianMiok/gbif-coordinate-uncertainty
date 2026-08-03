# =============================================================================
# Figure 1 -- publisher regimes, with the three fixes
#
# The exploratory version in 08_ had three defects:
#
#   1. The zero spike merged two OPPOSITE publisher behaviours. A dataset that
#      reports no radius at all has n_rep = 0, so usable = 0, and lands in the
#      same bar as a dataset reporting precise 15 m coordinates. One says "I do
#      not know", the other says "I know exactly". Separated here.
#   2. Bars counted datasets while the mean line was record-weighted, so the two
#      things on the same panel were not measured on the same footing. Both
#      versions are drawn.
#   3. Free y-axes made cross-panel comparison misleading. Panels now show
#      proportions on a common axis.
#
# Panels are ordered by the uncertainty coefficient, so the gradient in how
# strongly the publisher determines the outcome runs left to right.
#
# Reads results/c8_icc.rds and results/dataset_effect.rds. No network access.
# =============================================================================

z    <- readRDS("results/c8_icc.rds")
rows <- readRDS("results/dataset_effect.rds")
tabs <- z$tables

ord <- names(sort(vapply(rows, function(r) r$u, numeric(1)), decreasing = TRUE))
ord <- ord[ord %in% names(tabs)]

BR <- seq(0, 100, 5)
COL_NONE <- "#B03A2E"    # dataset reports no radius at all
COL_SOME <- "#5D6D7E"    # dataset reports radii, this is the usable share

classify <- function(d) {
  # a dataset is a NON-REPORTER if essentially none of its records carry a radius
  rep_rate <- d$n_rep / d$n
  list(none = rep_rate < 0.01, share = 100 * d$k / d$n)
}

panel <- function(nm, weight_by_records, ymax) {
  d <- tabs[[nm]]; cl <- classify(d)
  w <- if (weight_by_records) d$n else rep(1, nrow(d))
  h_none <- hist(cl$share[cl$none], breaks = BR, plot = FALSE)
  h_some <- hist(cl$share[!cl$none], breaks = BR, plot = FALSE)
  cn <- vapply(seq_len(length(BR) - 1), function(i)
    sum(w[cl$none  & cl$share >= BR[i] & cl$share < BR[i + 1] |
          cl$none  & i == length(BR) - 1 & cl$share == 100]), numeric(1))
  cs <- vapply(seq_len(length(BR) - 1), function(i)
    sum(w[!cl$none & cl$share >= BR[i] & cl$share < BR[i + 1] |
          !cl$none & i == length(BR) - 1 & cl$share == 100]), numeric(1))
  tot <- sum(cn) + sum(cs)
  m <- rbind(cn, cs) / tot
  bp <- barplot(m, beside = FALSE, col = c(COL_NONE, COL_SOME), border = NA,
                space = 0, ylim = c(0, ymax), axes = FALSE, xaxt = "n")
  axis(1, at = bp[seq(1, length(bp), 4)], labels = BR[seq(1, length(BR) - 1, 4)],
       cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
  p <- 100 * sum(d$k) / sum(d$n)
  abline(v = bp[max(1, ceiling(p / 5))], col = "black", lwd = 2.2, lty = 2)
  text(bp[max(1, ceiling(p / 5))], ymax * 0.96, sprintf(" mean %.0f%%", p),
       adj = 0, cex = 0.75, col = "black")
  title(sprintf("%s  (U = %.2f)", nm, rows[[nm]]$u), cex.main = 1.0, font.main = 1)
}

draw <- function(file, weight_by_records, ymax, lab) {
  png(file, width = 1200, height = 780, res = 110)
  op <- par(mfrow = c(3, 4), mar = c(3.0, 3.4, 2.2, 0.6), mgp = c(2, 0.55, 0),
            oma = c(4.6, 3.0, 3.4, 0.4))
  for (nm in ord) panel(nm, weight_by_records, ymax)
  mtext("% of a dataset's records that are usable for a support correction",
        side = 1, outer = TRUE, line = 1.2, cex = 0.85)
  mtext(lab, side = 2, outer = TRUE, line = 1.0, cex = 0.85)
  mtext("Publishers occupy regimes; the aggregate describes no dataset", side = 3, outer = TRUE,
        line = 1.6, cex = 1.05, font = 2, adj = 0)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", horiz = TRUE, inset = c(0, -0.02), xpd = NA, bty = "n",
         cex = 0.9, fill = c(COL_NONE, COL_SOME), border = NA,
         legend = c("dataset reports no radius at all", "dataset reports radii"))
  par(op); dev.off()
  cat("wrote", file, "\n")
}

draw("results/fig1_regimes_by_dataset.png", FALSE, 0.85,
     "proportion of datasets")
draw("results/fig1_regimes_by_record.png",  TRUE,  0.85,
     "proportion of records")

# --- the number the figure is making ----------------------------------------
cat("\n=== THE ZERO SPIKE, SPLIT ===\n")
cat(sprintf("%-20s %12s %12s %12s\n", "taxon",
            "ds at ~0%", "of which:", "non-reporters"))
for (nm in ord) {
  d <- tabs[[nm]]; cl <- classify(d)
  z0 <- cl$share < 1
  cat(sprintf("%-20s %11.1f%% %12s %11.1f%%\n", nm,
              100 * mean(z0), "",
              100 * sum(z0 & cl$none) / max(sum(z0), 1)))
}
cat("\nThe right-hand column is the share of the zero spike that is 'no radius\n")
cat("reported' rather than 'radius reported and small'. Where it is high, the\n")
cat("spike is publishers declining to state uncertainty; where it is low, the\n")
cat("spike is genuinely precise data for which no correction is needed. These\n")
cat("are opposite situations and the exploratory figure conflated them.\n")

cat("\n=== HOW FAR IS THE MEAN FROM ANY ACTUAL DATASET? ===\n")
cat(sprintf("%-20s %10s %14s\n", "taxon", "mean", "% ds within +-10pts"))
for (nm in ord) {
  d <- tabs[[nm]]
  p <- 100 * sum(d$k) / sum(d$n); sh <- 100 * d$k / d$n
  cat(sprintf("%-20s %9.1f%% %13.1f%%\n", nm, p, 100 * mean(abs(sh - p) <= 10)))
}
cat("\nThis is the argument against reporting an aggregate at all.\n")
