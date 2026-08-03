# =============================================================================
# How much does knowing the publisher tell you about whether a record is usable?
#
# Both variance models failed: the one-way ANOVA R^2 is optimistic with hundreds
# of groups, and the binomial GLMM degenerates because so many datasets sit at
# exactly 0% or 100% usable -- the separation IS the finding, so any model that
# chokes on separation is the wrong tool.
#
# The information-theoretic version does not have that problem. A dataset that
# is entirely usable, or entirely not, contributes ZERO conditional entropy,
# which is exactly the behaviour wanted here.
#
#   U = 1 - H(usable | dataset) / H(usable)          [uncertainty coefficient]
#
# U = 0 : the publisher tells you nothing
# U = 1 : the publisher tells you everything
#
# U is still inflated by having many groups, so it is compared against a
# PERMUTATION NULL: usable records reshuffled across datasets with dataset sizes
# and the overall usable count held fixed. The reported excess is what survives.
#
# Reads results/c8_icc.rds. NO network access.
# =============================================================================

set.seed(20260803)
NPERM <- 200

z <- readRDS("results/c8_icc.rds")
tabs <- z$tables
stopifnot(length(tabs) > 0)

H <- function(p) {                       # binary entropy, bits
  p <- pmin(pmax(p, 0), 1)
  out <- ifelse(p == 0 | p == 1, 0, -(p * log2(p) + (1 - p) * log2(1 - p)))
  out
}

uncertainty_coef <- function(n, k) {
  N <- sum(n); p <- sum(k) / N
  h_y <- H(p)
  if (h_y <= 0) return(NA_real_)
  h_y_given_x <- sum((n / N) * H(k / n))
  1 - h_y_given_x / h_y
}

gini <- function(x, w) {                 # weighted Gini of the per-dataset share
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  cx <- cumsum(x * w) / sum(x * w)
  1 - sum((cw - c(0, head(cw, -1))) * (cx + c(0, head(cx, -1))))
}

perm_null <- function(n, k, B = NPERM) {
  N <- sum(n); K <- sum(k)
  vapply(seq_len(B), function(i) {
    # multivariate hypergeometric: deal K usable records into datasets of size n
    rem_N <- N; rem_K <- K
    kk <- numeric(length(n))
    for (j in seq_along(n)) {
      kk[j] <- rhyper(1, m = rem_K, n = rem_N - rem_K, k = n[j])
      rem_K <- rem_K - kk[j]; rem_N <- rem_N - n[j]
    }
    uncertainty_coef(n, kk)
  }, numeric(1))
}

cat("=== HOW MUCH THE PUBLISHER DETERMINES USABILITY ===\n")
cat("U = 1 - H(usable | dataset) / H(usable), against a permutation null\n\n")
cat(sprintf("%-20s %5s %12s %7s %8s %8s %8s %7s %7s\n",
            "taxon", "n_ds", "records", "usable", "U", "U_null", "excess",
            "%ds~0", "%ds~1"))

rows <- list()
for (nm in names(tabs)) {
  d <- tabs[[nm]]
  n <- d$n; k <- d$k
  p <- sum(k) / sum(n)
  u <- uncertainty_coef(n, k)
  nul <- perm_null(n, k)
  sh <- k / n
  rows[[nm]] <- list(k = length(n), N = sum(n), p = p, u = u,
                     u0 = mean(nul), sd0 = sd(nul),
                     zero = mean(sh < 0.01), full = mean(sh > 0.99),
                     g = gini(pmax(sh, 1e-9), n))
  r <- rows[[nm]]
  cat(sprintf("%-20s %5d %12s %6.1f%% %7.3f %8.3f %8.3f %6.1f%% %6.1f%%\n",
              nm, r$k, format(r$N, big.mark = ","), 100 * r$p,
              r$u, r$u0, r$u - r$u0, 100 * r$zero, 100 * r$full))
}

us <- vapply(rows, function(r) r$u, numeric(1))
ex <- vapply(rows, function(r) r$u - r$u0, numeric(1))
u0 <- vapply(rows, function(r) r$u0, numeric(1))

cat(sprintf("\nU        : %.3f to %.3f (median %.3f)\n",
            min(us), max(us), median(us)))
cat(sprintf("null U   : %.3f to %.3f (median %.3f)  <- what many groups buy you\n",
            min(u0), max(u0), median(u0)))
cat(sprintf("excess   : %.3f to %.3f (median %.3f)\n",
            min(ex), max(ex), median(ex)))

cat("\n=== REGIMES, NOT GRADIENTS ===\n")
cat(sprintf("%-20s %12s %12s %10s\n", "taxon", "% ds ~0%", "% ds ~100%", "Gini"))
for (nm in names(rows)) {
  r <- rows[[nm]]
  cat(sprintf("%-20s %11.1f%% %11.1f%% %10.3f\n",
              nm, 100 * r$zero, 100 * r$full, r$g))
}
bim <- vapply(rows, function(r) r$zero + r$full, numeric(1))
cat(sprintf("\ndatasets in a pure regime (either extreme): %.1f%% to %.1f%% (median %.1f%%)\n",
            100 * min(bim), 100 * max(bim), 100 * median(bim)))

# --- figure ------------------------------------------------------------------
png("results/dataset_regimes.png", width = 1100, height = 700)
op <- par(mfrow = c(3, 4), mar = c(3.2, 3.2, 2.4, 0.8), mgp = c(2, 0.6, 0))
for (nm in names(tabs)) {
  d <- tabs[[nm]]
  hist(100 * d$k / d$n, breaks = seq(0, 100, 5), col = "grey65", border = "white",
       main = nm, xlab = "% usable in dataset", ylab = "datasets")
  abline(v = 100 * rows[[nm]]$p, col = "firebrick", lwd = 2)
}
par(op); dev.off()
cat("\nwrote results/dataset_regimes.png\n")

saveRDS(rows, "results/dataset_effect.rds")

cat("\n--- WHAT TO REPORT ---\n")
cat("Headline: the share of datasets sitting in a pure regime, and U with its\n")
cat("permutation null alongside. Both are assumption-free and neither breaks\n")
cat("under the separation that is the finding.\n")
cat("Do NOT report the ANOVA R2 or the GLMM ICC; both are in the log as failed.\n")
cat("If the excess over the null is small, the dataset effect is mostly an\n")
cat("artefact of having many groups and the claim must be weakened.\n")
