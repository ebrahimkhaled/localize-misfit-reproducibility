## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# check_multiplier_var.R -- a check of the multiplier variance: SLOPE member (eta projected empirically off 1)
# under a probability-scale-pure intercept departure pi = p0 + 0.9 w0 (n = 2000). The variance of the observed score
# U should match the robust variance computed from the recentred residuals e = r0 - w Z theta (ratio near 1); the
# uncorrected robust variance from r0 overestimates it (ratio about .9).
suppressMessages(library(splines))
source(file.path(.ROOT, "option_a/R/localize_groups.R"))
set.seed(1); n <- 2000; b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
out <- t(replicate(1500, {
  X <- matrix(rnorm(n * 5), n, 5); eta <- drop(cbind(1, X) %*% b0); p <- plogis(eta); w <- p * (1 - p)
  y <- rbinom(n, 1, p + 0.9 * w); r0 <- y - p; one <- matrix(1, n, 1)
  Bp <- as.numeric(.wperp(matrix(eta), one, w))
  ee <- w * as.numeric(.wperp(matrix(r0 / w), one, w))
  c(U = sum(Bp * r0), V_recentred = sum(Bp^2 * ee^2), V_raw = sum(Bp^2 * r0^2))
}))
cat(sprintf("var(U) / mean V: recentred %.3f, raw %.3f (1,500 data sets; SE of a variance ratio about %.3f)\n",
            var(out[, "U"]) / mean(out[, "V_recentred"]), var(out[, "U"]) / mean(out[, "V_raw"]), sqrt(2 / 1499)))
