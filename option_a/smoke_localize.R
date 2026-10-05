## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
source(file.path(.ROOT, "option_a/R/localize_groups.R"))
b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
gen <- function(n, fam, C) {
  X <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
  eta0 <- as.numeric(cbind(1, X) %*% b0)
  d <- switch(fam, null = 0, intercept = C, slope = (C - 1) * eta0, link = C * (eta0^2 - 1),
              ushape = C * (X[, 1]^2 - 1), inter = C * X[, 1] * X[, 3])
  list(X = X, p = plogis(eta0), y = rbinom(n, 1, plogis(eta0 + d)))
}
set.seed(1)
cat("== external validation, n = 1000, M = 499\n")
for (fam in c("null", "intercept", "slope", "link", "ushape", "inter")) {
  C <- c(null = 0, intercept = .4, slope = .6, link = .25, ushape = .8, inter = 1)[[fam]]
  t0 <- Sys.time(); g <- gen(1000, fam, C); r <- localize_external(g$y, g$p, g$X, M = 499)
  cat(sprintf("%-9s named: %-22s  %.1fs | %s\n", fam, paste(r$named, collapse = ","), as.numeric(difftime(Sys.time(), t0, units = "secs")),
              paste(sprintf("%s=%.3f", names(r$intersection)[1:4], r$intersection[1:4]), collapse = " ")))
}
cat("== in-sample checking, n = 1000, B = 99\n")
for (fam in c("null", "link", "ushape", "inter")) {
  C <- c(null = 0, link = .25, ushape = .8, inter = 1)[[fam]]
  g <- gen(1000, fam, C); d <- data.frame(g$X, y = g$y)
  fit <- glm(y ~ x1 + x2 + x3 + x4 + x5, data = d, family = binomial())
  t0 <- Sys.time(); r <- localize_insample(fit, g$X, B = 99)
  cat(sprintf("%-9s named: %-12s %.1fs B=%d | %s\n", fam, paste(r$named, collapse = ","), as.numeric(difftime(Sys.time(), t0, units = "secs")),
              r$B_used, paste(sprintf("%s=%.3f", names(r$intersection), r$intersection), collapse = " ")))
}
