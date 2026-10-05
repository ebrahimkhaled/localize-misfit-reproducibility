## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# check_T1.R -- numerical check of Lemma 1 / Theorem 1 (what the in-sample fit erases).
# (a) logit fit: sum(y - p) and sum((y - p) * eta) are zero to machine precision, so the logistic recalibration of
#     the fitted risks on the same data returns intercept 0 and slope 1 exactly, whatever the truth is;
# (b) probit fit (non-canonical): the same sums are not zero and the recalibration is not the identity, so the
#     intercept and slope groups keep information in-sample;
# (c) population version: at the pseudo-true logit fit of a misspecified truth, the probability residual has no
#     component along 1, x1, ..., xd (checked on one large sample, n = 200,000).
set.seed(20261005)
recal <- function(y, p) unname(coef(glm(y ~ qlogis(p), family = binomial())))
one <- function(n, truth) {
  X <- matrix(rnorm(n * 3), n, 3)
  eta <- -0.4 + X %*% c(0.8, -0.6, 0.5)
  lp <- switch(truth, logit = eta, curved = eta + 0.35 * eta^2, inter = eta + 0.9 * X[, 1] * X[, 2])
  y <- rbinom(n, 1, plogis(lp))
  fl <- glm(y ~ X, family = binomial("logit")); fp <- glm(y ~ X, family = binomial("probit"))
  pl <- fitted(fl); pp <- fitted(fp)
  c(logit_sum = sum(y - pl), logit_sum_eta = sum((y - pl) * qlogis(pl)),
    logit_recal_a = recal(y, pl)[1], logit_recal_b = recal(y, pl)[2],
    probit_sum = sum(y - pp), probit_sum_eta = sum((y - pp) * qlogis(pp)),
    probit_recal_a = recal(y, pp)[1], probit_recal_b = recal(y, pp)[2])
}
for (truth in c("logit", "curved", "inter")) {
  R <- t(replicate(200, one(1000, truth)))
  cat("\ntruth =", truth, "(200 datasets, n = 1000): max |.| and mean over datasets\n")
  print(round(rbind(max_abs = apply(abs(R), 2, max), mean = colMeans(R)), 6))
}
## (c) population erasure at the pseudo-true fit
n <- 2e5; X <- matrix(rnorm(n * 3), n, 3); eta <- -0.4 + X %*% c(0.8, -0.6, 0.5)
pi_true <- plogis(eta + 0.35 * eta^2 + 0.9 * X[, 1] * X[, 2])
y <- rbinom(n, 1, pi_true); f <- glm(y ~ X, family = binomial())
ps <- fitted(f); r <- as.numeric(pi_true) - ps                      # probability residual at (approximately) the pseudo-true fit
cat("\n(c) n = 2e5: E[(pi - p*) x~] (should be ~ n^-1/2 noise only):\n")
print(signif(colMeans(cbind(1, X) * r), 3))
cat("versus E[|pi - p*|] =", signif(mean(abs(r)), 3), "\n")
## machine-readable summary for the appendix table
set.seed(20261005)
S <- do.call(rbind, lapply(c("logit", "curved", "inter"), function(truth) {
  R <- t(replicate(200, one(1000, truth)))
  R[, c("logit_recal_b", "probit_recal_b")] <- R[, c("logit_recal_b", "probit_recal_b")] - 1   # distance from slope 1
  data.frame(truth = truth, stat = colnames(R), max_abs = apply(abs(R), 2, max), mean = colMeans(R))
}))
write.csv(S, "CHECK_T1.csv", row.names = FALSE)
