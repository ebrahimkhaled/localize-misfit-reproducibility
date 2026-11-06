# localize_size.R -- how large each part of the misfit is, on the risk scale, for frozen predictions p on validation data.
# The four nested updates of the calibration hierarchy are fitted on the data:
#   p_I  intercept update            logit p_I = a + eta0
#   p_S  logistic recalibration      logit p_S = a + b eta0
#   p_L  flexible recalibration      logit p_L = ns(eta0, 4)
#   p_C  revision of the covariates  logit p_L plus the main effects and the COV bases of localize_external() (squares and cubes of the
#        continuous covariates, their pairwise products and natural splines), fitted as one logistic model
# The SIZE of a part is the mean absolute change in predicted risk that its update makes beyond the simpler one
# (INTERCEPT: |p_I - p|, SLOPE: |p_S - p_I|, LINK: |p_L - p_S|, COV: |p_C - p_L|), in percentage points.
# An update fitted to noise also moves the predictions, so each size is reported with the size reached by chance: the
# mean over R outcome vectors drawn from the simpler update (p for INTERCEPT, p_I for SLOPE, ...), under which the
# part is absent. EXCESS = size - chance. Its 95% noise interval is size minus the 97.5% and 2.5% points of the
# chance draws: it removes the noise once, and under no misfit it contains 0 with probability 0.95. Also returned: the calibration intercept and slope with Wald 95% intervals,
# and the integrated calibration index (ICI; Austin and Steyerberg 2019) of p.
localize_size <- function(y, p, X, R = 200L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X); eta <- stats::qlogis(p); n <- length(y)
  C <- cbind(X, do.call(cbind, .group_bases(eta, X)$COV))
  fits <- function(yy, idx = seq_len(n)) {
    e <- eta[idx]; yy <- yy[idx]
    pI <- stats::fitted(stats::glm(yy ~ 1, offset = e, family = stats::binomial()))
    pS <- stats::fitted(stats::glm(yy ~ e, family = stats::binomial()))
    S4 <- splines::ns(e, df = 4)
    pL <- stats::fitted(stats::glm(yy ~ S4, family = stats::binomial()))
    Ck <- C[idx, , drop = FALSE]
    pC <- stats::fitted(suppressWarnings(stats::glm(yy ~ S4 + Ck, family = stats::binomial())))
    list(p = c(list(stats::plogis(e)), list(pI, pS, pL, pC)))
  }
  sizes <- function(f) vapply(1:4, function(k) 100 * mean(abs(f$p[[k + 1]] - f$p[[k]])), 0)
  obs_f <- fits(y)
  obs <- sizes(obs_f)
  ## chance level of part k: outcomes drawn from the simpler update, under which part k is absent
  draws <- vapply(1:4, function(k) {
    base <- obs_f$p[[k]]
    replicate(R, sizes(fits(stats::rbinom(n, 1L, base)))[k])
  }, numeric(R))
  chance <- colMeans(draws); q <- apply(draws, 2, stats::quantile, c(.025, .975))
  parts <- data.frame(part = c("INTERCEPT", "SLOPE", "LINK", "COV"), size = obs, chance = chance,
                      excess = obs - chance, lo = obs - q[2, ], hi = obs - q[1, ])
  ## the standard calibration measures
  ci_int <- stats::glm(y ~ 1, offset = eta, family = stats::binomial())
  ci_slo <- stats::glm(y ~ eta, family = stats::binomial())
  wald <- function(f, j) { e <- stats::coef(f)[j]; s <- sqrt(diag(stats::vcov(f)))[j]; c(e, e - 1.96 * s, e + 1.96 * s) }
  lo <- stats::loess(y ~ p, span = 0.75, degree = 2)
  ici <- mean(abs(stats::predict(lo) - p))
  list(parts = parts,
       calibration = rbind(intercept = wald(ci_int, 1), slope = wald(ci_slo, 2)),
       ici = 100 * ici)
}
