# localize_update.R -- external validation with each part tested against the model ALREADY UPDATED for the simpler parts
# ("updated-model reference"). The test statistics are those of localize_external(); only the reference changes:
#   INTERCEPT  against the frozen model p:                    Monte Carlo, y* ~ Bernoulli(p)
#   SLOPE      against the intercept update, a + eta0:         parametric bootstrap from the fitted update, refitted
#   LINK       against the logistic recalibration, a + b eta0: parametric bootstrap, refitted
#   COV        against the flexible recalibration, ns(eta0,4): parametric bootstrap, refitted
# The update is refitted to every draw and each part's bases are projected off that update's own directions with its
# own fitted weights, so the null of each test is "this part is zero, whatever the simpler parts are": a large error in
# overall risk or slope is absorbed by the update and cannot inflate the tests of the more demanding parts. This is the
# parametric bootstrap of the in-sample procedure (Web Appendix A, Lemma A.boot) with the update as the fitted model.
# With cluster and tau the draws carry a cluster random effect, so that a between-centre variance is part of the reference.
# The four parts are named by Holm's correction over their p-values (default) or by Bonferroni.
localize_external_update <- function(y, p, X, B = 199L, alpha = 0.05, naming = c("holm", "bonferroni"), cov_df = 5,
                                     cluster = NULL, tau = 0, main = TRUE) {
  naming <- match.arg(naming)
  X <- as.matrix(X); p <- pmin(pmax(p, 1e-10), 1 - 1e-10); eta0 <- stats::qlogis(p); n <- length(y)
  one <- matrix(1, n, 1)
  ## reference draws; with clusters, each draw adds a new random effect per cluster, N(0, tau^2) on the logit scale
  rdraw <- function(ph) {
    if (is.null(cluster) || tau <= 0) return(stats::rbinom(n, 1L, ph))
    cl <- as.integer(factor(cluster))
    stats::rbinom(n, 1L, stats::plogis(stats::qlogis(ph) + stats::rnorm(max(cl), 0, tau)[cl]))
  }
  cauchy_p <- function(P) { s <- colMeans(.cauchy(P)); (1 + sum(s[-1] >= s[1])) / length(s) }
  gb <- .group_bases(eta0, X)
  ## a fourth COV member: the main effects, so that a changed coefficient (of a binary covariate too) has a test of its own
  ## main = FALSE for pooled cross-validated risks, which are not frozen: their coefficients differ between folds, and a main-effects
  ## test sees that difference even when the model's form is right
  if (main && length(gb$COV)) gb$COV <- c(list(main = X), gb$COV)
  N4 <- splines::ns(eta0, df = 4)
  ## each part: the update it is tested against (fitted to an outcome vector), and its members at that fit
  parts <- list(
    SLOPE = list(fit = function(yy) stats::fitted(stats::glm(yy ~ 1, offset = eta0, family = stats::binomial())),
                 members = function(ph) list(.member(matrix(eta0), one, ph * (1 - ph)))),
    LINK = list(fit = function(yy) stats::fitted(stats::glm(yy ~ eta0, family = stats::binomial())),
                members = function(ph) lapply(gb$LINK, function(Bm) .member(Bm, cbind(one, eta0), ph * (1 - ph)))),
    COV = list(fit = function(yy) stats::fitted(stats::glm(yy ~ N4, family = stats::binomial())),
               members = function(ph) {
                 el <- stats::qlogis(ph)
                 Z <- cbind(one, el, splines::ns(el, df = .cov_df(cov_df, n)), N4)
                 lapply(gb$COV, function(Bm) .member(Bm, Z, ph * (1 - ph)))
               }))
  if (!length(gb$COV)) parts$COV <- NULL
  ## the member p-values of one part for one outcome vector (the update refitted to it)
  pv <- function(part, yy) {
    ph <- as.numeric(part$fit(yy)); m <- Filter(Negate(is.null), part$members(ph))
    vapply(m, function(f) f(matrix(yy - ph)), 0)
  }
  single <- c(INTERCEPT = {
    w <- p * (1 - p); f <- .member(one, NULL, w)
    R <- cbind(y - p, vapply(seq_len(B), function(b) rdraw(p) - p, numeric(n)))
    cauchy_p(matrix(f(R), nrow = 1))
  })
  for (g in names(parts)) {
    pr <- parts[[g]]; ph <- as.numeric(pr$fit(y)); p_obs <- pv(pr, y); k <- length(p_obs)
    if (anyNA(p_obs)) stop("the observed test of ", g, " could not be computed")
    draws <- vapply(seq_len(B), function(b) {
      yb <- rdraw(ph)
      tryCatch(pv(pr, yb), error = function(e) rep(NA_real_, k))
    }, numeric(k))
    ## keep a members-by-draws matrix also when the part has a single member (vapply then returns a vector)
    P <- cbind(p_obs, matrix(draws, nrow = k))
    P <- P[, colSums(is.na(P)) == 0, drop = FALSE]
    single[g] <- cauchy_p(P)
  }
  adj <- stats::p.adjust(single, naming)
  out <- list(named = names(single)[adj <= alpha], intersection = single, adjusted = adj, naming = naming,
              setting = "external validation, updated-model reference", alpha = alpha, n = n)
  class(out) <- "gof_localize"; out
}
