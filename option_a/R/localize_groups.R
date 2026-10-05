## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# localize_groups.R -- closed-testing localization of misfit over ORTHOGONAL groups (Option A).
#
# Each group asks one question and has its own Cauchy combination of score tests whose bases are confined to that
# group's part of the misfit (weighted-orthogonal to the earlier parts):
#   INTERCEPT  calibration in the large                      basis 1
#   SLOPE      calibration slope                             basis eta, orthogonal to 1
#   LINK       bends of the score-to-probability map         eta^2, eta^3 | Stukel's two terms | ns(eta, 4),
#                                                            each orthogonal to {1, eta}
#   COV        misfit within covariate patterns              squares and cubes | pairwise products | ns(x_j, 3),
#                                                            each orthogonal to {1, eta, ns(eta, 5)} (functions of eta)
# External validation (frozen predictions): all four groups; the reference is Monte Carlo, y* ~ Bernoulli(p), exact.
# In-sample checking (a fitted glm): INTERCEPT and SLOPE are identically zero after the fit, so only LINK and COV;
# bases are also projected off the model matrix; the reference is the parametric bootstrap with refits.
# A group is named when every intersection containing it rejects (closure); each intersection statistic is the mean
# Cauchy coordinate of its members' chi-square p-values, calibrated by its rank among the reference draws.

.cauchy <- function(p) tan((0.5 - pmin(pmax(p, 1e-12), 1 - 1e-12)) * pi)

## weighted residualization of the columns of B against Z (w = p(1-p))
.wperp <- function(B, Z, w) {
  if (is.null(Z) || !ncol(Z)) return(B)
  sw <- sqrt(w)
  ## weighted least-squares residuals through a pivoted QR, so redundant columns of Z (e.g. eta inside the span of
  ## an intercept plus a spline in eta) are dropped instead of making the projection singular
  qr.resid(qr(Z * sw, tol = 1e-9), B * sw) / sw
}

## a score-test member: returns a function of the residual matrix R (n x m) giving chi-square p-values (length m).
## The basis is projected with the model weights w and its variance is the model variance sum(w b b').
## Multiplier calibration (r0 = observed residuals given): the variance is robust, sum(e^2 b b'), with e the residuals
## after refitting the projected-off directions Z as an offset model (e = r0 - w Z theta), so that it estimates the
## variance of the empirically projected statistic; R then holds r0 in column 1 and Rademacher signs in the others,
## and the reference draws are sum(b * sign * e). The observed statistic is unchanged, since sum(b w Z) = 0.
.member <- function(B, Z, w, r0 = NULL) {
  Bp <- .wperp(as.matrix(B), Z, w)
  ## drop columns the projection emptied (judged by their norm: the intercept's column is constant but not empty)
  keep <- sqrt(colSums(Bp^2)) > 1e-8 * (1 + sqrt(colSums(as.matrix(B)^2)))
  Bp <- Bp[, keep, drop = FALSE]
  if (!ncol(Bp)) return(NULL)
  ee <- if (is.null(r0)) NULL else if (is.null(Z) || !ncol(Z)) r0 else w * as.numeric(.wperp(matrix(r0 / w), Z, w))
  V  <- if (is.null(ee)) crossprod(Bp, Bp * w) else crossprod(Bp, Bp * ee^2)
  e  <- eigen(V, symmetric = TRUE)
  k  <- sum(e$values > 1e-9 * max(e$values))
  Vi <- e$vectors[, 1:k, drop = FALSE] %*% (t(e$vectors[, 1:k, drop = FALSE]) / e$values[1:k])
  function(R) {
    U <- if (is.null(ee)) crossprod(Bp, R) else cbind(crossprod(Bp, R[, 1]), crossprod(Bp, R[, -1, drop = FALSE] * ee))
    stats::pchisq(colSums(U * (Vi %*% U)), k, lower.tail = FALSE)
  }
}

## the group members on a design: eta and the covariate matrix X (numeric, named columns)
.group_bases <- function(eta, X) {
  num <- vapply(seq_len(ncol(X)), function(j) length(unique(X[, j])) > 4, TRUE)
  Xn  <- X[, num, drop = FALSE]
  cs  <- function(v) (v - mean(v)) / stats::sd(v)
  sq  <- if (ncol(Xn)) do.call(cbind, lapply(seq_len(ncol(Xn)), function(j) { z <- cs(Xn[, j]); cbind(z^2, z^3) })) else NULL
  pr  <- if (ncol(X) > 1) do.call(cbind, utils::combn(ncol(X), 2, function(ij) cs(X[, ij[1]]) * cs(X[, ij[2]]), simplify = FALSE)) else NULL
  sp  <- if (ncol(Xn)) do.call(cbind, lapply(seq_len(ncol(Xn)), function(j) splines::ns(Xn[, j], df = 3))) else NULL
  ## with a single covariate every function of it is a function of eta, so there is no covariate-structure part
  ## (Remark A.1): the COV group is dropped instead of testing the spline's approximation error
  list(
    LINK = list(poly = cbind(eta^2, eta^3), stukel = cbind(0.5 * eta^2 * (eta >= 0), -0.5 * eta^2 * (eta < 0)),
                spline = splines::ns(eta, df = 4)),
    COV  = if (ncol(X) < 2) list() else Filter(Negate(is.null), list(poly = sq, products = pr, spline = sp)))
}

## the spline dimension used for "functions of eta" when projecting the COV bases: fixed (default 5) or growing with
## n ("auto", Proposition A.6: the approximation error must vanish faster than n^(-1/2))
.cov_df <- function(cov_df, n) if (identical(cov_df, "auto")) max(5L, as.integer(ceiling(n^(1 / 3)))) else as.integer(cov_df)

## build every member function for one design; Zfit = the model matrix to project off (in-sample) or NULL (external).
## dealias_cols (in-sample only, de-aliased link group of Web Appendix A): model columns whose conditional mean given
## the score is curved; the LINK bases are also projected off their spline fit on the score, so that in-sample LINK
## reads link shape only. The columns are chosen once on the observed fit and frozen for the bootstrap draws.
.build_members <- function(eta, X, w, Zfit = NULL, external = TRUE, r0 = NULL, cov_df = 5, dealias_E = NULL) {
  gb <- .group_bases(eta, X)
  one <- matrix(1, length(eta), 1)
  Z2  <- cbind(one, eta, Zfit)
  Z3  <- cbind(Z2, splines::ns(eta, df = .cov_df(cov_df, length(eta))))
  ## de-aliasing: dealias_E holds the spline fits of the selected model columns on the OBSERVED score, frozen, so the
  ## same directions are removed in the data and in every bootstrap draw
  if (!is.null(dealias_E) && ncol(dealias_E)) Z2 <- cbind(Z2, dealias_E)
  mem <- list()
  if (external) {
    mem$INTERCEPT <- list(cal = .member(one, NULL, w, r0))
    mem$SLOPE     <- list(slope = .member(matrix(eta), one, w, r0))
  }
  mem$LINK <- Filter(Negate(is.null), lapply(gb$LINK, function(B) .member(B, Z2, w, r0)))
  mem$COV  <- Filter(Negate(is.null), lapply(gb$COV,  function(B) .member(B, Z3, w, r0)))
  Filter(length, mem)
}

## columns of the model matrix whose weighted mean given the score is curved: weighted F test of a 3-df natural spline
## in eta against a straight line, level .05
.curved_cols <- function(Zfit, eta, w, level = 0.05) {
  n <- length(eta); one <- matrix(1, n, 1)
  L <- cbind(one, eta); S <- cbind(one, splines::ns(eta, df = 3))
  which(vapply(seq_len(ncol(Zfit)), function(j) {
    z <- Zfit[, j]
    r1 <- sum(w * .wperp(matrix(z), L, w)^2); r2 <- sum(w * .wperp(matrix(z), S, w)^2)
    if (r2 <= 0) return(FALSE)
    stats::pf(((r1 - r2) / 2) / (r2 / (n - 4)), 2, n - 4, lower.tail = FALSE) < level
  }, TRUE))
}

## closure over the groups, given member p-values for the observed data (column 1) and the reference draws.
## naming = "closure": a group is named when every intersection containing it rejects (Cauchy intersections);
## "holm" / "bonferroni": the same groups' own (single-group) p-values with Holm's or Bonferroni's correction, which are
## closed procedures with Bonferroni intersections and so share the strong familywise guarantee
.closure <- function(P, groups, alpha, naming = "closure") {
  names_g <- names(groups)
  pint <- list()
  for (r in seq_along(names_g)) for (S in utils::combn(names_g, r, simplify = FALSE)) {
    rows <- unlist(groups[S], use.names = FALSE)
    stat <- colMeans(.cauchy(P[rows, , drop = FALSE]))
    pint[[paste(S, collapse = "+")]] <- (1 + sum(stat[-1] >= stat[1])) / length(stat)
  }
  single <- unlist(pint[names_g])
  named <- switch(naming,
    closure    = vapply(names_g, function(g) all(unlist(pint[grepl(paste0("(^|\\+)", g, "($|\\+)"), names(pint))]) <= alpha), TRUE),
    holm       = stats::p.adjust(single, "holm") <= alpha,
    bonferroni = stats::p.adjust(single, "bonferroni") <= alpha)
  list(named = names_g[named], intersection = unlist(pint), naming = naming)
}

## normal scores, the robust option: each continuous column and the score are replaced by qnorm(rank / (n + 1)), so that
## no record can carry more than a bounded share of any statistic however extreme its covariate values are; columns
## with four or fewer distinct values (binary, small ordinal) are left as they are
.nscore <- function(v) stats::qnorm(rank(v, ties.method = "average") / (length(v) + 1))
.robust_cols <- function(X) {
  X <- as.matrix(X)
  for (j in seq_len(ncol(X))) if (length(unique(X[, j])) > 4) X[, j] <- .nscore(X[, j])
  X
}

## EXTERNAL validation of frozen predictions p (from any model) for outcomes y, covariates X.
## calibration = "montecarlo": y* ~ Bernoulli(p), exact under no misfit (Theorem A.3, local alternatives);
## calibration = "multiplier": robust variance and sign-flipped residuals, so that a group whose part of the misfit is
## zero on the probability scale keeps its level for departures of any size (Theorem A.5), asymptotically
localize_external <- function(y, p, X, M = 999L, alpha = 0.05, calibration = c("montecarlo", "multiplier"), cov_df = 5,
                              naming = c("closure", "holm", "bonferroni"), robust = FALSE) {
  calibration <- match.arg(calibration); naming <- match.arg(naming)
  ## predictions of exactly 0 or 1 would give an infinite score; keep them just inside (0, 1)
  X <- as.matrix(X); p <- pmin(pmax(p, 1e-10), 1 - 1e-10); eta <- stats::qlogis(p); w <- p * (1 - p)
  if (robust) { eta <- .nscore(eta); X <- .robust_cols(X) }   # bases built on normal scores (bounded influence)
  r0 <- y - p
  mem <- .build_members(eta, X, w, external = TRUE, r0 = if (calibration == "multiplier") r0 else NULL, cov_df = cov_df)
  R <- if (calibration == "montecarlo") {
    cbind(r0, vapply(seq_len(M), function(m) stats::rbinom(length(p), 1L, p) - p, numeric(length(p))))
  } else {
    ## Rademacher signs, shared by all members so that the draws keep the members' joint dependence
    cbind(r0, matrix(sample(c(-1, 1), length(p) * M, replace = TRUE), length(p), M))
  }
  P <- do.call(rbind, lapply(unlist(mem, recursive = FALSE), function(f) f(R)))
  groups <- lapply(names(mem), function(g) grep(paste0("^", g, "\\."), rownames(P)))
  names(groups) <- names(mem)
  out <- .closure(P, groups, alpha, naming)
  out$members <- P[, 1]; out$setting <- "external validation"; out$calibration <- calibration; out$robust <- robust; out
}

## IN-SAMPLE checking of a fitted binomial glm (logit link); covariates = the model frame's numeric predictors.
## robust = TRUE builds the bases on normal scores of the covariates, the model columns and the fitted score; the
## bootstrap refits the model on the original design, so the level does not depend on the choice of bases
localize_insample <- function(fit, X, B = 199L, alpha = 0.05, dealias = FALSE, cov_df = 5,
                              naming = c("closure", "holm", "bonferroni"), robust = FALSE) {
  naming <- match.arg(naming)
  X <- as.matrix(X); mm <- stats::model.matrix(fit)
  Xb <- if (robust) .robust_cols(X) else X
  Zb <- if (robust) .robust_cols(mm[, -1, drop = FALSE]) else mm[, -1, drop = FALSE]
  sc <- function(ph) if (robust) .nscore(stats::qlogis(ph)) else stats::qlogis(ph)
  ph0 <- as.numeric(stats::fitted(fit))
  dcols <- if (dealias) .curved_cols(Zb, sc(ph0), ph0 * (1 - ph0)) else NULL
  dE <- NULL
  if (length(dcols)) {   # fitted once, on the observed score, and frozen
    e0 <- sc(ph0); w0 <- ph0 * (1 - ph0)
    Zc <- Zb[, dcols, drop = FALSE]
    dE <- Zc - .wperp(Zc, cbind(1, splines::ns(e0, df = 3)), w0)
  }
  members_at <- function(f) {
    ph <- as.numeric(stats::fitted(f))
    .build_members(sc(ph), Xb, ph * (1 - ph), Zfit = Zb, external = FALSE, cov_df = cov_df, dealias_E = dE)
  }
  pv_of <- function(f, y) { m <- members_at(f); r <- matrix(y - as.numeric(stats::fitted(f)))
                            vapply(unlist(m, recursive = FALSE), function(g) g(r), 0) }
  y0 <- as.numeric(fit$y); ph <- as.numeric(stats::fitted(fit))
  P <- matrix(pv_of(fit, y0), ncol = 1)
  rn <- names(pv_of(fit, y0))
  for (b in seq_len(B)) {
    yb <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(mm, yb, family = stats::binomial())), error = function(e) NULL)
    if (is.null(fb)) next
    class(fb) <- c("glm", "lm"); fb$y <- yb
    P <- cbind(P, tryCatch(pv_of(fb, yb), error = function(e) rep(NA_real_, length(rn))))
  }
  P <- P[, colSums(is.na(P)) == 0, drop = FALSE]
  rownames(P) <- rn
  gn <- Filter(function(g) any(grepl(paste0("^", g, "\\."), rn)), c("LINK", "COV"))   # COV is absent with one covariate
  groups <- lapply(gn, function(g) grep(paste0("^", g, "\\."), rn)); names(groups) <- gn
  out <- .closure(P, groups, alpha, naming)
  out$members <- P[, 1]; out$B_used <- ncol(P) - 1L; out$setting <- "in-sample checking"; out$robust <- robust
  out$dealiased <- colnames(mm)[-1][dcols]; out
}
