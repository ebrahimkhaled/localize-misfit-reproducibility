# run_sim2.R -- the main simulation of the SiM version: realistic external-validation designs and the CONSEQUENCES of
# following each procedure's verdict for new patients.
# Designs: covariates "indep5" (five independent N(0,1)) or "mixed10" (five correlated normals, rho = 0.5; three binary,
# P = 0.3; two right-skewed, standardized lognormal); prevalence 5%, 20%, 40% (intercept set so that the frozen model
# has that mean risk); validation sample of 100 or 200 expected events (n = events / prevalence).
# Departures (true logit relative to the frozen score eta0 = a + x'b): none; intercept +0.4; slope 0.7 (the frozen
# model is overfitted: its risks too extreme); transport (intercept -1 and slope 0.8); curved link 0.25 (z^2 - 1), z the
# standardized score; omitted interaction 0.6 x1 x3; U-shaped effect 0.5 (x1^2 - 1).
# Per data set: (1) the verdicts of the localization procedure (updated-model reference, Holm; R/localize_update.R), the
# closed testing procedure of Vergouwe et al. (2017), and AIC over the updating ladder; (2) six candidate updates fitted
# on the validation data -- none, intercept update, logistic recalibration, flexible recalibration (ns(eta0, 3)), linear
# revision (all coefficients refitted, as in Vergouwe et al.), and covariate revision (the frozen score as an offset plus
# the nonlinear and interaction terms of the COV group, ridge-penalized, penalty by 5-fold cross-validation) -- and their
# error on 20,000 new patients from the same population: mean absolute error against the TRUE risk (percentage points)
# and expected log-loss. Output sim2_rows/<cell>_r<first>.csv, ten data sets per file, written as they finish.
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
OUT <- file.path(HERE, "sim2_rows"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

DESIGNS <- c("indep5", "mixed10"); PREV <- c(0.05, 0.20, 0.40); EVENTS <- c(100L, 200L)
DEPS <- c("null", "intercept", "slope", "transport", "link", "inter", "ushape")
CELLS <- expand.grid(dep = DEPS, events = EVENTS, prev = PREV, design = DESIGNS, stringsAsFactors = FALSE)
CELLS$id <- seq_len(nrow(CELLS))
CELLS$reps <- ifelse(CELLS$dep == "null", 500L, 200L)
if (SMOKE) CELLS$reps <- 1L
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 10))) {
  f <- file.path(OUT, sprintf("%s_p%02d_e%03d_%s_r%04d.csv", CELLS$design[q], round(100 * CELLS$prev[q]), CELLS$events[q],
                              CELLS$dep[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

## covariates and coefficients of the two designs
gen_X <- function(design, n) {
  if (design == "indep5") {
    X <- matrix(stats::rnorm(n * 5), n, 5)
  } else {
    Z <- matrix(stats::rnorm(n * 5), n, 5) %*% chol(0.5 + 0.5 * diag(5))
    B <- matrix(stats::rbinom(n * 3, 1L, 0.3), n, 3)
    L <- exp(matrix(stats::rnorm(n * 2, 0, 0.75), n, 2)); L <- (L - exp(0.75^2 / 2)) / sqrt((exp(0.75^2) - 1) * exp(0.75^2))
    X <- cbind(Z, B, L)
  }
  colnames(X) <- paste0("x", seq_len(ncol(X))); X
}
BETA <- list(indep5 = c(0.5, 0.8, 0.6, 0.4, 0.3), mixed10 = c(0.5, 0.6, 0.4, 0.3, 0.2, 0.5, 0.4, -0.3, 0.4, 0.3))
## intercept giving the target mean risk of the frozen model (computed once per design and prevalence)
intercept_for <- function(design, prev) {
  set.seed(20261105L); X <- gen_X(design, 200000L); lp <- as.numeric(X %*% BETA[[design]])
  stats::uniroot(function(a) mean(stats::plogis(a + lp)) - prev, c(-10, 5))$root
}
A0 <- with(expand.grid(design = DESIGNS, prev = PREV, stringsAsFactors = FALSE),
           stats::setNames(mapply(intercept_for, design, prev), paste(design, prev)))

true_logit <- function(dep, eta, X, a) {
  z <- (eta - mean(eta)) / stats::sd(eta)
  switch(dep, null = eta, intercept = eta + 0.4, slope = a + 0.7 * (eta - a), transport = eta - 1 + (0.8 - 1) * (eta - a),
         link = eta + 0.25 * (z^2 - 1), inter = eta + 0.6 * X[, 1] * X[, 3], ushape = eta + 0.5 * (X[, 1]^2 - 1))
}

## the nonlinear and interaction terms of the covariate revision, with the transformation fixed on the validation data
rev_basis <- function(Xv) {
  cont <- which(apply(Xv, 2, function(v) length(unique(v)) > 4))
  mu <- colMeans(Xv); sdv <- apply(Xv, 2, stats::sd)
  kn <- lapply(cont, function(j) splines::ns(Xv[, j], df = 3))
  function(Xn) {
    Zs <- sweep(sweep(Xn, 2, mu), 2, sdv, "/")
    sq <- do.call(cbind, lapply(cont, function(j) cbind(Zs[, j]^2, Zs[, j]^3)))
    pr <- do.call(cbind, utils::combn(ncol(Xn), 2, function(ij) Zs[, ij[1]] * Zs[, ij[2]], simplify = FALSE))
    sp <- do.call(cbind, lapply(seq_along(cont), function(k) stats::predict(kn[[k]], Xn[, cont[k]])))
    cbind(sq, pr, sp)
  }
}

vergouwe_choice <- function(a, b, X, y, alpha = 0.05) {
  lp_old <- as.numeric(a + X %*% b); ll <- function(lp) sum(y * lp - log1p(exp(lp)))
  int <- stats::coef(stats::glm(y ~ 1, offset = lp_old, family = stats::binomial()))
  rc <- stats::coef(stats::glm(y ~ lp_old, family = stats::binomial()))
  rv <- suppressWarnings(stats::glm(y ~ X, family = stats::binomial()))
  l0 <- ll(lp_old); l1 <- ll(int + lp_old); l2 <- ll(rc[1] + rc[2] * lp_old); l3 <- as.numeric(stats::logLik(rv)); p <- ncol(X)
  if (stats::pchisq(2 * (l3 - l0), p + 1, lower.tail = FALSE) >= alpha) "none"
  else if (stats::pchisq(2 * (l3 - l1), p, lower.tail = FALSE) >= alpha) "intercept"
  else if (stats::pchisq(2 * (l3 - l2), p - 1, lower.tail = FALSE) >= alpha) "recalibration" else "revision_linear"
}

one_chunk <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  cc <- j$cc; a <- A0[[paste(cc$design, cc$prev)]]; b <- BETA[[cc$design]]
  n <- as.integer(ceiling(cc$events / cc$prev)); nt <- 20000L
  G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
  ACT <- c(INTERCEPT = "intercept", SLOPE = "recalibration", LINK = "flexible", COV = "revision_cov")
  R <- do.call(rbind, lapply(j$reps, function(r) {
    set.seed(54000000L + 10000L * cc$id + r)
    X <- gen_X(cc$design, n); eta <- as.numeric(a + X %*% b); y <- stats::rbinom(n, 1L, stats::plogis(true_logit(cc$dep, eta, X, a)))
    Xt <- gen_X(cc$design, nt); etat <- as.numeric(a + Xt %*% b); pit <- stats::plogis(true_logit(cc$dep, etat, Xt, a))
    t0 <- proc.time()[["elapsed"]]
    loc <- localize_external_update(y, stats::plogis(eta), X, B = 199L)
    t_loc <- proc.time()[["elapsed"]] - t0
    top <- if (length(loc$named)) ACT[[G[max(match(loc$named, G))]]] else "none"
    ver <- vergouwe_choice(a, b, X, y)
    ## the six candidate updates, fitted on the validation data, predicted for the new patients
    N3 <- splines::ns(eta, df = 3)
    fI <- stats::glm(y ~ 1, offset = eta, family = stats::binomial())
    fS <- stats::glm(y ~ eta, family = stats::binomial())
    fL <- stats::glm(y ~ N3, family = stats::binomial())
    fR <- suppressWarnings(stats::glm(y ~ X, family = stats::binomial()))
    rb <- rev_basis(X); Cv <- rb(X); Ct <- rb(Xt)
    cvf <- glmnet::cv.glmnet(Cv, y, family = "binomial", alpha = 0, offset = eta, nfolds = 5)
    pred <- list(none = stats::plogis(etat),
                 intercept = stats::plogis(stats::coef(fI)[1] + etat),
                 recalibration = stats::plogis(stats::coef(fS)[1] + stats::coef(fS)[2] * etat),
                 flexible = stats::plogis(as.numeric(cbind(1, stats::predict(N3, etat)) %*% stats::coef(fL))),
                 revision_linear = stats::plogis(as.numeric(cbind(1, Xt) %*% stats::coef(fR))),
                 revision_cov = as.numeric(stats::predict(cvf, Ct, newoffset = etat, s = "lambda.min", type = "response")))
    ll <- c(none = 0, intercept = 1, recalibration = 2, flexible = 4, revision_linear = ncol(X) + 1)
    lv <- c(none = sum(stats::dbinom(y, 1, stats::plogis(eta), log = TRUE)), intercept = as.numeric(stats::logLik(fI)),
            recalibration = as.numeric(stats::logLik(fS)), flexible = as.numeric(stats::logLik(fL)),
            revision_linear = as.numeric(stats::logLik(fR)))
    aic <- names(which.min(-2 * lv + 2 * ll))
    err <- vapply(pred, function(q) 100 * mean(abs(q - pit)), 0)
    lls <- vapply(pred, function(q) { q <- pmin(pmax(q, 1e-12), 1 - 1e-12); -mean(pit * log(q) + (1 - pit) * log1p(-q)) }, 0)
    data.frame(design = cc$design, prev = cc$prev, events = cc$events, dep = cc$dep, rep = r, n = n, y_events = sum(y),
               loc = top, vergouwe = ver, aic = aic,
               t(stats::setNames(G %in% loc$named, paste0("named_", G))),
               t(stats::setNames(as.numeric(loc$intersection[G]), paste0("p_", G))),
               t(stats::setNames(err, paste0("err_", names(err)))), t(stats::setNames(lls, paste0("ll_", names(lls)))),
               secs_loc = t_loc)
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}

cat(sprintf("sim2%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  ## the biggest validation sets first, so that the slowest chunks do not finish last
  jobs <- jobs[order(-vapply(jobs, function(j) ceiling(j$cc$events / j$cc$prev), 0))]
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "18")))
  parallel::clusterExport(cl, c("HERE", "one_chunk", "A0", "BETA", "gen_X", "true_logit", "rev_basis", "vergouwe_choice"))
  invisible(parallel::parLapplyLB(cl, jobs, function(j) tryCatch(one_chunk(j), error = function(e) paste("ERROR", j$f, conditionMessage(e))),
                                  chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
