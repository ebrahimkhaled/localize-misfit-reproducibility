# run_sim4.R -- the frozen model has ESTIMATED coefficients. The designs of run_sim3.R, with one change: the frozen model is a
# main-effects logistic model fitted on a development sample of 500 or 2,000 patients from the baseline population
# (true logit a + x'b), so it carries development error along its covariates and no update on the ladder is exactly right.
# Validation data and new patients follow the settings of run_sim3.R, defined on the TRUE score: none (the same population),
# transport, curved link, omitted interaction. Prevalence 20%; validation sample of 100, 200 or 400 expected events.
# Per data set: the verdicts of the three rules, and the error of the six candidate updates for 20,000 new patients.
# Output sim4_rows/<cell>_r<first>.csv, ten data sets per file, written as they finish.
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
OUT <- file.path(HERE, "sim4_rows"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

DESIGNS <- c("indep5", "mixed10"); PREV <- 0.20; EVENTS <- c(100L, 200L, 400L); NDEV <- c(500L, 2000L)
DEPS <- c("null", "transport", "link", "inter")
CELLS <- expand.grid(dep = DEPS, events = EVENTS, ndev = NDEV, design = DESIGNS, stringsAsFactors = FALSE)
CELLS$prev <- PREV; CELLS$id <- seq_len(nrow(CELLS)); CELLS$reps <- if (SMOKE) 1L else 200L
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 10))) {
  f <- file.path(OUT, sprintf("%s_d%04d_e%03d_%s_r%04d.csv", CELLS$design[q], CELLS$ndev[q], CELLS$events[q], CELLS$dep[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

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
intercept_for <- function(design, prev) {
  set.seed(20261105L); X <- gen_X(design, 200000L); lp <- as.numeric(X %*% BETA[[design]])
  stats::uniroot(function(a) mean(stats::plogis(a + lp)) - prev, c(-10, 5))$root
}
A0 <- stats::setNames(vapply(DESIGNS, intercept_for, 0, prev = PREV), DESIGNS)

true_logit <- function(dep, eta, X, a) {
  z <- (eta - mean(eta)) / stats::sd(eta)
  switch(dep, null = eta, transport = eta - 1 + (0.8 - 1) * (eta - a), link = eta + 0.25 * (z^2 - 1), inter = eta + 0.6 * X[, 1] * X[, 3])
}

rev_basis <- function(Xv) {
  cont <- which(apply(Xv, 2, function(v) length(unique(v)) > 4))
  mu <- colMeans(Xv); sdv <- apply(Xv, 2, stats::sd)
  kn <- lapply(cont, function(j) splines::ns(Xv[, j], df = 3))
  function(Xn) {
    Zs <- sweep(sweep(Xn, 2, mu), 2, sdv, "/")
    sq <- do.call(cbind, lapply(cont, function(j) cbind(Zs[, j]^2, Zs[, j]^3)))
    pr <- do.call(cbind, utils::combn(ncol(Xn), 2, function(ij) Zs[, ij[1]] * Zs[, ij[2]], simplify = FALSE))
    sp <- do.call(cbind, lapply(seq_along(cont), function(k) stats::predict(kn[[k]], Xn[, cont[k]])))
    cbind(Zs, sq, pr, sp)
  }
}

vergouwe_choice <- function(lp_old, X, y, alpha = 0.05) {
  ll <- function(lp) sum(y * lp - log1p(exp(lp)))
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
  cc <- j$cc; a <- A0[[cc$design]]; b <- BETA[[cc$design]]
  n <- as.integer(ceiling(cc$events / cc$prev)); nt <- 20000L
  G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
  ACT <- c(INTERCEPT = "intercept", SLOPE = "recalibration", LINK = "flexible", COV = "revision_cov")
  R <- do.call(rbind, lapply(j$reps, function(r) {
    set.seed(56000000L + 10000L * cc$id + r)
    ## the development sample and the frozen model
    Xd <- gen_X(cc$design, cc$ndev); yd <- stats::rbinom(cc$ndev, 1L, stats::plogis(as.numeric(a + Xd %*% b)))
    bh <- stats::coef(stats::glm(yd ~ Xd, family = stats::binomial()))
    score <- function(Xn) as.numeric(bh[1] + Xn %*% bh[-1])
    ## validation data and new patients: truth on the true score, predictions from the frozen model
    X <- gen_X(cc$design, n); y <- stats::rbinom(n, 1L, stats::plogis(true_logit(cc$dep, as.numeric(a + X %*% b), X, a)))
    Xt <- gen_X(cc$design, nt); pit <- stats::plogis(true_logit(cc$dep, as.numeric(a + Xt %*% b), Xt, a))
    eta <- score(X); etat <- score(Xt)
    loc <- localize_external_update(y, stats::plogis(eta), X, B = 199L)
    top <- if (length(loc$named)) ACT[[G[max(match(loc$named, G))]]] else "none"
    ver <- vergouwe_choice(eta, X, y)
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
    data.frame(design = cc$design, ndev = cc$ndev, events = cc$events, dep = cc$dep, rep = r, n = n, y_events = sum(y),
               loc = top, vergouwe = ver, aic = aic,
               t(stats::setNames(G %in% loc$named, paste0("named_", G))),
               t(stats::setNames(as.numeric(loc$intersection[G]), paste0("p_", G))),
               t(stats::setNames(err, paste0("err_", names(err)))), t(stats::setNames(lls, paste0("ll_", names(lls)))))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}

cat(sprintf("sim4%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  jobs <- jobs[order(-vapply(jobs, function(j) ceiling(j$cc$events / j$cc$prev), 0))]
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "16")))
  parallel::clusterExport(cl, c("HERE", "one_chunk", "A0", "BETA", "gen_X", "true_logit", "rev_basis", "vergouwe_choice"))
  res <- parallel::parLapplyLB(cl, jobs, function(j) tryCatch(one_chunk(j), error = function(e) paste("ERROR", j$f, conditionMessage(e))),
                               chunk.size = 1L)
  parallel::stopCluster(cl)
  bad <- grep("^ERROR", unlist(res), value = TRUE); cat(length(bad), "failed chunks\n"); if (length(bad)) print(utils::head(bad))
}
cat("done", format(Sys.time()), "\n")
