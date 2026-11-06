# run_davis.R -- two further rivals on the SAME data sets as the main simulation (run_sim2.R) and the focused comparison
# (run_rivals_update.R); every data set is checked against its stored event count.
#   davis     the nonparametric testing procedure of Davis et al. (2019, JAMIA 26:1448-1457): candidates no update,
#             intercept correction, logistic recalibration, flexible recalibration, refit. Stage 1: B1 = 100 bootstrap
#             samples of the update set; each update is fitted in-bag and applied to the out-of-bag patients, whose
#             predictions are pooled. Stage 2: B2 = 1,000 bootstrap samples of the pooled set give each method's Brier
#             score. Mr is the method with the lowest median score; going from the simplest method up, the first method
#             whose 95% percentile interval for the paired difference with Mr contains 0 is recommended.
#   ctp_flex  closed testing of Vergouwe et al. (2017) with a flexible step, our implementation: likelihood-ratio tests of
#             the largest model (all coefficients re-estimated plus a spline in the score) against the original, the
#             intercept update, the logistic recalibration and the flexible recalibration, in that order, each at 0.05.
# Output sim2d_rows/ (main study, same file names as sim2_rows/) and rows_davis/ (focused comparison).
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
OUT1 <- file.path(HERE, "sim2d_rows"); dir.create(OUT1, showWarnings = FALSE)
OUT2 <- file.path(HERE, "rows_davis"); dir.create(OUT2, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

## main study: designs exactly as run_sim2.R
DESIGNS <- c("indep5", "mixed10"); PREV <- c(0.05, 0.20, 0.40); EVENTS <- c(100L, 200L)
DEPS <- c("null", "intercept", "slope", "transport", "link", "inter", "ushape")
CELLS <- expand.grid(dep = DEPS, events = EVENTS, prev = PREV, design = DESIGNS, stringsAsFactors = FALSE)
CELLS$id <- seq_len(nrow(CELLS))
## the setting added later keeps the ids of all earlier cells, so their seeds and data sets are unchanged
NEWC <- expand.grid(dep = "effect", events = EVENTS, prev = PREV, design = DESIGNS, stringsAsFactors = FALSE)
NEWC$id <- nrow(CELLS) + seq_len(nrow(NEWC)); CELLS <- rbind(CELLS, NEWC)
CELLS$reps <- ifelse(CELLS$dep == "null", 500L, 200L)
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
A0 <- with(expand.grid(design = DESIGNS, prev = PREV, stringsAsFactors = FALSE),
           stats::setNames(mapply(intercept_for, design, prev), paste(design, prev)))
true_logit <- function(dep, eta, X, a) {
  z <- (eta - mean(eta)) / stats::sd(eta)
  switch(dep, null = eta, intercept = eta + 0.4, slope = a + 0.7 * (eta - a), transport = eta - 1 + (0.8 - 1) * (eta - a),
         link = eta + 0.25 * (z^2 - 1), inter = eta + 0.6 * X[, 1] * X[, 3], ushape = eta + 0.5 * (X[, 1]^2 - 1),
         effect = if (ncol(X) == 5L) eta + 0.5 * X[, 1] else eta + 0.8 * X[, 6])
}

METHODS <- c("none", "intercept", "recalibration", "flexible", "revision_linear")
fitb <- function(Z, y, off = NULL) stats::glm.fit(Z, y, offset = off, family = stats::binomial())$coefficients

davis_choice <- function(y, eta, X, B1 = 100L, B2 = 1000L, alpha = 0.05) {
  n <- length(y); N <- splines::ns(eta, df = 3)
  Nb <- cbind(1, unclass(N)[, 1:3, drop = FALSE]); Xb <- cbind(1, X); one <- matrix(1, n, 1); Sb <- cbind(1, eta)
  Hy <- vector("list", B1); Hp <- vector("list", B1)
  for (b in seq_len(B1)) {
    u <- sample.int(n, n, replace = TRUE); h <- setdiff(seq_len(n), u)
    yu <- y[u]
    ci <- suppressWarnings(fitb(one[u, , drop = FALSE], yu, eta[u]))
    cs <- suppressWarnings(fitb(Sb[u, ], yu)); cf <- suppressWarnings(fitb(Nb[u, ], yu)); cr <- suppressWarnings(fitb(Xb[u, ], yu))
    cf[is.na(cf)] <- 0; cr[is.na(cr)] <- 0; cs[is.na(cs)] <- 0
    Hp[[b]] <- cbind(stats::plogis(eta[h]), stats::plogis(ci + eta[h]), stats::plogis(as.numeric(Sb[h, ] %*% cs)),
                     stats::plogis(as.numeric(Nb[h, ] %*% cf)), stats::plogis(as.numeric(Xb[h, ] %*% cr)))
    Hy[[b]] <- y[h]
  }
  P <- do.call(rbind, Hp); Y <- unlist(Hy); E <- (P - Y)^2
  S <- matrix(0, B2, 5); m <- nrow(E)
  for (b in seq_len(B2)) S[b, ] <- colMeans(E[sample.int(m, n, replace = TRUE), , drop = FALSE])
  r <- which.min(apply(S, 2, stats::median))
  for (k in seq_len(r)) {
    if (k == r) return(METHODS[r])
    ci <- stats::quantile(S[, k] - S[, r], c(alpha / 2, 1 - alpha / 2), names = FALSE)
    if (ci[1] <= 0 && ci[2] >= 0) return(METHODS[k])
  }
  METHODS[r]
}

ctp_flex_choice <- function(y, eta, X, alpha = 0.05) {
  ll <- function(lp) sum(y * lp - log1p(exp(lp)))
  N <- unclass(splines::ns(eta, df = 3))[, 1:3, drop = FALSE]
  l0 <- ll(eta); l1 <- ll(stats::coef(stats::glm(y ~ 1, offset = eta, family = stats::binomial())) + eta)
  l2 <- as.numeric(stats::logLik(stats::glm(y ~ eta, family = stats::binomial())))
  l3 <- as.numeric(stats::logLik(stats::glm(y ~ N, family = stats::binomial())))
  big <- suppressWarnings(stats::glm(y ~ X + N, family = stats::binomial()))
  l4 <- as.numeric(stats::logLik(big)); k4 <- big$rank
  pv <- function(l, k) stats::pchisq(2 * (l4 - l), k4 - k, lower.tail = FALSE)
  if (pv(l0, 0) >= alpha) "none" else if (pv(l1, 1) >= alpha) "intercept" else if (pv(l2, 2) >= alpha) "recalibration"
  else if (pv(l3, 4) >= alpha) "flexible" else "revision_linear"
}

main_chunk <- function(j) {
  suppressMessages(library(splines))
  cc <- j$cc; a <- A0[[paste(cc$design, cc$prev)]]; b <- BETA[[cc$design]]; n <- as.integer(ceiling(cc$events / cc$prev))
  st <- utils::read.csv(j$stored)
  R <- do.call(rbind, lapply(j$reps, function(r) {
    set.seed(54000000L + 10000L * cc$id + r)
    X <- gen_X(cc$design, n); eta <- as.numeric(a + X %*% b); y <- stats::rbinom(n, 1L, stats::plogis(true_logit(cc$dep, eta, X, a)))
    ev <- st$y_events[st$rep == r]
    if (length(ev) != 1L || as.integer(ev) != sum(y)) stop("data set mismatch ", j$f, " rep ", r)
    set.seed(59000000L + 10000L * cc$id + r)
    data.frame(design = cc$design, prev = cc$prev, events = cc$events, dep = cc$dep, rep = r,
               davis = davis_choice(y, eta, X), ctp_flex = ctp_flex_choice(y, eta, X))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}

## focused comparison: as run_rivals_update.R
ext <- data.frame(fam = c("null", "intercept", "intercept", "slope", "slope", "slope_pure", "link", "link_pure", "link_pure",
                          "ushape", "thresh", "inter", "cov_pure", "cov_pure"), C = c(0, .2, .4, .8, .6, .6, .25, .25, .5, .8, 4, 1, 1, 2))
ins <- data.frame(fam = c("null", "link", "link_pure", "link_pure", "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C = c(0, .25, .25, .5, .8, 4, 1, 1, 2))
FC <- rbind(merge(cbind(setting = "external", ext), data.frame(n = c(500L, 1000L))),
            merge(cbind(setting = "insample", ins), data.frame(n = c(500L, 1000L))))
FC$reps <- ifelse(FC$fam == "null", 1000L, 300L); FC$id <- seq_len(nrow(FC))
keep <- paste(c("null", "intercept", "slope", "link", "ushape", "thresh", "inter"), c(0, .4, .6, .25, .8, 4, 1))
FC <- FC[FC$setting == "external" & paste(FC$fam, FC$C) %in% keep & FC$n == 1000L, ]
focus_chunk <- function(j) {
  suppressMessages(library(splines))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  stored <- utils::read.csv(file.path(HERE, "rows", sprintf("external_%s_C%.2f_n%04d_r%04d.csv", j$cc$fam, j$cc$C, j$cc$n, j$reps[1])))
  R <- do.call(rbind, lapply(j$reps, function(r) {
    n <- j$cc$n; set.seed(41000000L + j$cc$id * 10000L + r)
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); C <- j$cc$C
    d <- switch(j$cc$fam, null = 0, intercept = C, slope = (C - 1) * eta0, link = C * (eta0^2 - 1),
                ushape = C * (X[, 1]^2 - 1), thresh = C * pmax(X[, 2] - 1, 0), inter = C * X[, 1] * X[, 3])
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
    ev <- stored$events[stored$rep == r]
    if (length(ev) != 1L || ev != sum(y)) stop("data set mismatch ", j$f)
    set.seed(60000000L + j$cc$id * 10000L + r)
    data.frame(fam = j$cc$fam, C = C, n = n, rep = r, davis = davis_choice(y, eta0, X), ctp_flex = ctp_flex_choice(y, eta0, X))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}

jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 10))) {
  nm <- sprintf("%s_p%02d_e%03d_%s_r%04d.csv", CELLS$design[q], round(100 * CELLS$prev[q]), CELLS$events[q], CELLS$dep[q], ch[1])
  f <- file.path(OUT1, nm)
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(kind = "main", cc = CELLS[q, ], reps = if (SMOKE) ch[1] else ch, f = f, stored = file.path(HERE, "sim3_rows", nm))
}
for (q in seq_len(nrow(FC))) for (ch in split(seq_len(FC$reps[q]), ceiling(seq_len(FC$reps[q]) / 20))) {
  f <- file.path(OUT2, sprintf("%s_C%.2f_n%04d_r%04d.csv", FC$fam[q], FC$C[q], FC$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(kind = "focus", cc = FC[q, ], reps = if (SMOKE) ch[1] else ch, f = f)
}
if (SMOKE) {
  j <- jobs[[which(vapply(jobs, function(x) x$kind == "main" && x$cc$prev == 0.05 && x$cc$events == 200L && x$cc$design == "mixed10", TRUE))[1]]]
  j$f <- tempfile(fileext = ".csv"); t0 <- Sys.time(); main_chunk(j); print(utils::read.csv(j$f)); print(Sys.time() - t0)
  j <- jobs[[which(vapply(jobs, function(x) x$kind == "focus", TRUE))[1]]]
  j$f <- tempfile(fileext = ".csv"); main <- NULL; focus_chunk(j); print(utils::read.csv(j$f))
  quit("no")
}
cat(sprintf("davis | %d chunks | %s\n", length(jobs), format(Sys.time())))
if (length(jobs)) {
  jobs <- jobs[order(-vapply(jobs, function(j) if (j$kind == "main") ceiling(j$cc$events / j$cc$prev) else 1000, 0))]
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "20")))
  parallel::clusterExport(cl, c("HERE", "main_chunk", "focus_chunk", "A0", "BETA", "gen_X", "true_logit", "davis_choice", "ctp_flex_choice", "fitb", "METHODS"))
  r <- parallel::parLapplyLB(cl, jobs, function(j) tryCatch(if (j$kind == "main") main_chunk(j) else focus_chunk(j),
                                                        error = function(e) paste("ERROR", j$f, conditionMessage(e))), chunk.size = 1L)
  parallel::stopCluster(cl)
  print(grep("ERROR", unlist(r), value = TRUE))
}
cat("done", format(Sys.time()), "\n")
