# run_vergouwe.R -- the rivals that choose an update for a frozen model, on the SAME external-validation data sets as
# run_sim_groups.R (same seeds; every data set is checked against the stored event count before it is used):
#   vergouwe  the closed testing procedure of Vergouwe et al. (2017, Stat Med 36:4529-4539), translated from their
#             Appendix R code to base glm: likelihood-ratio tests of the revised model (all p coefficients refitted)
#             against the original (p+1 df), recalibration in the large (p df) and recalibration (p-1 df), in that
#             order, each at alpha. Decision: original / intercept / recalibration / revision.
#   aic, bic  the update with the smallest AIC (BIC) among original, intercept, recalibration, flexible recalibration
#             (natural spline in the score, 3 df) and revision.
# Output vergouwe_rows/<fam>_C<C>_n<n>_r<first>.csv, 20 data sets per file, written as they finish (resumable).
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
OUT <- file.path(HERE, "vergouwe_rows"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

## the cells and their ids exactly as in run_sim_groups.R (the id enters the seed)
ext <- data.frame(fam = c("null", "intercept", "intercept", "slope", "slope", "slope_pure", "link", "link_pure", "link_pure",
                          "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .2, .4, .8, .6, .6, .25, .25, .5, .8, 4, 1, 1, 2))
ins <- data.frame(fam = c("null", "link", "link_pure", "link_pure", "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .25, .25, .5, .8, 4, 1, 1, 2))
CELLS <- rbind(merge(cbind(setting = "external", ext), data.frame(n = c(500L, 1000L))),
               merge(cbind(setting = "insample", ins), data.frame(n = c(500L, 1000L))))
CELLS$reps <- ifelse(CELLS$fam == "null", 1000L, 300L)
CELLS$id <- seq_len(nrow(CELLS))
CELLS <- CELLS[CELLS$setting == "external", ]
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 20))) {
  f <- file.path(OUT, sprintf("%s_C%.2f_n%04d_r%04d.csv", CELLS$fam[q], CELLS$C[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

## Vergouwe et al. (2017), Appendix: log-likelihoods of the four models, then the three ordered LR tests
vergouwe_ctp <- function(coefs, X, y, alpha = 0.05) {
  lp_old <- as.numeric(X %*% coefs[-1])
  ll <- function(lp) sum(y * lp - log1p(exp(lp)))
  int <- stats::coef(stats::glm(y ~ 1, offset = lp_old, family = stats::binomial()))
  rc <- stats::coef(stats::glm(y ~ lp_old, family = stats::binomial()))
  rv <- stats::coef(suppressWarnings(stats::glm(y ~ X, family = stats::binomial())))
  l0 <- ll(coefs[1] + lp_old); l1 <- ll(int + lp_old); l2 <- ll(rc[1] + rc[2] * lp_old); l3 <- ll(cbind(1, X) %*% rv)
  p <- ncol(X)
  pB <- stats::pchisq(2 * (l3 - l0), p + 1, lower.tail = FALSE)
  pC <- stats::pchisq(2 * (l3 - l1), p, lower.tail = FALSE)
  pD <- stats::pchisq(2 * (l3 - l2), p - 1, lower.tail = FALSE)
  choice <- if (pB >= alpha) "original" else if (pC >= alpha) "intercept" else if (pD >= alpha) "recalibration" else "revision"
  list(choice = choice, pB = pB, pC = pC, pD = pD, ll = c(l0, l1, l2, l3))
}

one_chunk <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  stored <- utils::read.csv(file.path(HERE, "rows", sprintf("external_%s_C%.2f_n%04d_r%04d.csv", j$cc$fam, j$cc$C, j$cc$n, j$reps[1])))
  R <- do.call(rbind, lapply(j$reps, function(r) {
    n <- j$cc$n; set.seed(41000000L + j$cc$id * 10000L + r)
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- stats::plogis(eta0); w0 <- p0 * (1 - p0)
    proj_off <- function(h, Z) as.numeric(.wperp(matrix(h), Z, w0))
    C <- j$cc$C
    d <- switch(j$cc$fam, null = 0, intercept = C, slope = (C - 1) * eta0,
                slope_pure = (C - 1) * proj_off(eta0, matrix(1, n, 1)), link = C * (eta0^2 - 1),
                link_pure = C * proj_off(eta0^2, cbind(1, eta0)),
                ushape = C * (X[, 1]^2 - 1), thresh = C * pmax(X[, 2] - 1, 0), inter = C * X[, 1] * X[, 3],
                cov_pure = C * proj_off(X[, 1] * X[, 3], cbind(1, eta0, splines::ns(eta0, df = 5))))
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
    ## fingerprint: the same data set as the localization study
    ev <- stored$events[stored$rep == r]
    if (length(ev) != 1L || ev != sum(y)) stop(sprintf("data set mismatch: %s C=%.2f n=%d rep %d", j$cc$fam, C, n, r))
    v <- vergouwe_ctp(b0, X, y)
    ## AIC / BIC over the update ladder (k = number of re-estimated parameters)
    flex <- stats::glm(y ~ splines::ns(eta0, df = 3), family = stats::binomial())
    ll <- c(original = v$ll[1], intercept = v$ll[2], recalibration = v$ll[3],
            flexible = as.numeric(stats::logLik(flex)), revision = v$ll[4])
    k <- c(original = 0, intercept = 1, recalibration = 2, flexible = 4, revision = 6)
    data.frame(fam = j$cc$fam, C = C, n = n, rep = r, events = sum(y), vergouwe = v$choice,
               pB = v$pB, pC = v$pC, pD = v$pD,
               aic = names(which.min(-2 * ll + 2 * k)), bic = names(which.min(-2 * ll + log(n) * k)))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("vergouwe%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "8")))
  parallel::clusterExport(cl, c("HERE", "one_chunk", "vergouwe_ctp"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk, chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
