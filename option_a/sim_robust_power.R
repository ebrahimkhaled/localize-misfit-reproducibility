## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# sim_robust_power.R -- the ROBUST option (normal scores) on the SAME data sets as the main simulation (run_sim_groups.R; same cell
# table, same seeds, same data generation), n = 1,000. The procedure's own results are in rows/ already; this script
# adds, per data set, the p-values of
#   external validation: Hosmer-Lemeshow with frozen predictions (deciles, chi-square 10 df), Spiegelhalter's z test,
#                        the GiViTI calibration belt (external);
#   in-sample checking:  Hosmer-Lemeshow (10 groups, 8 df) and le Cessie-van Houwelingen (ebrahim.gof), the GiViTI
#                        calibration belt (internal).
# Output rows_compare/<setting>_<fam>_C<C>_n<n>_r<first>.csv (resumable).
suppressMessages(library(splines))
HERE <- file.path(.ROOT, "option_a")
OUT <- file.path(HERE, "rows_robustpower"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

## the cell table of run_sim_groups.R, rebuilt identically so that the cell ids (and so the seeds) match
ext <- data.frame(fam = c("null", "intercept", "intercept", "slope", "slope", "slope_pure", "link", "link_pure", "link_pure",
                          "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .2, .4, .8, .6, .6, .25, .25, .5, .8, 4, 1, 1, 2))
ins <- data.frame(fam = c("null", "link", "link_pure", "link_pure", "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .25, .25, .5, .8, 4, 1, 1, 2))
CELLS <- rbind(merge(cbind(setting = "external", ext), data.frame(n = c(500L, 1000L))),
               merge(cbind(setting = "insample", ins), data.frame(n = c(500L, 1000L))))
CELLS$reps <- ifelse(CELLS$fam == "null", 1000L, 300L)
CELLS$id <- seq_len(nrow(CELLS))
KEEP <- data.frame(fam = c("null", "intercept", "slope", "link", "ushape", "thresh", "inter"), C = c(0, .4, .6, .25, .8, 4, 1))
CELLS <- CELLS[CELLS$n == 1000L & paste(CELLS$fam, CELLS$C) %in% paste(KEEP$fam, KEEP$C), ]
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 20))) {
  f <- file.path(OUT, sprintf("%s%s_%s_C%.2f_n%04d_r%04d.csv", if (SMOKE) "SMOKE_" else "", CELLS$setting[q], CELLS$fam[q],
                              CELLS$C[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

one_chunk <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  hl_ext <- function(y, p) {   # frozen predictions: deciles of p, chi-square with 10 df
    g <- cut(p, unique(stats::quantile(p, 0:10 / 10)), include.lowest = TRUE)
    o <- tapply(y, g, sum); e <- tapply(p, g, sum); m <- tapply(p, g, length); pb <- e / m
    stats::pchisq(sum((o - e)^2 / (m * pb * (1 - pb))), nlevels(g), lower.tail = FALSE)
  }
  spieg <- function(y, p) { z <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * p * (1 - p))); 2 * stats::pnorm(-abs(z)) }
  belt <- function(y, p, devel) tryCatch(givitiR::givitiCalibrationBelt(y, p, devel = devel)$p.value, error = function(e) NA_real_)
  R <- do.call(rbind, lapply(j$reps, function(r) {
    n <- j$cc$n; set.seed(41000000L + j$cc$id * 10000L + r)        # exactly the seed of run_sim_groups.R
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
    G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
    out <- if (j$cc$setting == "external") localize_external(y, p0, X, M = 999L, robust = TRUE) else {
      dat <- data.frame(X, y = y)
      fit <- suppressWarnings(stats::glm(y ~ x1 + x2 + x3 + x4 + x5, data = dat, family = stats::binomial()))
      localize_insample(fit, X, B = 199L, robust = TRUE)
    }
    data.frame(setting = j$cc$setting, fam = j$cc$fam, C = C, n = n, rep = r, events = sum(y),
               global_p = out$intersection[[length(out$intersection)]], t(stats::setNames(G %in% out$named, paste0("named_", G))))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}

cat(sprintf("compare-detect%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("CD_WORKERS", "6")))
  parallel::clusterExport(cl, c("HERE", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, function(j) tryCatch(one_chunk(j), error = function(e) conditionMessage(e)), chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
