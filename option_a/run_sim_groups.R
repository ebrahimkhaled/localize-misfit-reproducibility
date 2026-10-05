## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# run_sim_groups.R -- operating characteristics of the orthogonal-group localization (localize_groups.R).
# Five N(0,1) covariates, eta0 = -0.5 + 0.5 x1 + 0.8 x2 + 0.6 x3 + 0.4 x4 + 0.3 x5; true logit = eta0 + d.
# Departures d (C = strength):
#   null       0
#   intercept  C                                   (external only: absorbed in-sample)
#   slope      (C - 1) eta0                        (external only; also shifts the mean risk)
#   slope_pure (C - 1) [eta0 projected off 1]      (external only; slope alone)
#   link       C (eta0^2 - 1)                      natural curvature (has a linear part too)
#   link_pure  C * [eta0^2 projected off {1, eta0}]           only curvature of the score
#   ushape     C (x1^2 - 1)
#   thresh     C max(x2 - 1, 0)
#   inter      C x1 x3
#   cov_pure   C * [x1 x3 projected off functions of eta0 (1, eta0, ns(eta0, 5))]   only within-pattern misfit
# The projections use the true weights p0(1 - p0) on each dataset, so "pure" means pure for that design.
# EXTERNAL validation: the frozen model is eta0 (nothing refitted), Monte Carlo reference M = 999.
# IN-SAMPLE checking: the analyst fits logit ~ x1 + ... + x5; parametric bootstrap B = 199 with refits.
# Output rows/<setting>_<fam>_C<C>_n<n>_r<first>.csv (20 replicates per file; resumable).
HERE <- file.path(.ROOT, "option_a")
OUT <- file.path(HERE, "rows"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

ext <- data.frame(fam = c("null", "intercept", "intercept", "slope", "slope", "slope_pure", "link", "link_pure", "link_pure",
                          "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .2, .4, .8, .6, .6, .25, .25, .5, .8, 4, 1, 1, 2))
ins <- data.frame(fam = c("null", "link", "link_pure", "link_pure", "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .25, .25, .5, .8, 4, 1, 1, 2))
CELLS <- rbind(merge(cbind(setting = "external", ext), data.frame(n = c(500L, 1000L))),
               merge(cbind(setting = "insample", ins), data.frame(n = c(500L, 1000L))))
CELLS$reps <- ifelse(CELLS$fam == "null", 1000L, 300L)
if (SMOKE) CELLS$reps <- 2L
CELLS$id <- seq_len(nrow(CELLS))
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 20))) {
  f <- file.path(OUT, sprintf("%s_%s_C%.2f_n%04d_r%04d.csv", CELLS$setting[q], CELLS$fam[q], CELLS$C[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

one_chunk <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
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
    t0 <- proc.time()[["elapsed"]]
    out <- if (j$cc$setting == "external") localize_external(y, p0, X, M = 999L) else {
      dat <- data.frame(X, y = y)
      fit <- suppressWarnings(stats::glm(y ~ x1 + x2 + x3 + x4 + x5, data = dat, family = stats::binomial()))
      localize_insample(fit, X, B = 199L)
    }
    g <- c("INTERCEPT", "SLOPE", "LINK", "COV")
    top <- out$intersection[[length(out$intersection)]]
    data.frame(setting = j$cc$setting, fam = j$cc$fam, C = C, n = n, rep = r, events = sum(y), global_p = top,
               t(stats::setNames(g %in% out$named, paste0("named_", g))),
               t(stats::setNames(out$intersection[intersect(g, names(out$intersection))], paste0("p_", intersect(g, names(out$intersection))))),
               secs = proc.time()[["elapsed"]] - t0, check.names = FALSE)
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("sim-groups%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", if (SMOKE) "6" else "12")))
  parallel::clusterExport(cl, c("HERE", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk, chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
