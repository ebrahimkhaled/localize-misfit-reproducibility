## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# check_robust.R -- Theorem A.5: the multiplier calibration keeps strong familywise error control for departures of
# FIXED size, where the Monte Carlo calibration is guaranteed only locally (Proposition A.1(iii)).
# Departures are pure on the PROBABILITY scale and need no truncation: pi = p0 + s * w0 * g with g bounded on the
# sample and |s g| <= 0.95, so p0 - w0 < pi < p0 + w0 stays inside (0, 1). Purity is exact for the software's
# hypotheses on each data set (g is weighted-orthogonal to every member basis of the other groups):
#   link_b  g = a spline-in-eta direction of ns(eta, 5) orthogonal to {1, eta}: only LINK has a part
#           (COV bases are projected off ns(eta, 5), so they cannot read it)
#   cov_b   g = tanh(x1 x3) projected off {1, eta, ns(eta, 5)} and off every LINK basis: only COV has a part
# Both calibrations are run on the same data sets. Null cells check the multiplier's level. External validation,
# five-covariate design of the paper, M = 999.
HERE <- file.path(.ROOT, "option_a")
OUT <- file.path(HERE, "theory", "rows_robust"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
CELLS <- rbind(expand.grid(fam = c("link_b", "cov_b"), n = c(1000L, 4000L, 16000L), stringsAsFactors = FALSE),
               expand.grid(fam = c("null", "transport"), n = c(1000L, 4000L), stringsAsFactors = FALSE))
CELLS$reps <- ifelse(CELLS$n == 16000L, 200L, 400L)
if (SMOKE) { CELLS <- CELLS[CELLS$n == 1000L, ]; CELLS$reps <- 2L }
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 20))) {
  f <- file.path(OUT, sprintf("%s%s_n%05d_r%04d.csv", if (SMOKE) "SMOKE_" else "", CELLS$fam[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}
one_chunk <- function(j) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  key <- sum(utf8ToInt(paste("robust", j$cc$fam, j$cc$n)) * seq_len(nchar(paste("robust", j$cc$fam, j$cc$n))))
  R <- do.call(rbind, lapply(j$reps, function(r) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(71000000L + key * 1000L + r)
    n <- j$cc$n
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- stats::plogis(eta0); w0 <- p0 * (1 - p0)
    po <- function(h, Z) as.numeric(.wperp(matrix(h), Z, w0))
    g <- switch(j$cc$fam, null = , transport = rep(0, n),
                link_b = po(splines::ns(eta0, df = 5)[, 3], cbind(1, eta0)),
                cov_b = { gb <- .group_bases(eta0, X)
                          po(tanh(X[, 1] * X[, 3]), cbind(1, eta0, splines::ns(eta0, df = 5), do.call(cbind, gb$LINK))) })
    s <- if (j$cc$fam %in% c("null", "transport")) 0 else 0.95 / max(abs(g))
    ## transport: a gross miscalibration of the kind met when a model is transported (about 12% observed against 40%
    ## predicted), logit pi = -2 + 0.8 eta0; a function of the score only, so COV has no part and any COV named is false
    pt <- if (j$cc$fam == "transport") stats::plogis(-2 + 0.8 * eta0) else p0 + s * w0 * g
    y <- stats::rbinom(n, 1L, pt)
    G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
    do.call(rbind, lapply(c("montecarlo", "multiplier"), function(cal) {
      out <- localize_external(y, p0, X, M = 999L, calibration = cal)
      data.frame(fam = j$cc$fam, n = n, rep = r, calibration = cal, s = s,
                 t(stats::setNames(G %in% out$named, paste0("named_", G))))
    }))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}
cat(sprintf("check-robust%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("CL_WORKERS", "6")))
  parallel::clusterExport(cl, c("HERE", "one_chunk"))
  ord <- order(-sapply(jobs, function(j) j$cc$n))
  invisible(parallel::parLapplyLB(cl, jobs[ord], function(j) tryCatch(one_chunk(j), error = function(e) conditionMessage(e)), chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
