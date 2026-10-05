## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# check_dealias.R -- Propositions A.2(iii) and A.7 in a design where E(X | score) is NOT linear (in-sample checking).
# Covariates: x1 ~ N(0,1), x2 = Exp(1) + 0.5 x1^2 (skewed, curved in x1), x3 ~ Bernoulli(0.4); model logit = eta0 =
# -1 + 0.6 x1 + 0.5 x2 + 0.5 x3 (main effects, correctly specified under the null).
# Departures (logit scale, true logit = eta0 + C h):
#   null      none
#   cov       h = x1 x2 projected off {1, eta0, ns(eta0, 5)}: pure covariate structure in the external sense; in-sample
#             part of it is aliased into LINK because E(x | eta0) is curved here
#   link      h = eta0^2 projected off {1, eta0}: pure link shape
# Each data set is analysed by the default in-sample procedure and by the de-aliased one (LINK bases also projected off
# the spline fit of each model column on the score). Reported: rate LINK named, COV named, any named.
HERE <- file.path(.ROOT, "option_a")
OUT <- file.path(HERE, "theory", "rows_dealias"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
CELLS <- data.frame(fam = c("null", "cov", "link"), C = c(0, 0.6, 0.3), n = 1000L, reps = 200L)
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 10))) {
  f <- file.path(OUT, sprintf("%s%s_n%05d_r%04d.csv", if (SMOKE) "SMOKE_" else "", CELLS$fam[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}
one_chunk <- function(j) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R"))
  key <- sum(utf8ToInt(paste("dealias", j$cc$fam, j$cc$n)) * seq_len(nchar(paste("dealias", j$cc$fam, j$cc$n))))
  R <- do.call(rbind, lapply(j$reps, function(r) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(72000000L + key * 1000L + r)
    n <- j$cc$n
    x1 <- stats::rnorm(n); x2 <- stats::rexp(n) + 0.5 * x1^2; x3 <- stats::rbinom(n, 1, 0.4)
    X <- cbind(x1 = x1, x2 = x2, x3 = x3)
    eta0 <- -1 + 0.6 * x1 + 0.5 * x2 + 0.5 * x3; p0 <- stats::plogis(eta0); w0 <- p0 * (1 - p0)
    po <- function(h, Z) as.numeric(.wperp(matrix(h), Z, w0))
    h <- switch(j$cc$fam, null = 0, cov = po(x1 * x2, cbind(1, eta0, splines::ns(eta0, df = 5))),
                link = po(eta0^2, cbind(1, eta0)))
    if (j$cc$fam != "null") h <- h / sqrt(stats::weighted.mean(h^2, w0))
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + j$cc$C * h))
    dat <- data.frame(X, y = y)
    fit <- suppressWarnings(stats::glm(y ~ x1 + x2 + x3, data = dat, family = stats::binomial()))
    do.call(rbind, lapply(c(FALSE, TRUE), function(da) {
      out <- localize_insample(fit, X, B = 199L, dealias = da)
      data.frame(fam = j$cc$fam, n = n, rep = r, dealias = da,
                 named_LINK = "LINK" %in% out$named, named_COV = "COV" %in% out$named)
    }))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}
cat(sprintf("check-dealias%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("CL_WORKERS", "6")))
  parallel::clusterExport(cl, c("HERE", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, function(j) tryCatch(one_chunk(j), error = function(e) conditionMessage(e)), chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
