## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# run_groups_corruption.R -- the corrupted-records check for the orthogonal-group procedure of the Biometrics paper
# (option_a/R/localize_groups.R, in-sample, default and de-aliased). Same stored datasets as run_localize_corruption.R,
# rebuilt from their seeds and checked against the stored events and fitted slope: x ~ U(-3, 3), binary d, the CORRECT
# logistic model in x and d, and after the outcomes were drawn x multiplied by 4 or 8 in a few records. Any group
# named is a false one. Bootstrap B = 199 under Mersenne-Twister, seed + 7 (as in the stored study).
# Output rows_groups/<cell>_r<rep>.csv (one row per dataset and version; resumable).
BAT  <- file.path(.ROOT, "data/battery")
HERE <- file.path(.ROOT, "localize_robust")
SRC  <- file.path(.ROOT, "option_a/R/localize_groups.R")
OUT  <- file.path(HERE, "rows_groups"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
CELLS <- rbind(
  data.frame(block = "9b", file = c("clean", "x4", "x8"), n = 500L, k = c(0L, 3L, 3L), mult = c(1, 4, 8), reps = 100L),
  data.frame(block = "9", file = c("logit_clean_n1000", "logit_C1_r010_n1000"), n = 1000L, k = c(0L, 10L),
             mult = c(1, 4), reps = 100L))
jobs <- list()
for (q in seq_len(nrow(CELLS))) {
  S <- utils::read.csv(gzfile(file.path(BAT, CELLS$block[q], paste0(CELLS$file[q], "_pvalues.csv.gz"))))
  S <- S[order(S$rep), intersect(c("rep", "seed", "events", "b.x"), names(S))]
  for (i in seq_len(if (SMOKE) 1L else CELLS$reps[q])) {
    f <- file.path(OUT, sprintf("%s%s_%s_r%04d.csv", if (SMOKE) "SMOKE_" else "", CELLS$block[q], CELLS$file[q], S$rep[i]))
    if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], s = S[i, ], f = f)
  }
}
one_job <- function(j) {
  suppressMessages(library(splines)); source(SRC)
  n <- j$cc$n
  RNGkind("L'Ecuyer-CMRG"); set.seed(j$s$seed)
  x <- stats::runif(n, -3, 3); d <- stats::rbinom(n, 1, 0.5)
  y <- stats::rbinom(n, 1, stats::plogis(0.6 * x + 0.5 * d))
  if (j$cc$mult != 1) { idx <- sample.int(n, j$cc$k); x[idx] <- j$cc$mult * x[idx] }
  dat <- data.frame(y = y, x = x, d = d)
  fit <- suppressWarnings(stats::glm(y ~ x + d, data = dat, family = stats::binomial()))
  stopifnot(sum(y) == j$s$events, is.null(j$s$b.x) || is.na(j$s$b.x) || abs(stats::coef(fit)[["x"]] - j$s$b.x) < 1e-8)
  out <- do.call(rbind, lapply(c(FALSE, TRUE), function(da) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(j$s$seed + 7L)
    r <- localize_insample(fit, cbind(x = x, d = d), B = 199L, dealias = da)
    data.frame(block = j$cc$block, file = j$cc$file, rep = j$s$rep, dealias = da,
               named_LINK = "LINK" %in% r$named, named_COV = "COV" %in% r$named,
               p_LINK = unname(r$intersection["LINK"]), p_COV = unname(r$intersection["COV"]))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(out, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}
cat(sprintf("groups-corruption%s | %d jobs | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("LC_WORKERS", "4")))
  parallel::clusterExport(cl, c("SRC", "one_job"))
  invisible(parallel::parLapplyLB(cl, jobs, function(j) tryCatch(one_job(j), error = function(e) conditionMessage(e)), chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
