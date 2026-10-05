## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# check_local.R -- numerical check of Theorems 3 and 4 (strong familywise error and consistent localization).
# The theorems are LOCAL: a departure of size s_n = C * sqrt(1000 / n) on the logit scale. The prediction:
#   * a group whose part of the departure is zero is named with probability -> at most alpha as n grows
#     (the wrong-group rate falls towards .05 even for the cells where it was .9-1.0 at fixed size);
#   * a group whose part is not zero keeps roughly constant power (local alternative).
# The contrast track keeps the size FIXED (s = C at every n): there the theorem promises nothing for the wrong
# groups (non-local leak), and the right group's power goes to 1 (consistency).
# Departures are pure for the design (weighted-orthogonal on each dataset, as in run_sim_groups.R):
#   slope_pure  (C) [eta0 off 1]            -> only SLOPE has a part        (external only)
#   link_pure   (C) [eta0^2 off {1, eta0}]  -> only LINK
#   cov_pure    (C) [x1 x3 off 1, eta0, ns(eta0, 5)] -> only COV
#   combo       link_pure + cov_pure        -> LINK and COV, never INTERCEPT or SLOPE
# Each dataset is written to disk as it finishes (one file per 20 replicates, resumable).
HERE <- file.path(.ROOT, "option_a")
OUT <- file.path(HERE, "theory", "rows_local"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

fam_C <- data.frame(fam = c("null", "slope_pure", "link_pure", "cov_pure", "combo"), C = c(0, -0.4, 0.5, 2, 0.5))
ext <- merge(merge(fam_C, data.frame(track = c("local", "fixed"))), data.frame(n = c(500L, 2000L, 8000L)))
ext <- ext[!(ext$fam == "null" & ext$track == "fixed"), ]
ext$setting <- "external"; ext$reps <- 400L
ins <- merge(fam_C[fam_C$fam %in% c("null", "link_pure", "cov_pure", "combo"), ], data.frame(n = c(500L, 2000L, 4000L)))
ins$track <- "local"; ins$setting <- "insample"; ins$reps <- 200L
## MODERATE departures (added after the first results): sizes at which the own group is still found almost always at
## n = 1,000, followed far enough (to n = 32,000) for the local regime to be reached; external validation only
mod <- merge(data.frame(fam = c("slope_mod", "link_mod", "cov_mod"), C = c(-0.25, 0.2, 0.7), track = "local"),
             data.frame(n = c(1000L, 4000L, 16000L, 32000L)))
mod$setting <- "external"; mod$reps <- ifelse(mod$n == 32000L, 200L, 400L)
CELLS <- rbind(ext[, c("setting", "fam", "C", "track", "n", "reps")], ins[, c("setting", "fam", "C", "track", "n", "reps")],
               mod[, c("setting", "fam", "C", "track", "n", "reps")])
if (SMOKE) { CELLS <- CELLS[CELLS$n == 500 | (CELLS$setting == "insample" & CELLS$n == 2000), ]; CELLS$reps <- 2L }
CELLS$id <- seq_len(nrow(CELLS))

jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 20))) {
  f <- file.path(OUT, sprintf("%s%s_%s_%s_n%05d_r%04d.csv", if (SMOKE) "SMOKE_" else "", CELLS$setting[q], CELLS$fam[q],
                              CELLS$track[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

one_chunk <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  ## seed from the cell's identity (not its row number), so adding cells never changes existing datasets
  key <- sum(utf8ToInt(paste(j$cc$setting, j$cc$fam, j$cc$track, j$cc$n)) * seq_len(nchar(paste(j$cc$setting, j$cc$fam, j$cc$track, j$cc$n))))
  R <- do.call(rbind, lapply(j$reps, function(r) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(63000000L + key * 1000L + r)
    n <- j$cc$n
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- stats::plogis(eta0); w0 <- p0 * (1 - p0)
    po <- function(h, Z) as.numeric(.wperp(matrix(h), Z, w0))
    s <- j$cc$C * if (j$cc$track == "local") sqrt(1000 / n) else 1
    lk <- po(eta0^2, cbind(1, eta0)); cv <- po(X[, 1] * X[, 3], cbind(1, eta0, splines::ns(eta0, df = 5)))
    d <- switch(j$cc$fam, null = 0, slope_pure = , slope_mod = s * po(eta0, matrix(1, n, 1)), link_pure = , link_mod = s * lk,
                cov_pure = , cov_mod = s * cv, combo = s * lk + 4 * s * cv)
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
    t0 <- proc.time()[["elapsed"]]
    out <- if (j$cc$setting == "external") localize_external(y, p0, X, M = 999L) else {
      dat <- data.frame(X, y = y)
      fit <- suppressWarnings(stats::glm(y ~ x1 + x2 + x3 + x4 + x5, data = dat, family = stats::binomial()))
      localize_insample(fit, X, B = 199L)
    }
    g <- c("INTERCEPT", "SLOPE", "LINK", "COV")
    data.frame(setting = j$cc$setting, fam = j$cc$fam, C = j$cc$C, track = j$cc$track, n = n, s = s, rep = r,
               t(stats::setNames(g %in% out$named, paste0("named_", g))), secs = proc.time()[["elapsed"]] - t0)
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("check-local%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("CL_WORKERS", "10")))
  parallel::clusterExport(cl, c("HERE", "one_chunk"), envir = environment())
  ## longest jobs first (in-sample, large n) so the tail is short
  ord <- order(-(sapply(jobs, function(j) j$cc$n) * ifelse(sapply(jobs, function(j) j$cc$setting) == "insample", 20, 1)))
  invisible(parallel::parLapplyLB(cl, jobs[ord], function(j) tryCatch(one_chunk(j), error = function(e) NA), chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
