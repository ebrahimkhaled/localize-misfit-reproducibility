# run_rivals_update.R -- the updated-model reference on the SAME data sets as the rival comparison (run_vergouwe.R and
# run_sim_groups.R: same seeds, every data set checked against the stored event count). External validation, the seven
# rows of the rivals table, n = 500 and 1,000. Output rows_update/<fam>_C<C>_n<n>_r<first>.csv, 20 data sets per file.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
OUT <- file.path(HERE, "rows_update"); dir.create(OUT, showWarnings = FALSE)

ext <- data.frame(fam = c("null", "intercept", "intercept", "slope", "slope", "slope_pure", "link", "link_pure", "link_pure",
                          "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .2, .4, .8, .6, .6, .25, .25, .5, .8, 4, 1, 1, 2))
ins <- data.frame(fam = c("null", "link", "link_pure", "link_pure", "ushape", "thresh", "inter", "cov_pure", "cov_pure"),
                  C   = c(0, .25, .25, .5, .8, 4, 1, 1, 2))
CELLS <- rbind(merge(cbind(setting = "external", ext), data.frame(n = c(500L, 1000L))),
               merge(cbind(setting = "insample", ins), data.frame(n = c(500L, 1000L))))
CELLS$reps <- ifelse(CELLS$fam == "null", 1000L, 300L)
CELLS$id <- seq_len(nrow(CELLS))
keep <- paste(c("null", "intercept", "slope", "link", "ushape", "thresh", "inter"), c(0, .4, .6, .25, .8, 4, 1))
CELLS <- CELLS[CELLS$setting == "external" & paste(CELLS$fam, CELLS$C) %in% keep, ]
stopifnot(nrow(CELLS) == 14L)
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[q]), ceiling(seq_len(CELLS$reps[q]) / 20))) {
  f <- file.path(OUT, sprintf("%s_C%.2f_n%04d_r%04d.csv", CELLS$fam[q], CELLS$C[q], CELLS$n[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}

one_chunk <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3); G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
  stored <- utils::read.csv(file.path(HERE, "rows", sprintf("external_%s_C%.2f_n%04d_r%04d.csv", j$cc$fam, j$cc$C, j$cc$n, j$reps[1])))
  R <- do.call(rbind, lapply(j$reps, function(r) {
    n <- j$cc$n; set.seed(41000000L + j$cc$id * 10000L + r)
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); C <- j$cc$C
    d <- switch(j$cc$fam, null = 0, intercept = C, slope = (C - 1) * eta0, link = C * (eta0^2 - 1),
                ushape = C * (X[, 1]^2 - 1), thresh = C * pmax(X[, 2] - 1, 0), inter = C * X[, 1] * X[, 3])
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
    ev <- stored$events[stored$rep == r]
    if (length(ev) != 1L || ev != sum(y)) stop(sprintf("data set mismatch: %s C=%.2f n=%d rep %d", j$cc$fam, C, n, r))
    set.seed(47000000L + j$cc$id * 10000L + r)
    up <- localize_external_update(y, stats::plogis(eta0), X, B = 199L)
    data.frame(fam = j$cc$fam, C = C, n = n, rep = r, t(stats::setNames(G %in% up$named, paste0("named_", G))),
               t(stats::setNames(as.numeric(up$intersection[G]), paste0("p_", G))))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}
cat(sprintf("rivals update | %d chunks | %s\n", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "10")))
  parallel::clusterExport(cl, c("HERE", "one_chunk"))
  r <- parallel::parLapplyLB(cl, jobs, function(j) tryCatch(one_chunk(j), error = function(e) paste("ERROR", j$f, conditionMessage(e))), chunk.size = 1L)
  parallel::stopCluster(cl)
  print(grep("ERROR", unlist(r), value = TRUE))
}
cat("done", format(Sys.time()), "\n")
