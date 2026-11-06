# check_update.R -- the updated-model reference (R/localize_update.R) against the Monte Carlo reference (localize_external,
# Holm) on the SAME data sets, in the design of run_sim_groups.R, external validation, n = 1000, 200 data sets per row.
# The rows include the two cases that the Monte Carlo reference gets wrong: a fixed-size slope error and a gross,
# transport-type miscalibration (logit pi = -2 + 0.8 eta0: an intercept and a slope error only, repaired exactly by
# a logistic recalibration). Output theory/update_rows/<fam>_r<first>.csv (written as it goes) and theory/CHECK_UPDATE.csv.
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
OUT <- file.path(HERE, "theory", "update_rows"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
DEP <- list(null = function(X, e) 0, intercept = function(X, e) 0.4, slope = function(X, e) (0.6 - 1) * e,
            transport = function(X, e) -2 + (0.8 - 1) * e - (-0.5) * 0,  # logit pi = -2 + 0.8 eta0 (see header)
            link = function(X, e) 0.25 * (e^2 - 1), ushape = function(X, e) 0.8 * (X[, 1]^2 - 1),
            thresh = function(X, e) 4 * pmax(X[, 2] - 1, 0), inter = function(X, e) X[, 1] * X[, 3])
NEED <- c(null = "none", intercept = "INTERCEPT", slope = "SLOPE", transport = "SLOPE", link = "LINK", ushape = "COV",
          thresh = "COV", inter = "COV")
REPS <- if (SMOKE) 2L else 200L
jobs <- list()
for (f in names(DEP)) for (ch in split(seq_len(REPS), ceiling(seq_len(REPS) / 10))) {
  out <- file.path(OUT, sprintf("%s_r%04d.csv", f, ch[1]))
  if (!file.exists(out)) jobs[[length(jobs) + 1]] <- list(fam = f, reps = ch, f = out)
}
one <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
  R <- do.call(rbind, lapply(j$reps, function(r) {
    set.seed(53000000L + 10000L * match(j$fam, names(DEP)) + r)
    X <- matrix(stats::rnorm(5000), 1000, 5); e <- as.numeric(cbind(1, X) %*% b0)
    lp <- if (j$fam == "transport") -2 + 0.8 * e else e + DEP[[j$fam]](X, e)
    y <- stats::rbinom(1000, 1L, stats::plogis(lp))
    mc <- localize_external(y, stats::plogis(e), X, M = 999L, naming = "holm")
    up <- localize_external_update(y, stats::plogis(e), X, B = 199L)
    data.frame(fam = j$fam, rep = r, events = sum(y),
               t(stats::setNames(G %in% mc$named, paste0("mc_", G))), t(stats::setNames(G %in% up$named, paste0("up_", G))),
               t(stats::setNames(as.numeric(up$intersection[G]), paste0("p_up_", G))))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}
cat(sprintf("check_update%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "10")))
  parallel::clusterExport(cl, c("HERE", "one", "DEP", "b0"))
  invisible(parallel::parLapplyLB(cl, jobs, one, chunk.size = 1L)); parallel::stopCluster(cl)
}
d <- do.call(rbind, lapply(list.files(OUT, pattern = "\\.csv$", full.names = TRUE), utils::read.csv))
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
top <- function(M) apply(M, 1, function(v) if (any(v)) G[max(which(v))] else "none")
res <- do.call(rbind, lapply(split(d, d$fam), function(g) {
  mcM <- as.matrix(g[, paste0("mc_", G)]); upM <- as.matrix(g[, paste0("up_", G)])
  data.frame(fam = g$fam[1], reps = nrow(g), needs = NEED[[g$fam[1]]],
             t(stats::setNames(colMeans(mcM), paste0("mc_", G))), right_mc = mean(top(mcM) == NEED[[g$fam[1]]]),
             t(stats::setNames(colMeans(upM), paste0("up_", G))), right_up = mean(top(upM) == NEED[[g$fam[1]]]))
}))
res <- res[match(names(DEP), res$fam), ]
utils::write.csv(res, file.path(HERE, "theory", "CHECK_UPDATE.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
