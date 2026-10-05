## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# run_purity_prob.R -- is the leak of large "pure" departures a scale effect? Departures built pure on the PROBABILITY
# scale: p_true = p0 + C * w0 * g, with g a weighted-orthogonal direction of one group (w0 = p0(1 - p0)), so the mean
# residual p_true - p0 = w0 * C * g lies in that group's part exactly, at any strength C. Theory: a group whose part
# is zero keeps its null law, so the wrong-group rate stays near alpha whatever C is. External validation, n = 1,000,
# frozen model eta0, Monte Carlo reference M = 999; 300 datasets per cell.
suppressMessages(library(splines))
HERE <- file.path(.ROOT, "option_a")
OUT <- file.path(HERE, "rows_prob"); dir.create(OUT, showWarnings = FALSE)
CELLS <- data.frame(fam = rep(c("slope_p", "link_p", "cov_p"), each = 3), C = c(.6, 1.2, 2.4, .5, 1, 2, .5, 1, 2))
CELLS$id <- seq_len(nrow(CELLS))
jobs <- list()
for (q in seq_len(nrow(CELLS))) for (ch in split(1:300, ceiling(1:300 / 20))) {
  f <- file.path(OUT, sprintf("%s_C%.2f_r%04d.csv", CELLS$fam[q], CELLS$C[q], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], reps = ch, f = f)
}
one_chunk <- function(j) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R"))
  b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  R <- do.call(rbind, lapply(j$reps, function(r) {
    n <- 1000L; set.seed(52000000L + j$cc$id * 10000L + r)
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- stats::plogis(eta0); w0 <- p0 * (1 - p0)
    perp <- function(h, Z) { g <- as.numeric(.wperp(matrix(h), Z, w0)); g / sqrt(stats::weighted.mean(g^2, w0)) }
    g <- switch(j$cc$fam,
                slope_p = perp(eta0, matrix(1, n, 1)),
                link_p  = perp(eta0^2, cbind(1, eta0)),
                cov_p   = perp(X[, 1] * X[, 3], cbind(1, eta0, splines::ns(eta0, df = 5))))
    pt <- pmin(pmax(p0 + j$cc$C * w0 * g, 0.001), 0.999)
    y <- stats::rbinom(n, 1L, pt)
    out <- localize_external(y, p0, X, M = 999L)
    G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
    data.frame(fam = j$cc$fam, C = j$cc$C, rep = r, clipped = mean(pt <= 0.001 | pt >= 0.999),
               t(stats::setNames(G %in% out$named, paste0("named_", G))))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f); j$f
}
cl <- parallel::makeCluster(12L); parallel::clusterExport(cl, c("HERE", "one_chunk"))
invisible(parallel::parLapplyLB(cl, jobs, one_chunk, chunk.size = 1L)); parallel::stopCluster(cl)
d <- do.call(rbind, lapply(list.files(OUT, "\\.csv$", full.names = TRUE), utils::read.csv))
own <- c(slope_p = "SLOPE", link_p = "LINK", cov_p = "COV")
d$own <- d[cbind(seq_len(nrow(d)), match(paste0("named_", own[d$fam]), names(d)))]
oth <- sapply(seq_len(nrow(d)), function(i) any(unlist(d[i, setdiff(paste0("named_", c("INTERCEPT", "SLOPE", "LINK", "COV")), paste0("named_", own[d$fam[i]]))])))
d$wrong <- oth
S <- aggregate(cbind(own, wrong, clipped) ~ fam + C, d, mean)
print(S, row.names = FALSE, digits = 3)
utils::write.csv(S, file.path(HERE, "PURITY_PROB.csv"), row.names = FALSE)
