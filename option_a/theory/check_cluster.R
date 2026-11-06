# check_cluster.R -- the error rate of the updated-model reference when centres differ in baseline risk. A correct frozen
# model apart from a centre random effect u ~ N(0, tau^2): 40 centres, 25 patients each (the burn design), and 20 centres
# of 100; tau = 0, 0.3, 0.6. Plain reference against the cluster-aware reference with tau estimated by glmer on each
# data set. 400 data sets per row. Output theory/CHECK_CLUSTER.csv.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
b0 <- c(-1.7, 0.5, 0.8, 0.6, 0.4, 0.3)
CELLS <- expand.grid(tau = c(0, 0.3, 0.6), centres = c(40L, 20L)); CELLS$id <- seq_len(nrow(CELLS))
jobs <- do.call(c, lapply(seq_len(nrow(CELLS)), function(q) lapply(split(1:400, ceiling((1:400) / 20)), function(r) list(cc = CELLS[q, ], reps = r))))
one <- function(j) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  do.call(rbind, lapply(j$reps, function(r) {
    set.seed(56000000L + 10000L * j$cc$id + r); n <- 1000L; K <- j$cc$centres
    X <- matrix(stats::rnorm(n * 5), n, 5); e <- as.numeric(cbind(1, X) %*% b0)
    cl <- rep(seq_len(K), length.out = n); u <- stats::rnorm(K, 0, j$cc$tau)[cl]
    y <- stats::rbinom(n, 1L, stats::plogis(e + u)); p <- stats::plogis(e)
    pl <- localize_external_update(y, p, X, B = 199L)
    g <- suppressMessages(suppressWarnings(lme4::glmer(y ~ 1 + offset(e) + (1 | cl), family = stats::binomial())))
    th <- as.numeric(attr(lme4::VarCorr(g)$cl, "stddev"))
    ca <- localize_external_update(y, p, X, B = 199L, cluster = cl, tau = th)
    data.frame(tau = j$cc$tau, centres = K, rep = r, tau_hat = th, t(stats::setNames(G %in% pl$named, paste0("plain_", G))),
               t(stats::setNames(G %in% ca$named, paste0("clust_", G))))
  }))
}
cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "6")))
parallel::clusterExport(cl, c("HERE", "one", "b0", "G"))
d <- do.call(rbind, parallel::parLapplyLB(cl, jobs, one, chunk.size = 1L)); parallel::stopCluster(cl)
res <- do.call(rbind, lapply(split(d, list(d$tau, d$centres)), function(g) data.frame(
  tau = g$tau[1], centres = g$centres[1], reps = nrow(g), tau_hat = mean(g$tau_hat),
  t(colMeans(g[, paste0("plain_", G)])), plain_any = mean(rowSums(g[, paste0("plain_", G)]) > 0),
  t(colMeans(g[, paste0("clust_", G)])), clust_any = mean(rowSums(g[, paste0("clust_", G)]) > 0))))
utils::write.csv(res, file.path(HERE, "theory", "CHECK_CLUSTER.csv"), row.names = FALSE); print(res, row.names = FALSE, digits = 3)
