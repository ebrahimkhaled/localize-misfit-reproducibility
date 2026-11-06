# check_size_cover.R -- coverage of the 95% noise interval of localize_size() under no misfit (true size of every part
# zero), 1,000 data sets at n = 500 and at n = 1,000, R = 200 chance draws. Output theory/CHECK_SIZE_COVER.csv.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
jobs <- expand.grid(n = c(500L, 1000L), block = 1:50)
one <- function(k) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_size.R"))
  n <- jobs$n[k]
  do.call(rbind, lapply(1:20, function(i) {
    r <- (jobs$block[k] - 1L) * 20L + i; set.seed(57000000L + 10000L * (n == 1000L) + r)
    X <- matrix(stats::rnorm(n * 5), n, 5); e <- as.numeric(cbind(1, X) %*% b0)
    y <- stats::rbinom(n, 1L, stats::plogis(e))
    s <- localize_size(y, stats::plogis(e), X, R = 200L, seed = 58000000L + r)
    data.frame(n = n, rep = r, s$parts)
  }))
}
cl <- parallel::makeCluster(20L); parallel::clusterExport(cl, c("HERE", "jobs", "one", "b0"))
d <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(jobs)), one, chunk.size = 1L)); parallel::stopCluster(cl)
res <- do.call(rbind, lapply(split(d, list(d$n, d$part)), function(g)
  data.frame(n = g$n[1], part = g$part[1], reps = nrow(g), cover = mean(g$lo <= 0 & g$hi >= 0), excess = mean(g$excess), se = sd(g$excess) / sqrt(nrow(g)))))
utils::write.csv(res, file.path(HERE, "theory", "CHECK_SIZE_COVER.csv"), row.names = FALSE); print(res, row.names = FALSE, digits = 3)
