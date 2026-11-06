# check_size.R -- does localize_size() put the size where the misfit is? Simulation design of run_sim_groups.R,
# external validation, n = 1000; 40 data sets per departure (R = B = 50 for speed). Output theory/CHECK_SIZE.csv.
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_size.R"))
b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
deps <- list(null = function(X, e) 0, intercept = function(X, e) 0.4, slope = function(X, e) -0.4 * e,
             link = function(X, e) 0.25 * (e^2 - 1), ushape = function(X, e) 0.8 * (X[, 1]^2 - 1),
             inter = function(X, e) X[, 1] * X[, 3])
one <- function(k) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_size.R"))
  fam <- jobs$fam[k]; r <- jobs$rep[k]; set.seed(51000000L + 1000L * match(fam, names(deps)) + r)
  X <- matrix(rnorm(5000), 1000, 5); e <- as.numeric(cbind(1, X) %*% b0)
  y <- rbinom(1000, 1L, plogis(e + deps[[fam]](X, e)))
  s <- localize_size(y, plogis(e), X, R = 100L)
  data.frame(fam = fam, rep = r, part = s$parts$part, size = s$parts$size, chance = s$parts$chance,
             excess = s$parts$excess, lo = s$parts$lo, hi = s$parts$hi)
}
jobs <- expand.grid(fam = names(deps), rep = 1:40, stringsAsFactors = FALSE)
cl <- parallel::makeCluster(8L); parallel::clusterExport(cl, c("HERE", "jobs", "deps", "b0", "one"))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(jobs)), one)); parallel::stopCluster(cl)
utils::write.csv(res, file.path(HERE, "theory", "CHECK_SIZE.csv"), row.names = FALSE)
agg <- aggregate(cbind(size, chance, excess, covered0 = (lo <= 0)) ~ fam + part, res, mean)
agg$part <- factor(agg$part, c("INTERCEPT", "SLOPE", "LINK", "COV"))
print(agg[order(agg$fam, agg$part), ], row.names = FALSE, digits = 3)
