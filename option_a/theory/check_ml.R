# check_ml.R -- a machine-learning frozen predictor: a neural network (nnet, 5 hidden units) trained once on 3,000 patients with a
# nonlinear true risk, then frozen. Validation outcomes (n = 1,000, 500 data sets per setting) are drawn with logit risk: the network's
# own score (no misfit); -1 + 0.8 score (transport, a logistic recalibration repairs it); score + 0.8 x1 x3 (covariate misfit).
# Output theory/CHECK_ML.csv.
suppressMessages({ library(splines); library(nnet) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
set.seed(64000000L)
Xd <- matrix(rnorm(3000 * 5), 3000, 5); yd <- rbinom(3000, 1L, plogis(-.5 + .5 * Xd[, 1] + .8 * Xd[, 2]^2 - .4 + .6 * Xd[, 3] + .4 * Xd[, 4] * Xd[, 5]))
nn <- nnet(Xd, yd, size = 5, decay = .01, maxit = 1000, entropy = TRUE, trace = FALSE)
jobs <- expand.grid(rep = 1:500, s = c("none", "transport", "cov"), stringsAsFactors = FALSE)
one <- function(k) {
  suppressMessages({ library(splines); library(nnet) })
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  set.seed(65000000L + k); n <- 1000L; X <- matrix(rnorm(n * 5), n, 5); p <- pmin(pmax(as.numeric(predict(nn, X)), 1e-4), 1 - 1e-4); e <- qlogis(p); s <- jobs$s[k]
  lp <- switch(s, none = e, transport = -1 + .8 * e, cov = e + .8 * X[, 1] * X[, 3])
  y <- rbinom(n, 1L, plogis(lp))
  r <- tryCatch(localize_external_update(y, p, X, B = 199L), error = function(x) NULL)
  if (is.null(r)) return(data.frame(s = s, ok = FALSE, INTERCEPT = NA, SLOPE = NA, LINK = NA, COV = NA))
  data.frame(s = s, ok = TRUE, t(setNames(G %in% r$named, G)))
}
cl <- parallel::makeCluster(21L); parallel::clusterExport(cl, c("HERE", "G", "nn", "jobs", "one"))
d <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(jobs)), one, chunk.size = 5L)); parallel::stopCluster(cl)
ok <- d[d$ok, ]; res <- do.call(rbind, lapply(split(ok, ok$s), function(g) data.frame(setting = g$s[1], reps = nrow(g), t(colMeans(g[, G])), any = mean(rowSums(g[, G]) > 0))))
utils::write.csv(res, file.path(HERE, "theory", "CHECK_ML.csv"), row.names = FALSE); print(res, row.names = FALSE, digits = 3); cat("done\n")
