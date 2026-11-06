# check_size_null.R -- under no misfit the observed intercept size and its chance level have the same law; compare
# their means over 400 data sets (R = 20 chance draws each, no bootstrap), n = 1000.
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_size.R"))
b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
res <- t(vapply(1:400, function(r) {
  set.seed(52000000L + r)
  X <- matrix(rnorm(5000), 1000, 5); e <- as.numeric(cbind(1, X) %*% b0); y <- rbinom(1000, 1L, plogis(e))
  s <- localize_size(y, plogis(e), X, R = 20L, B = 2L)
  c(obs = s$parts$size[1], chance = s$parts$chance[1], obs_cov = s$parts$size[4], chance_cov = s$parts$chance[4])
}, numeric(4)))
print(round(rbind(mean = colMeans(res), se = apply(res, 2, sd) / sqrt(nrow(res))), 3))
