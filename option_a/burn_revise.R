## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# burn_revise.R -- the burn application, revision step. The main-effects model B0 is named COV in-sample
# (burn_localize.R). Revisions, each checked in-sample (B = 199):
#   B1: B0 with a natural spline in age (3 df)
#   B2: B1 plus the age x inhalation-injury interaction
#   B3: B2 with a natural spline in total burn surface area (3 df)
# B3 was found by looking at these data, so it is also validated across the same 20 random facility splits as B0
# (develop on half of the facilities, validate on the other half, M = 999). Output burn_revise.csv, burn_revise_splits.csv.
suppressMessages({ library(splines); library(aplore3) })
HERE <- file.path(.ROOT, "option_a")
source(file.path(HERE, "R", "localize_groups.R"))
b <- aplore3::burn1000
d <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"),
                facility = b$facility)
VARS <- c("age", "tbsa", "male", "white", "inh", "flame")
FM <- list(B0 = y ~ age + tbsa + male + white + inh + flame,
           B1 = y ~ ns(age, 3) + tbsa + male + white + inh + flame,
           B2 = y ~ ns(age, 3) + tbsa + male + white + inh + flame + age:inh,
           B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)
ins <- do.call(rbind, lapply(names(FM), function(m) {
  fit <- glm(FM[[m]], data = d, family = binomial()); set.seed(20261020L)
  r <- localize_insample(fit, as.matrix(d[, VARS]), B = 199L)
  data.frame(model = m, aic = AIC(fit), named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
             LINK = unname(r$intersection["LINK"]), COV = unname(r$intersection["COV"]))
}))
print(ins, row.names = FALSE)
utils::write.csv(ins, file.path(HERE, "burn_revise.csv"), row.names = FALSE)
## the same 20 facility splits as burn_localize.R (seeds 20262000 + k for the split, 20263000 + k for the reference)
fac <- sort(unique(d$facility)); G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
spl <- do.call(rbind, lapply(1:20, function(k) do.call(rbind, lapply(c("B0", "B3"), function(m) {
  set.seed(20262000L + k); dev_f <- sample(fac, length(fac) %/% 2)
  dv <- d[d$facility %in% dev_f, ]; v <- d[!d$facility %in% dev_f, ]
  f <- suppressWarnings(glm(FM[[m]], data = dv, family = binomial()))
  p <- as.numeric(predict(f, newdata = v, type = "response"))
  set.seed(20263000L + k); r <- localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L)
  data.frame(split = k, model = m, n_val = nrow(v), deaths = sum(v$y),
             t(stats::setNames(G %in% r$named, paste0("named_", G))), any = length(r$named) > 0)
}))))
utils::write.csv(spl, file.path(HERE, "burn_revise_splits.csv"), row.names = FALSE)
cat("\nshare of the 20 facility splits naming each group:\n")
print(aggregate(cbind(named_INTERCEPT, named_SLOPE, named_LINK, named_COV, any) ~ model, spl, mean))
