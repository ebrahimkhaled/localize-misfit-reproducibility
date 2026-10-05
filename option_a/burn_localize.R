## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# burn_localize.R -- second application: mortality after burn injury, 1,000 patients from 40 burn facilities (the
# burn1000 data of Hosmer, Lemeshow and Sturdivant, R package aplore3). Model B0: death ~ age + tbsa + gender + race +
# inhalation injury + flame, linear on the logit scale (main effects).
#   (1) in-sample checking on all 1,000 patients, default and de-aliased link group (B = 199)
#   (2) cross-centre external validation: develop on a random half of the facilities, validate the frozen model on the
#       patients of the other half (M = 999); the reported split, and 20 further facility splits to show how stable
#       the verdict is
# Output burn_localize.csv (analyses) and burn_splits.csv (the 20 splits).
suppressMessages({ library(splines); library(aplore3) })
HERE <- file.path(.ROOT, "option_a")
source(file.path(HERE, "R", "localize_groups.R"))
b <- aplore3::burn1000
d <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"),
                facility = b$facility)
VARS <- c("age", "tbsa", "male", "white", "inh", "flame")
F0 <- y ~ age + tbsa + male + white + inh + flame
res <- list()
add <- function(analysis, n, r, extra = list()) {
  g <- intersect(c("INTERCEPT", "SLOPE", "LINK", "COV"), names(r$intersection))
  x <- data.frame(analysis = analysis, n = n, named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
                  as.data.frame(t(stats::setNames(as.numeric(r$intersection[g]), g))), check.names = FALSE)
  for (k in names(extra)) x[[k]] <- extra[[k]]
  res[[length(res) + 1]] <<- x; print(x, row.names = FALSE); flush.console()
}
cat("patients", nrow(d), "deaths", sum(d$y), "facilities", length(unique(d$facility)), "\n")
## (1) in-sample
fit <- glm(F0, data = d, family = binomial())
for (da in c(FALSE, TRUE)) {
  set.seed(20261020L)
  r <- localize_insample(fit, as.matrix(d[, VARS]), B = 199L, dealias = da)
  add(if (da) "in-sample, de-aliased" else "in-sample", nrow(d), r, list(dealiased = paste(r$dealiased, collapse = " ")))
}
## (2) cross-centre external validation
fac <- sort(unique(d$facility))
split_once <- function(seed) {
  set.seed(seed); dev_f <- sample(fac, length(fac) %/% 2)
  dv <- d[d$facility %in% dev_f, ]; v <- d[!d$facility %in% dev_f, ]
  f <- glm(F0, data = dv, family = binomial())
  list(dv = dv, v = v, p = as.numeric(predict(f, newdata = v, type = "response")))
}
s0 <- split_once(20261021L)
set.seed(20261022L)
add("external, cross-centre", nrow(s0$v), localize_external(s0$v$y, s0$p, as.matrix(s0$v[, VARS]), M = 999L),
    list(dev_n = nrow(s0$dv), dev_deaths = sum(s0$dv$y), deaths = sum(s0$v$y), observed = mean(s0$v$y), predicted = mean(s0$p),
         slope = unname(coef(glm(s0$v$y ~ qlogis(s0$p), family = binomial()))[2])))
spl <- do.call(rbind, lapply(1:20, function(k) {
  s <- split_once(20262000L + k); set.seed(20263000L + k)
  r <- localize_external(s$v$y, s$p, as.matrix(s$v[, VARS]), M = 999L)
  data.frame(split = k, n_val = nrow(s$v), deaths = sum(s$v$y), named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
             t(stats::setNames(c("INTERCEPT", "SLOPE", "LINK", "COV") %in% r$named, paste0("named_", c("INTERCEPT", "SLOPE", "LINK", "COV")))))
}))
print(spl)
cat("\nshare of 20 facility splits naming each group:\n"); print(colMeans(spl[, grep("^named_", names(spl))]))
cols <- unique(unlist(lapply(res, names)))
res <- lapply(res, function(r) { for (cc in setdiff(cols, names(r))) r[[cc]] <- NA; r[cols] })
utils::write.csv(do.call(rbind, res), file.path(HERE, "burn_localize.csv"), row.names = FALSE)
utils::write.csv(spl, file.path(HERE, "burn_splits.csv"), row.names = FALSE)
cat("done\n")
