## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# support_closure_full.R -- re-runs the split-first SUPPORT analysis with the same seeds as support_split_first.R and keeps
# the FULL result of every call (all 15 intersection p-values of the closure), for Figure 4 of the paper. The single-part
# p-values are checked against support_split_first.csv before anything is saved.
suppressMessages(library(splines))
HERE <- file.path(.ROOT, "option_a")
source(file.path(HERE, "R", "localize_groups.R"))
s <- utils::read.csv(file.path(.ROOT, "data/support2.csv"),
                     na.strings = c("NaN", "NA", ""))
d <- stats::na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                               crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                               acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
VARS <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co", "male")
F0 <- y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male
F2 <- y ~ age + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) + ns(wblc, 3) + num.co + male
v <- d[d$acute == 0, ]
set.seed(20261011L); A <- sort(sample(nrow(v), floor(nrow(v) / 2))); B <- setdiff(seq_len(nrow(v)), A)
out <- list()
for (m in c("M0", "M2")) {
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d[d$acute == 1, ], family = binomial()))
  p <- as.numeric(predict(fit, newdata = v, type = "response")); e <- qlogis(p)
  set.seed(20261012L)
  out[[m]][["verdict on A"]] <- localize_external(v$y[A], p[A], as.matrix(v[A, VARS]), M = 999L)
  upd <- list(
    "intercept update" = function() { f <- glm(v$y[A] ~ offset(e[A]), family = binomial()); plogis(coef(f)[1] + e[B]) },
    "logistic recalibration" = function() { f <- glm(v$y[A] ~ e[A], family = binomial()); plogis(coef(f)[1] + coef(f)[2] * e[B]) },
    "flexible recalibration" = function() { f <- glm(v$y[A] ~ ns(e[A], df = 3), family = binomial())
                                            plogis(cbind(1, predict(ns(e[A], df = 3), e[B])) %*% coef(f)) },
    "model revision" = function() { f <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = v[A, ], family = binomial()))
                                    predict(f, newdata = v[B, ], type = "response") })
  for (rn in names(upd)) {
    pB <- as.numeric(upd[[rn]]())
    set.seed(20261013L)
    out[[m]][[rn]] <- localize_external(v$y[B], pB, as.matrix(v[B, VARS]), M = 999L)
  }
}
## fingerprint: every single-part p-value must equal the stored table
ref <- utils::read.csv(file.path(HERE, "support_split_first.csv"))
for (i in seq_len(nrow(ref))) {
  r <- out[[ref$model[i]]][[if (ref$stage[i] == "verdict on A") "verdict on A" else ref$rung[i]]]
  stopifnot(isTRUE(all.equal(as.numeric(r$intersection[c("INTERCEPT", "SLOPE", "LINK", "COV")]),
                             as.numeric(ref[i, c("INTERCEPT", "SLOPE", "LINK", "COV")]), tolerance = 1e-9)))
}
saveRDS(out, file.path(HERE, "support_closure_full.rds"))
cat("saved; fingerprints match\n"); print(out$M0[["verdict on A"]]$intersection)
