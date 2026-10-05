## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# support_split_first.R -- the update-and-retest protocol done "split first" on the transported SUPPORT models:
# the chronic patients are split at random into A and B; the verdict is reached on A alone (Monte Carlo reference),
# the model is updated on A at the rung that verdict names (highest part: INTERCEPT -> intercept update, SLOPE ->
# logistic recalibration, LINK -> flexible recalibration, COV -> re-estimation), and the updated model is tested on B,
# which neither choice saw. All rungs are also tested on B for comparison. Output support_split_first.csv.
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
RUNG <- c(INTERCEPT = "intercept update", SLOPE = "logistic recalibration", LINK = "flexible recalibration", COV = "model revision")
res <- list()
row <- function(model, stage, rung, n, r, extra = list()) {
  g <- c("INTERCEPT", "SLOPE", "LINK", "COV")
  x <- data.frame(model = model, stage = stage, rung = rung, n = n,
                  named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
                  as.data.frame(t(stats::setNames(as.numeric(r$intersection[g]), g))), check.names = FALSE)
  for (k in names(extra)) x[[k]] <- extra[[k]]
  res[[length(res) + 1]] <<- x; print(x, row.names = FALSE); flush.console()
}
for (m in c("M0", "M2")) {
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d[d$acute == 1, ], family = binomial()))
  p <- as.numeric(predict(fit, newdata = v, type = "response")); e <- qlogis(p)
  set.seed(20261012L)
  vA <- localize_external(v$y[A], p[A], as.matrix(v[A, VARS]), M = 999L)
  top <- if (length(vA$named)) tail(intersect(names(RUNG), vA$named), 1) else NA
  row(m, "verdict on A", "--", length(A), vA, list(observed = mean(v$y[A]), predicted = mean(p[A])))
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
    row(m, if (!is.na(top) && rn == RUNG[[top]]) "retest on B (prescribed rung)" else "retest on B", rn, length(B),
        localize_external(v$y[B], pB, as.matrix(v[B, VARS]), M = 999L), list(observed = mean(v$y[B]), predicted = mean(pB)))
  }
}
cols <- unique(unlist(lapply(res, names)))
res <- lapply(res, function(r) { for (cc in setdiff(cols, names(r))) r[[cc]] <- NA; r[cols] })
utils::write.csv(do.call(rbind, res), file.path(HERE, "support_split_first.csv"), row.names = FALSE)
cat("done\n")
