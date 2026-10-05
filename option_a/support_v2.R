## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# support_v2.R -- the SUPPORT application extended with the results added to Web Appendix A on 2026-10-05:
#   (a) in-sample checking with the de-aliased link group (SUPPORT has a binary covariate and skewed laboratory values,
#       so E(X | score) is not linear and in-sample LINK can carry covariate structure)
#   (b) external validation with the multiplier calibration (error control for departures of any size)
#   (c) planning: the smallest miscalibration in the large and slope error each external analysis can detect with
#       power .80, from the covariates and frozen predictions alone (Corollary "Power and sample size")
#   (d) updating: the transported models are updated on a random half of the chronic patients at the rung the verdict
#       names (and one rung lower), and the updated models are validated on the other half (Proposition "Updating")
#   (e) any prediction model: a neural network (nnet, 5 hidden units) developed on the acute patients and validated on
#       the chronic ones, and on the random split
# Same data, models and splits as support_localize.R. Output support_v2.csv (one row per analysis).
suppressMessages({ library(splines); library(nnet) })
HERE <- file.path(.ROOT, "option_a")
source(file.path(HERE, "R", "localize_groups.R"))
s <- utils::read.csv(file.path(.ROOT, "data/support2.csv"),
                     na.strings = c("NaN", "NA", ""))
d <- data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma")))
d <- stats::na.omit(d)
VARS <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co", "male")
F0 <- y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male
F2 <- y ~ age + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) + ns(wblc, 3) + num.co + male
res <- list()
add <- function(analysis, model, n, r, extra = list()) {
  g <- intersect(c("INTERCEPT", "SLOPE", "LINK", "COV"), names(r$intersection))
  row <- data.frame(analysis = analysis, model = model, n = n,
                    named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
                    as.data.frame(t(stats::setNames(as.numeric(r$intersection[g]), g))), check.names = FALSE)
  for (k in names(extra)) row[[k]] <- extra[[k]]
  res[[length(res) + 1]] <<- row
  print(row, row.names = FALSE); flush.console()
}
plan <- function(p) {   # Corollary: detectable miscalibration in the large and slope error at power .80, alpha .05
  w <- p * (1 - p); eta <- qlogis(p); eb <- sum(w * eta) / sum(w); n <- length(p)
  list(detect_delta = sqrt(7.85 * mean(w) / n), detect_slope = sqrt(7.85 / (n * mean(w * (eta - eb)^2))))
}

## (a) in-sample, default and de-aliased
for (m in c("M0", "M2")) for (da in c(FALSE, TRUE)) {
  set.seed(20261005L)
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d, family = binomial()))
  add(if (da) "in-sample, de-aliased" else "in-sample", m, nrow(d), localize_insample(fit, as.matrix(d[, VARS]), B = 199L, dealias = da))
}

## (b, c) external validation: both calibrations, with the planning numbers
set.seed(20261006L)
dev_r <- sort(sample(nrow(d), floor(nrow(d) / 2)))
splits <- list(random = list(dev = dev_r, val = setdiff(seq_len(nrow(d)), dev_r)),
               transported = list(dev = which(d$acute == 1), val = which(d$acute == 0)))
frozen <- list()
for (sp in names(splits)) for (m in c("M0", "M2")) {
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d[splits[[sp]]$dev, ], family = binomial()))
  v <- d[splits[[sp]]$val, ]
  p <- as.numeric(predict(fit, newdata = v, type = "response"))
  frozen[[paste(sp, m)]] <- list(v = v, p = p)
  for (cal in c("montecarlo", "multiplier")) {
    set.seed(20261007L)
    add(paste0("external ", sp, ", ", cal), m, nrow(v), localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L, calibration = cal),
        c(list(observed = mean(v$y), predicted = mean(p)), plan(p)))
  }
}

## (d) updating the transported models on half of the chronic patients, validating on the other half
set.seed(20261008L)
for (m in c("M0", "M2")) {
  v <- frozen[[paste("transported", m)]]$v; p <- frozen[[paste("transported", m)]]$p
  upd <- sort(sample(nrow(v), floor(nrow(v) / 2))); tst <- setdiff(seq_len(nrow(v)), upd)
  e <- qlogis(p)
  rungs <- list(
    "intercept update"        = function(i) glm(v$y[i] ~ offset(e[i]), family = binomial()),
    "logistic recalibration"  = function(i) glm(v$y[i] ~ e[i], family = binomial()),
    "flexible recalibration"  = function(i) glm(v$y[i] ~ ns(e[i], df = 3), family = binomial()),
    "model revision"          = function(i) suppressWarnings(glm(if (m == "M0") F0 else F2, data = v[i, ], family = binomial())))
  for (rn in names(rungs)) {
    f <- rungs[[rn]](upd)
    ## predictions for the test half from the updating model, evaluated at the test half's scores (or covariates)
    pt <- switch(rn,
      "intercept update" = plogis(coef(f)[1] + e[tst]),
      "logistic recalibration" = plogis(coef(f)[1] + coef(f)[2] * e[tst]),
      "flexible recalibration" = { B <- ns(e[upd], df = 3); plogis(cbind(1, predict(B, e[tst])) %*% coef(f)) },
      "model revision" = predict(f, newdata = v[tst, ], type = "response"))
    for (cal in c("montecarlo", "multiplier")) {
      set.seed(20261009L)
      add(paste0("updated: ", rn, ", ", cal), m, length(tst),
          localize_external(v$y[tst], as.numeric(pt), as.matrix(v[tst, VARS]), M = 999L, calibration = cal),
          list(observed = mean(v$y[tst]), predicted = mean(pt)))
    }
  }
}

## (e) a neural network as the frozen model
for (sp in names(splits)) {
  dv <- d[splits[[sp]]$dev, ]; v <- d[splits[[sp]]$val, ]
  mu <- colMeans(dv[, VARS]); sdv <- apply(dv[, VARS], 2, sd)
  Z <- function(a) scale(as.matrix(a[, VARS]), mu, sdv)
  set.seed(20261010L)
  nn <- nnet::nnet(Z(dv), dv$y, size = 5, decay = 0.01, maxit = 1000, entropy = TRUE, trace = FALSE)
  p <- as.numeric(predict(nn, Z(v)))
  for (cal in c("montecarlo", "multiplier")) {
    set.seed(20261007L)
    add(paste0("external ", sp, ", ", cal), "neural network", nrow(v),
        localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L, calibration = cal),
        c(list(observed = mean(v$y), predicted = mean(p)), plan(p)))
  }
}

cols <- unique(unlist(lapply(res, names)))
res <- lapply(res, function(r) { for (cc in setdiff(cols, names(r))) r[[cc]] <- NA; r[cols] })
utils::write.csv(do.call(rbind, res), file.path(HERE, "support_v2.csv"), row.names = FALSE)
cat("done\n")
