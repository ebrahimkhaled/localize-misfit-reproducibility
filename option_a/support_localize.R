## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# support_localize.R -- the orthogonal-group localization on the SUPPORT in-hospital mortality data (public SUPPORT2
# release, UCI). Models: M0, linear in age, mean arterial pressure, heart rate, respiratory rate, temperature,
# creatinine, sodium, white cell count, number of comorbidities and sex; M2, natural splines (3 df) in the seven
# physiological variables (the ladder model of the DeepGOF-1 application).
#   (1) in-sample checking of M0 and M2 on all complete cases (bootstrap B = 199)
#   (2) external validation, random split: develop on a random half, validate the frozen model on the other half
#   (3) external validation, transported: develop on the acute disease classes (ARF/MOSF, coma), validate on the
#       chronic ones (COPD/CHF/cirrhosis, cancer)
# External validation uses the Monte Carlo reference (M = 999). Output support_localize.csv
suppressMessages(library(splines))
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
cat("complete cases:", nrow(d), " deaths:", sum(d$y), " acute:", sum(d$acute), "\n")
res <- list()
add <- function(setting, model, n, r) {
  res[[length(res) + 1]] <<- data.frame(setting = setting, model = model, n = n,
    named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
    t(r$intersection[intersect(c("INTERCEPT", "SLOPE", "LINK", "COV"), names(r$intersection))]),
    global = r$intersection[[length(r$intersection)]], check.names = FALSE)
  print(utils::tail(res, 1)[[1]], row.names = FALSE); flush.console()
}
## (1) in-sample
for (m in c("M0", "M2")) {
  set.seed(20261005L)
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d, family = binomial()))
  add("in-sample", m, nrow(d), localize_insample(fit, as.matrix(d[, VARS]), B = 199L))
}
## (2) random split and (3) transported
set.seed(20261006L)
dev_r <- sort(sample(nrow(d), floor(nrow(d) / 2)))
splits <- list(random = list(dev = dev_r, val = setdiff(seq_len(nrow(d)), dev_r)),
               transported = list(dev = which(d$acute == 1), val = which(d$acute == 0)))
for (sp in names(splits)) for (m in c("M0", "M2")) {
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d[splits[[sp]]$dev, ], family = binomial()))
  v <- d[splits[[sp]]$val, ]
  p <- as.numeric(predict(fit, newdata = v, type = "response"))
  set.seed(20261007L)
  add(paste("external,", sp), m, nrow(v), localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L))
  cat(sprintf("   observed %.3f vs predicted %.3f; calibration slope %.2f\n", mean(v$y), mean(p),
              coef(glm(v$y ~ qlogis(p), family = binomial()))[2]))
}
## in-sample rows have no INTERCEPT/SLOPE columns: fill them with NA before stacking
cols <- unique(unlist(lapply(res, names)))
res <- lapply(res, function(r) { for (cc in setdiff(cols, names(r))) r[[cc]] <- NA; r[cols] })
utils::write.csv(do.call(rbind, res), file.path(HERE, "support_localize.csv"), row.names = FALSE)
