## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# support_timing.R -- wall time of the procedure on the SUPPORT data (one core): in-sample M0 with B = 199 bootstrap
# refits (n = 8,873), and external validation of the transported M0 with M = 999 Monte Carlo draws (n = 4,130).
suppressMessages(library(splines))
source(file.path(.ROOT, "option_a/R/localize_groups.R"))
s <- utils::read.csv(file.path(.ROOT, "data/support2.csv"),
                     na.strings = c("NaN", "NA", ""))
d <- stats::na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                               crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                               acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
VARS <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co", "male")
F0 <- y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male
set.seed(1)
fit <- suppressWarnings(glm(F0, data = d, family = binomial()))
t1 <- system.time(localize_insample(fit, as.matrix(d[, VARS]), B = 199L))[["elapsed"]]
fa <- suppressWarnings(glm(F0, data = d[d$acute == 1, ], family = binomial())); v <- d[d$acute == 0, ]
p <- as.numeric(predict(fa, newdata = v, type = "response"))
t2 <- system.time(localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L))[["elapsed"]]
cat(sprintf("in-sample M0 (n = %d, B = 199): %.0f s; external transported M0 (n = %d, M = 999): %.1f s\n", nrow(d), t1, nrow(v), t2))
