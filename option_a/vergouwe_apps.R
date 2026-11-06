# vergouwe_apps.R -- Vergouwe et al.'s (2017) closed testing procedure on the two applications, beside our verdicts:
# SUPPORT (the same developments and validation sets as support_v2.R / support_split_first.R) and the burn data (the
# same 20 facility splits as burn_revise.R). For a spline model the 'revision' refits the model's own design columns
# (the splines keep their development knots). Output vergouwe_apps.csv.
suppressMessages({ library(splines); library(aplore3) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"

## the procedure for a fitted glm 'fit' applied to new data 'nd': lp_old is the frozen linear predictor, X the
## model's own design columns (no intercept) evaluated on nd
ctp <- function(fit, nd, alpha = 0.05) {
  y <- nd$y
  tt <- stats::delete.response(stats::terms(fit))
  X <- stats::model.matrix(tt, stats::model.frame(tt, nd))[, -1, drop = FALSE]
  lp_old <- as.numeric(stats::predict(fit, newdata = nd, type = "link"))
  ll <- function(lp) sum(y * lp - log1p(exp(lp)))
  int <- stats::coef(stats::glm(y ~ 1, offset = lp_old, family = stats::binomial()))
  rc <- stats::coef(stats::glm(y ~ lp_old, family = stats::binomial()))
  rv <- suppressWarnings(stats::glm(y ~ X, family = stats::binomial()))
  l0 <- ll(lp_old); l1 <- ll(int + lp_old); l2 <- ll(rc[1] + rc[2] * lp_old); l3 <- as.numeric(stats::logLik(rv))
  p <- ncol(X)
  pB <- stats::pchisq(2 * (l3 - l0), p + 1, lower.tail = FALSE)
  pC <- stats::pchisq(2 * (l3 - l1), p, lower.tail = FALSE)
  pD <- stats::pchisq(2 * (l3 - l2), p - 1, lower.tail = FALSE)
  choice <- if (pB >= alpha) "original" else if (pC >= alpha) "intercept update" else if (pD >= alpha) "recalibration" else "revision"
  data.frame(n = nrow(nd), deaths = sum(y), vergouwe = choice, pB = pB, pC = pC, pD = pD)
}

## ---- SUPPORT
s <- utils::read.csv("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data/support2.csv",
                     na.strings = c("NaN", "NA", ""))
d <- stats::na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                               crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                               acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
FS <- list(M0 = y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male,
           M2 = y ~ age + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) + ns(wblc, 3) + num.co + male)
v <- d[d$acute == 0, ]
set.seed(20261011L); A <- sort(sample(nrow(v), floor(nrow(v) / 2)))          # half A of support_split_first.R
out <- list()
for (m in names(FS)) {
  fit <- suppressWarnings(stats::glm(FS[[m]], data = d[d$acute == 1, ], family = stats::binomial()))
  out[[length(out) + 1]] <- cbind(data = "SUPPORT", model = m, setting = "transported, all chronic patients", ctp(fit, v))
  out[[length(out) + 1]] <- cbind(data = "SUPPORT", model = m, setting = "transported, half A", ctp(fit, v[A, ]))
}
## random split: development on one half, validation on the other (the split of support_v2.R)
set.seed(20261006L); dev_r <- sort(sample(nrow(d), floor(nrow(d) / 2)))
for (m in names(FS)) {
  fit <- suppressWarnings(stats::glm(FS[[m]], data = d[dev_r, ], family = stats::binomial()))
  out[[length(out) + 1]] <- cbind(data = "SUPPORT", model = m, setting = "random split", ctp(fit, d[-dev_r, ]))
}

## ---- burn: the 20 facility splits of burn_revise.R
b <- aplore3::burn1000
bd <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                 white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"),
                 facility = b$facility)
FB <- list(B0 = y ~ age + tbsa + male + white + inh + flame,
           B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)
fac <- sort(unique(bd$facility))
for (k in 1:20) for (m in names(FB)) {
  set.seed(20262000L + k); dev_f <- sample(fac, length(fac) %/% 2)
  fit <- suppressWarnings(stats::glm(FB[[m]], data = bd[bd$facility %in% dev_f, ], family = stats::binomial()))
  out[[length(out) + 1]] <- cbind(data = "burn", model = m, setting = sprintf("facility split %02d", k),
                                  ctp(fit, bd[!bd$facility %in% dev_f, ]))
}
res <- do.call(rbind, out)
utils::write.csv(res, file.path(HERE, "vergouwe_apps.csv"), row.names = FALSE)
print(res[res$data == "SUPPORT", ], row.names = FALSE)
cat("\nburn, share of the 20 facility splits for each choice:\n")
print(with(res[res$data == "burn", ], round(prop.table(table(model, vergouwe), 1), 2)))
