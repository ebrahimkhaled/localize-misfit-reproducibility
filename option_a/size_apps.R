# size_apps.R -- the size of each part (localize_size(), percentage points of risk) and the standard calibration
# measures, for the external analyses the paper reports: SUPPORT random split, SUPPORT transported (all chronic
# patients, half A, and half B after each update fitted on A; the splits of support_split_first.R), and the burn data
# over the 20 facility splits of burn_revise.R. Output size_apps.csv (one row per analysis and part).
suppressMessages({ library(splines); library(aplore3) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"

s <- utils::read.csv("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data/support2.csv",
                     na.strings = c("NaN", "NA", ""))
d <- stats::na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                               crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                               acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
VS <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co", "male")
FS <- list(M0 = y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male,
           M2 = y ~ age + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) + ns(wblc, 3) + num.co + male)
b <- aplore3::burn1000
bd <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                 white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"),
                 facility = b$facility)
VB <- c("age", "tbsa", "male", "white", "inh", "flame")
FB <- list(B0 = y ~ age + tbsa + male + white + inh + flame,
           B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)

## the list of analyses: data set, frozen predictions, covariates
jobs <- list()
add <- function(data, model, setting, y, p, X) jobs[[length(jobs) + 1]] <<- list(data = data, model = model, setting = setting, y = y, p = p, X = X)
set.seed(20261006L); dev_r <- sort(sample(nrow(d), floor(nrow(d) / 2)))
v <- d[d$acute == 0, ]
set.seed(20261011L); A <- sort(sample(nrow(v), floor(nrow(v) / 2))); Bh <- setdiff(seq_len(nrow(v)), A)
for (m in names(FS)) {
  f <- suppressWarnings(glm(FS[[m]], data = d[dev_r, ], family = binomial()))
  add("SUPPORT", m, "random split", d$y[-dev_r], predict(f, d[-dev_r, ], type = "response"), as.matrix(d[-dev_r, VS]))
  f <- suppressWarnings(glm(FS[[m]], data = d[d$acute == 1, ], family = binomial()))
  p <- as.numeric(predict(f, v, type = "response")); e <- qlogis(p)
  add("SUPPORT", m, "transported, all chronic", v$y, p, as.matrix(v[, VS]))
  add("SUPPORT", m, "transported, half A", v$y[A], p[A], as.matrix(v[A, VS]))
  ## updates fitted on A, predictions on B (as support_split_first.R)
  fi <- glm(v$y[A] ~ offset(e[A]), family = binomial())
  add("SUPPORT", m, "half B after intercept update", v$y[Bh], plogis(coef(fi)[1] + e[Bh]), as.matrix(v[Bh, VS]))
  fr <- glm(v$y[A] ~ e[A], family = binomial())
  add("SUPPORT", m, "half B after logistic recalibration", v$y[Bh], plogis(coef(fr)[1] + coef(fr)[2] * e[Bh]), as.matrix(v[Bh, VS]))
  ff <- glm(v$y[A] ~ ns(e[A], df = 3), family = binomial())
  add("SUPPORT", m, "half B after flexible recalibration", v$y[Bh],
      as.numeric(plogis(cbind(1, predict(ns(e[A], df = 3), e[Bh])) %*% coef(ff))), as.matrix(v[Bh, VS]))
}
fac <- sort(unique(bd$facility))
for (k in 1:20) for (m in names(FB)) {
  set.seed(20262000L + k); dev_f <- sample(fac, length(fac) %/% 2)
  f <- suppressWarnings(glm(FB[[m]], data = bd[bd$facility %in% dev_f, ], family = binomial()))
  vb <- bd[!bd$facility %in% dev_f, ]
  add("burn", m, sprintf("facility split %02d", k), vb$y, as.numeric(predict(f, vb, type = "response")), as.matrix(vb[, VB]))
}

one <- function(k) {
  suppressMessages(library(splines))
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_size.R"))
  j <- jobs[[k]]
  s <- localize_size(j$y, as.numeric(j$p), j$X, R = 400L, seed = 20261104L + k)
  data.frame(data = j$data, model = j$model, setting = j$setting, n = length(j$y), events = sum(j$y), s$parts,
             cal_int = s$calibration["intercept", 1], cal_int_lo = s$calibration["intercept", 2], cal_int_hi = s$calibration["intercept", 3],
             cal_slope = s$calibration["slope", 1], cal_slope_lo = s$calibration["slope", 2], cal_slope_hi = s$calibration["slope", 3],
             ici = s$ici)
}
cl <- parallel::makeCluster(8L); parallel::clusterExport(cl, c("HERE", "jobs", "one"))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_along(jobs), one, chunk.size = 1L)); parallel::stopCluster(cl)
utils::write.csv(res, file.path(HERE, "size_apps.csv"), row.names = FALSE)
print(res[res$data == "SUPPORT", c("model", "setting", "part", "size", "chance", "excess", "lo", "hi", "cal_int", "cal_slope", "ici")],
      row.names = FALSE, digits = 3)
