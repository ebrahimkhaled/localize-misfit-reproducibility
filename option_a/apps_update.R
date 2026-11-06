# apps_update.R -- the two applications with the updated-model reference (R/localize_update.R), Holm, B = 999, on the
# data, models and splits of support_v2.R, support_split_first.R and burn_revise.R. Adds what the Statistics in Medicine
# version needs: (1) SUPPORT: random split, transported, half A, half B after each update, neural network; case mix
# against effect differences (membership model, spread of the score); (2) burn: 20 facility splits; correct-model
# baseline; internal-external cross-validation (each patient's risk from a model fitted without that patient's
# facility); facility variance; a cluster-aware reference; random-effects calibration across facilities.
# Output apps_update.csv, burn_update_splits.csv, burn_update_null.csv, burn_iecv.csv, support_casemix.csv.
suppressMessages({ library(splines); library(aplore3); library(nnet) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
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

jobs <- list()
add <- function(data, model, setting, y, p, X) jobs[[length(jobs) + 1]] <<- list(data = data, model = model, setting = setting, y = y, p = as.numeric(p), X = X)
set.seed(20261006L); dev_r <- sort(sample(nrow(d), floor(nrow(d) / 2)))
v <- d[d$acute == 0, ]
set.seed(20261011L); A <- sort(sample(nrow(v), floor(nrow(v) / 2))); Bh <- setdiff(seq_len(nrow(v)), A)
for (m in names(FS)) {
  f <- suppressWarnings(glm(FS[[m]], data = d[dev_r, ], family = binomial()))
  add("SUPPORT", m, "random split", d$y[-dev_r], predict(f, d[-dev_r, ], type = "response"), as.matrix(d[-dev_r, VS]))
  f <- suppressWarnings(glm(FS[[m]], data = d[d$acute == 1, ], family = binomial()))
  p <- as.numeric(predict(f, v, type = "response")); e <- qlogis(p)
  add("SUPPORT", m, "transported", v$y, p, as.matrix(v[, VS]))
  add("SUPPORT", m, "transported, half A", v$y[A], p[A], as.matrix(v[A, VS]))
  fi <- glm(v$y[A] ~ offset(e[A]), family = binomial())
  add("SUPPORT", m, "half B, intercept update", v$y[Bh], plogis(coef(fi)[1] + e[Bh]), as.matrix(v[Bh, VS]))
  fr <- glm(v$y[A] ~ e[A], family = binomial())
  add("SUPPORT", m, "half B, logistic recalibration", v$y[Bh], plogis(coef(fr)[1] + coef(fr)[2] * e[Bh]), as.matrix(v[Bh, VS]))
  ff <- glm(v$y[A] ~ ns(e[A], df = 3), family = binomial())
  add("SUPPORT", m, "half B, flexible recalibration", v$y[Bh],
      plogis(cbind(1, predict(ns(e[A], df = 3), e[Bh])) %*% coef(ff)), as.matrix(v[Bh, VS]))
  fv <- suppressWarnings(glm(FS[[m]], data = v[A, ], family = binomial()))
  add("SUPPORT", m, "half B, model re-estimated", v$y[Bh], predict(fv, v[Bh, ], type = "response"), as.matrix(v[Bh, VS]))
}
## the neural network of support_v2.R
Zs <- function(q, ref) scale(as.matrix(q[, VS]), center = colMeans(ref[, VS]), scale = apply(ref[, VS], 2, sd))
for (sp in c("random split", "transported")) {
  dv <- if (sp == "random split") d[dev_r, ] else d[d$acute == 1, ]
  vv <- if (sp == "random split") d[-dev_r, ] else v
  set.seed(20261010L)
  nn <- nnet::nnet(Zs(dv, dv), dv$y, size = 5, decay = 0.01, maxit = 1000, entropy = TRUE, trace = FALSE)
  add("SUPPORT", "neural network", sp, vv$y, predict(nn, Zs(vv, dv)), as.matrix(vv[, VS]))
}
## burn: pooled internal-external cross-validation, each patient's risk from a model without that patient's facility
fac <- sort(unique(bd$facility))
iecv <- lapply(FB, function(fm) {
  p <- numeric(nrow(bd))
  for (h in fac) { i <- bd$facility == h
    f <- suppressWarnings(glm(fm, data = bd[!i, ], family = binomial())); p[i] <- predict(f, bd[i, ], type = "response") }
  p
})
for (m in names(FB)) add("burn", m, "internal-external", bd$y, iecv[[m]], as.matrix(bd[, VB]))

one <- function(k) {
  suppressMessages({ library(splines) })
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  j <- jobs[[k]]; set.seed(20261105L + k)
  r <- localize_external_update(j$y, j$p, j$X, B = 999L)
  data.frame(data = j$data, model = j$model, setting = j$setting, n = length(j$y), events = sum(j$y),
             observed = mean(j$y), predicted = mean(j$p), named = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
             t(stats::setNames(as.numeric(r$intersection[G]), G)), t(stats::setNames(as.numeric(r$adjusted[G]), paste0("adj_", G))))
}
cl <- parallel::makeCluster(7L)
parallel::clusterExport(cl, c("HERE", "jobs", "one", "G"))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_along(jobs), one, chunk.size = 1L))
utils::write.csv(res, file.path(HERE, "apps_update.csv"), row.names = FALSE)
print(res[, c("data", "model", "setting", "n", "events", "named", G)], row.names = FALSE, digits = 3)
flush.console()

## burn: the 20 facility splits and the correct-model baseline
splitrun <- function(k) {
  suppressMessages({ library(splines) })
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  do.call(rbind, lapply(names(FB), function(m) {
    set.seed(20262000L + k); dev_f <- sample(fac, length(fac) %/% 2)
    f <- suppressWarnings(glm(FB[[m]], data = bd[bd$facility %in% dev_f, ], family = binomial()))
    vb <- bd[!bd$facility %in% dev_f, ]; p <- as.numeric(predict(f, vb, type = "response"))
    set.seed(20263000L + k); r <- localize_external_update(vb$y, p, as.matrix(vb[, VB]), B = 499L)
    data.frame(split = k, model = m, n = nrow(vb), deaths = sum(vb$y), t(stats::setNames(G %in% r$named, paste0("named_", G))),
               t(stats::setNames(as.numeric(r$intersection[G]), paste0("p_", G))))
  }))
}
nulljobs <- expand.grid(model = names(FB), sim = 1:25, stringsAsFactors = FALSE)
nullrun <- function(k) {
  suppressMessages({ library(splines) })
  source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  m <- nulljobs$model[k]; sm <- nulljobs$sim[k]
  ptrue <- fitted(glm(FB[[m]], data = bd, family = binomial()))
  set.seed(20264000L + 100L * sm + (m == "B3")); dd <- bd; dd$y <- rbinom(nrow(bd), 1L, ptrue)
  do.call(rbind, lapply(1:20, function(sp) {
    set.seed(20262000L + sp); dev_f <- sample(fac, length(fac) %/% 2)
    f <- suppressWarnings(glm(FB[[m]], data = dd[dd$facility %in% dev_f, ], family = binomial()))
    vb <- dd[!dd$facility %in% dev_f, ]; p <- as.numeric(predict(f, vb, type = "response"))
    set.seed(20265000L + 1000L * sm + sp); r <- localize_external_update(vb$y, p, as.matrix(vb[, VB]), B = 199L)
    data.frame(model = m, sim = sm, split = sp, t(stats::setNames(G %in% r$named, paste0("named_", G))))
  }))
}
parallel::clusterExport(cl, c("bd", "FB", "VB", "fac", "nulljobs", "splitrun", "nullrun"))
sp <- do.call(rbind, parallel::parLapplyLB(cl, 1:20, splitrun, chunk.size = 1L))
utils::write.csv(sp, file.path(HERE, "burn_update_splits.csv"), row.names = FALSE)
nl <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(nulljobs)), nullrun, chunk.size = 1L))
utils::write.csv(nl, file.path(HERE, "burn_update_null.csv"), row.names = FALSE)
parallel::stopCluster(cl)
for (m in names(FB)) cat(m, "splits naming each part:", round(colMeans(sp[sp$model == m, paste0("named_", G)]), 2),
                         "| baseline:", round(colMeans(nl[nl$model == m, paste0("named_", G)]), 3), "\n")

## burn: facility variance, cluster-aware reference, calibration across facilities (internal-external risks)
out <- list()
for (m in names(FB)) {
  e <- qlogis(iecv[[m]]); dd <- data.frame(y = bd$y, e = e, facility = factor(bd$facility))
  g0 <- lme4::glmer(y ~ 1 + offset(e) + (1 | facility), data = dd, family = binomial())
  tau <- as.numeric(attr(lme4::VarCorr(g0)$facility, "stddev"))
  g1 <- suppressWarnings(lme4::glmer(y ~ e + (1 | facility), data = dd, family = binomial()))
  gl <- glm(y ~ e, data = dd, family = binomial()); ci <- suppressMessages(confint.default(gl))
  set.seed(20261120L)
  r <- localize_external_update(bd$y, iecv[[m]], as.matrix(bd[, VB]), B = 999L, cluster = bd$facility, tau = tau)
  out[[m]] <- data.frame(model = m, tau = tau, cal_int_re = lme4::fixef(g0)[1], cal_slope_re = lme4::fixef(g1)[2],
                         slope_lo = ci[2, 1], slope_hi = ci[2, 2],
                         named_cluster = if (length(r$named)) paste(r$named, collapse = " + ") else "none",
                         t(stats::setNames(as.numeric(r$intersection[G]), G)))
}
out <- do.call(rbind, out); utils::write.csv(out, file.path(HERE, "burn_iecv.csv"), row.names = FALSE); print(out, row.names = FALSE, digits = 3)

## SUPPORT: case mix against effect differences (Debray et al. 2015): membership model and the spread of the score
fm <- glm(acute ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male + y, data = d, family = binomial())
cm <- data.frame(membership_c = as.numeric(pROC::auc(pROC::roc(d$acute, fitted(fm), quiet = TRUE))))
for (m in names(FS)) {
  f <- suppressWarnings(glm(FS[[m]], data = d[d$acute == 1, ], family = binomial()))
  lp1 <- predict(f, d[d$acute == 1, ]); lp0 <- predict(f, v); fv <- suppressWarnings(glm(FS[[m]], data = v, family = binomial()))
  cm[[paste0(m, "_mean_dev")]] <- mean(lp1); cm[[paste0(m, "_mean_val")]] <- mean(lp0)
  cm[[paste0(m, "_sd_dev")]] <- sd(lp1); cm[[paste0(m, "_sd_val")]] <- sd(lp0)
  cm[[paste0(m, "_c_dev")]] <- as.numeric(pROC::auc(pROC::roc(d$y[d$acute == 1], lp1, quiet = TRUE)))
  cm[[paste0(m, "_c_val")]] <- as.numeric(pROC::auc(pROC::roc(v$y, lp0, quiet = TRUE)))
}
utils::write.csv(cm, file.path(HERE, "support_casemix.csv"), row.names = FALSE); print(t(cm), digits = 3)
cat("done", format(Sys.time()), "\n")
