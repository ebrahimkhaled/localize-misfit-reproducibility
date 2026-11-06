# burn_iecv2.R -- internal-external cross-validation of the burn models, with the whole pipeline inside the simulation.
# Risks from a model fitted without a patient's facility depend on the outcomes of the other facilities, so given the
# covariates they carry information on the patient's own facility. Two consequences are handled here.
# (1) baseline: outcomes simulated from each model fitted to all patients (form right by construction), the 40 refits
#     repeated, and the verdict computed: the error rate of the verdict on these pooled risks.
# (2) sizes: the chance level of each part comes from outcomes drawn from the simpler update with the 40 refits repeated on
#     each draw, so the chance level contains the same dependence as the observed size.
# Output burn_iecv_baseline.csv, burn_iecv_size.csv.
suppressMessages({ library(splines); library(aplore3) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
b <- aplore3::burn1000
bd <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                 white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"), facility = b$facility)
VB <- c("age", "tbsa", "male", "white", "inh", "flame"); X <- as.matrix(bd[, VB]); G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
FB <- list(B0 = y ~ age + tbsa + male + white + inh + flame, B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)
fac <- sort(unique(bd$facility))
iecv <- function(m, yy) { dd <- bd; dd$y <- yy; p <- numeric(nrow(dd))
  for (h in fac) { i <- dd$facility == h; f <- suppressWarnings(glm(FB[[m]], data = dd[!i, ], family = binomial())); p[i] <- predict(f, dd[i, ], type = "response") }
  pmin(pmax(p, 1e-6), 1 - 1e-6) }
base_one <- function(k) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  m <- names(FB)[1 + k %% 2]; set.seed(62000000L + k)
  pt <- fitted(glm(FB[[m]], data = bd, family = binomial())); yy <- rbinom(nrow(bd), 1L, pt)
  r <- localize_external_update(yy, iecv(m, yy), X, B = 199L, main = FALSE)
  data.frame(model = m, t(setNames(G %in% r$named, G)), t(setNames(as.numeric(r$intersection[G]), paste0("p_", G))))
}
fits <- function(yy, p, C) { e <- qlogis(p)
  pI <- fitted(glm(yy ~ 1, offset = e, family = binomial())); pS <- fitted(glm(yy ~ e, family = binomial())); S4 <- ns(e, df = 4)
  pL <- fitted(glm(yy ~ S4, family = binomial())); pC <- fitted(suppressWarnings(glm(yy ~ S4 + C, family = binomial())))
  list(p, pI, pS, pL, pC) }
sizes <- function(f) vapply(1:4, function(k) 100 * mean(abs(f[[k + 1]] - f[[k]])), 0)
size_one <- function(j) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R"))
  m <- SJ$model[j]; k <- SJ$part[j]; set.seed(63000000L + j)
  p <- P[[m]]; C <- CC[[m]]; base <- suppressWarnings(fits(bd$y, p, C))[[k]]
  yy <- rbinom(nrow(bd), 1L, base); ps <- iecv(m, yy); Cs <- do.call(cbind, .group_bases(qlogis(ps), X)$COV)
  data.frame(model = m, part = k, size = tryCatch(suppressWarnings(sizes(fits(yy, ps, Cs)))[k], error = function(e) NA_real_))
}
source(file.path(HERE, "R", "localize_groups.R"))
P <- lapply(setNames(names(FB), names(FB)), function(m) iecv(m, bd$y))
CC <- lapply(P, function(p) do.call(cbind, .group_bases(qlogis(p), X)$COV))
SJ <- expand.grid(rep = 1:200, part = 1:4, model = names(FB), stringsAsFactors = FALSE)
cl <- parallel::makeCluster(21L)
parallel::clusterExport(cl, c("HERE", "bd", "X", "G", "FB", "fac", "iecv", "base_one", "fits", "sizes", "size_one", "P", "CC", "SJ"))
d <- do.call(rbind, parallel::parLapplyLB(cl, 1:600, base_one, chunk.size = 3L))
s <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(SJ)), size_one, chunk.size = 4L)); parallel::stopCluster(cl)
utils::write.csv(d, file.path(HERE, "burn_iecv_baseline_p.csv"), row.names = FALSE)
bl <- do.call(rbind, lapply(split(d, d$model), function(g) data.frame(model = g$model[1], reps = nrow(g), t(colMeans(g[, G])), any = mean(rowSums(g[, G]) > 0))))
utils::write.csv(bl, file.path(HERE, "burn_iecv_baseline.csv"), row.names = FALSE); print(bl, row.names = FALSE, digits = 3)
res <- do.call(rbind, lapply(names(FB), function(m) { o <- suppressWarnings(sizes(fits(bd$y, P[[m]], CC[[m]])))
  do.call(rbind, lapply(1:4, function(k) { v <- s$size[s$model == m & s$part == k]; v <- v[!is.na(v)]
    data.frame(model = m, part = G[k], size = o[k], chance = mean(v), excess = o[k] - mean(v), lo = o[k] - quantile(v, .975, names = FALSE), hi = o[k] - quantile(v, .025, names = FALSE), draws = length(v)) })) }))
utils::write.csv(res, file.path(HERE, "burn_iecv_size.csv"), row.names = FALSE); print(res, row.names = FALSE, digits = 3)
## the verdicts on the observed risks, with and without the centre effect
source(file.path(HERE, "R", "localize_update.R"))
out <- do.call(rbind, lapply(names(FB), function(m) {
  e <- qlogis(P[[m]]); dd <- data.frame(y = bd$y, e = e, facility = factor(bd$facility))
  g0 <- lme4::glmer(y ~ 1 + offset(e) + (1 | facility), data = dd, family = binomial()); tau <- as.numeric(attr(lme4::VarCorr(g0)$facility, "stddev"))
  g1 <- suppressWarnings(lme4::glmer(y ~ e + (1 | facility), data = dd, family = binomial())); ci <- suppressMessages(confint.default(glm(y ~ e, data = dd, family = binomial())))
  set.seed(20261120L); a <- localize_external_update(bd$y, P[[m]], X, B = 999L, main = FALSE)
  set.seed(20261120L); b <- localize_external_update(bd$y, P[[m]], X, B = 999L, cluster = bd$facility, tau = tau, main = FALSE)
  data.frame(model = m, tau = tau, cal_int_re = lme4::fixef(g0)[1], cal_slope_re = lme4::fixef(g1)[2], slope_lo = ci[2, 1], slope_hi = ci[2, 2],
             named_plain = paste(a$named, collapse = " + "), named_cluster = paste(b$named, collapse = " + "), t(setNames(as.numeric(a$intersection[G]), paste0("plain_", G))), t(setNames(as.numeric(b$intersection[G]), G)))
}))
utils::write.csv(out, file.path(HERE, "burn_iecv.csv"), row.names = FALSE); print(out, row.names = FALSE, digits = 3)
cat("done", format(Sys.time()), "\n")
