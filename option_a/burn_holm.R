# burn_holm.R -- the burn cross-centre analyses again, keeping each part's own p-value, so that Holm's correction (the
# default naming rule from the SiM version on) can be applied: (1) the 20 facility splits of burn_revise.R; (2) the
# correct-model baseline of burn_null.R. Same seeds, so the closure verdicts must equal the stored ones (checked).
# Output burn_holm_splits.csv, burn_holm_null.csv.
suppressMessages({ library(splines); library(aplore3) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
source(file.path(HERE, "R", "localize_groups.R"))
b <- aplore3::burn1000
d <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"),
                facility = b$facility)
VARS <- c("age", "tbsa", "male", "white", "inh", "flame")
FM <- list(B0 = y ~ age + tbsa + male + white + inh + flame,
           B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)
fac <- sort(unique(d$facility)); G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
holm_named <- function(p, alpha = 0.05) names(p)[stats::p.adjust(p, "holm") <= alpha]
rec <- function(r) {
  p <- r$intersection[G]
  c(stats::setNames(as.numeric(p), paste0("p_", G)),
    stats::setNames(G %in% r$named, paste0("closure_", G)),
    stats::setNames(G %in% holm_named(p), paste0("holm_", G)))
}

## (1) the 20 facility splits
spl <- do.call(rbind, lapply(1:20, function(k) do.call(rbind, lapply(c("B0", "B3"), function(m) {
  set.seed(20262000L + k); dev_f <- sample(fac, length(fac) %/% 2)
  dv <- d[d$facility %in% dev_f, ]; v <- d[!d$facility %in% dev_f, ]
  f <- suppressWarnings(glm(FM[[m]], data = dv, family = binomial()))
  p <- as.numeric(predict(f, newdata = v, type = "response"))
  set.seed(20263000L + k); r <- localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L)
  data.frame(split = k, model = m, n_val = nrow(v), deaths = sum(v$y), t(rec(r)))
}))))
old <- utils::read.csv(file.path(HERE, "burn_revise_splits.csv"))
for (g in G) stopifnot(all(as.logical(spl[[paste0("closure_", g)]]) == old[[paste0("named_", g)]]))
utils::write.csv(spl, file.path(HERE, "burn_holm_splits.csv"), row.names = FALSE)
cat("splits: closure verdicts match burn_revise_splits.csv\n")

## (2) the correct-model baseline (as burn_null.R)
jobs <- expand.grid(model = names(FM), sim = 1:25, stringsAsFactors = FALSE)
one <- function(k) {
  suppressMessages({ library(splines) }); source(file.path(HERE, "R", "localize_groups.R"))
  m <- jobs$model[k]; s <- jobs$sim[k]
  ptrue <- fitted(glm(FM[[m]], data = d, family = binomial()))
  set.seed(20264000L + 100L * s + (m == "B3")); dd <- d; dd$y <- rbinom(nrow(d), 1L, ptrue)
  do.call(rbind, lapply(1:20, function(sp) {
    set.seed(20262000L + sp); dev_f <- sample(fac, length(fac) %/% 2)
    dv <- dd[dd$facility %in% dev_f, ]; v <- dd[!dd$facility %in% dev_f, ]
    f <- suppressWarnings(glm(FM[[m]], data = dv, family = binomial()))
    p <- as.numeric(predict(f, newdata = v, type = "response"))
    set.seed(20265000L + 1000L * s + sp); r <- localize_external(v$y, p, as.matrix(v[, VARS]), M = 999L)
    data.frame(model = m, sim = s, split = sp, t(rec(r)))
  }))
}
cl <- parallel::makeCluster(6L)
parallel::clusterExport(cl, c("HERE", "jobs", "one", "d", "FM", "fac", "G", "VARS", "rec", "holm_named"))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(jobs)), one, chunk.size = 1L)); parallel::stopCluster(cl)
old <- utils::read.csv(file.path(HERE, "burn_null.csv"))
res <- res[order(res$model, res$sim, res$split), ]; old <- old[order(old$model, old$sim, old$split), ]
for (g in G) stopifnot(all(as.logical(res[[paste0("closure_", g)]]) == old[[paste0("named_", g)]]))
utils::write.csv(res, file.path(HERE, "burn_holm_null.csv"), row.names = FALSE)
cat("baseline: closure verdicts match burn_null.csv\n")
for (rule in c("closure", "holm")) {
  cat("\n", rule, ": share naming COV  real splits:",
      sprintf("B0 %.2f  B3 %.2f", mean(as.logical(spl[spl$model == "B0", paste0(rule, "_COV")])),
              mean(as.logical(spl[spl$model == "B3", paste0(rule, "_COV")]))),
      " | correct-model baseline:",
      sprintf("B0 %.3f  B3 %.3f", mean(as.logical(res[res$model == "B0", paste0(rule, "_COV")])),
              mean(as.logical(res[res$model == "B3", paste0(rule, "_COV")]))), "\n")
}
