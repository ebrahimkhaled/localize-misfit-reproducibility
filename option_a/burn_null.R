## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# burn_null.R -- how often does the cross-centre validation name COV when the model is RIGHT? Outcomes are simulated from
# B0 (and from B3) fitted to all 1,000 burn patients, so the model form is correct by construction; each simulated data
# set is then put through the same 20 facility splits as burn_revise.R (develop on half of the facilities, validate the
# frozen model on the other half, M = 999). The share of splits naming COV here is what the development half's own
# estimation error produces; compare it with the 70% (B0) and 75% (B3) of the real outcomes. Output burn_null.csv.
suppressMessages({ library(splines); library(aplore3) })
HERE <- file.path(.ROOT, "option_a")
source(file.path(HERE, "R", "localize_groups.R"))
b <- aplore3::burn1000
d <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"),
                facility = b$facility)
VARS <- c("age", "tbsa", "male", "white", "inh", "flame")
FM <- list(B0 = y ~ age + tbsa + male + white + inh + flame,
           B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)
fac <- sort(unique(d$facility)); G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
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
    data.frame(model = m, sim = s, split = sp, t(stats::setNames(G %in% r$named, paste0("named_", G))))
  }))
}
cl <- parallel::makeCluster(6L); parallel::clusterExport(cl, c("HERE", "jobs", "one", "d", "FM", "fac", "G", "VARS"))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(jobs)), one, chunk.size = 1L)); parallel::stopCluster(cl)
utils::write.csv(res, file.path(HERE, "burn_null.csv"), row.names = FALSE)
print(aggregate(cbind(named_INTERCEPT, named_SLOPE, named_LINK, named_COV) ~ model, res, mean))
