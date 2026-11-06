# check_audit.R -- stress tests of the updated-model reference with the final code (R/localize_update.R):
# coef     one predictor's effect changes: COV named against closed testing's revision (five settings, n = 1,000 and 4,000)
# ladder   exact members of the ladder with extreme parameters: any part above the simplest sufficient update is an error
# outside  a function of the score outside the four-df spline family: how often COV is named
# few      a correct model with few events
# dev      a model of correct form estimated on a development sample
# sizetrue the COV size against its true value, omitted interaction
# Output theory/CHECK_AUDIT.csv, theory/CHECK_SIZE_TRUE.csv.
suppressMessages(library(splines))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
gen_X <- function(design, n) {
  if (design == "indep5") X <- matrix(rnorm(n * 5), n, 5) else {
    Z <- matrix(rnorm(n * 5), n, 5) %*% chol(0.5 + 0.5 * diag(5)); Bn <- matrix(rbinom(n * 3, 1L, 0.3), n, 3)
    L <- exp(matrix(rnorm(n * 2, 0, 0.75), n, 2)); L <- (L - exp(0.75^2 / 2)) / sqrt((exp(0.75^2) - 1) * exp(0.75^2)); X <- cbind(Z, Bn, L) }
  colnames(X) <- paste0("x", seq_len(ncol(X))); X }
BETA <- list(indep5 = c(0.5, 0.8, 0.6, 0.4, 0.3), mixed10 = c(0.5, 0.6, 0.4, 0.3, 0.2, 0.5, 0.4, -0.3, 0.4, 0.3))
set.seed(1); A0 <- sapply(c("indep5", "mixed10"), function(d) { X <- gen_X(d, 2e5); lp <- as.numeric(X %*% BETA[[d]]); uniroot(function(a) mean(plogis(a + lp)) - .2, c(-10, 5))$root })
ver <- function(lp, X, y, alpha = .05) { ll <- function(l) sum(y * l - log1p(exp(l)))
  i <- coef(glm(y ~ 1, offset = lp, family = binomial())); r <- coef(glm(y ~ lp, family = binomial())); v <- suppressWarnings(glm(y ~ X, family = binomial()))
  l0 <- ll(lp); l1 <- ll(i + lp); l2 <- ll(r[1] + r[2] * lp); l3 <- as.numeric(logLik(v)); p <- ncol(X)
  if (pchisq(2 * (l3 - l0), p + 1, lower.tail = FALSE) >= alpha) "none" else if (pchisq(2 * (l3 - l1), p, lower.tail = FALSE) >= alpha) "intercept"
  else if (pchisq(2 * (l3 - l2), p - 1, lower.tail = FALSE) >= alpha) "recalibration" else "revision" }
jobs <- list(); addj <- function(exp, tag, reps, ...) for (r in seq_len(reps)) jobs[[length(jobs) + 1]] <<- list(exp = exp, tag = tag, r = r, a = list(...))
for (s in list(list("binary x6 +0.8", "mixed10", 6, .8), list("binary x7 +0.8", "mixed10", 7, .8), list("binary x6 +0.5", "mixed10", 6, .5),
               list("x2 effect removed", "indep5", 2, -.8), list("x1 +0.5", "mixed10", 1, .5)))
  for (n in c(1000L, 4000L)) addj("coef", paste(s[[1]], n), 100, design = s[[2]], j = s[[3]], d = s[[4]], n = n)
for (sc in c("flat", "nearflat", "steep_spline", "extreme_transport")) addj("ladder", sc, 300, sc = sc, n = 1000L)
for (sh in c("z2", "step", "kink")) for (n in c(2000L, 8000L)) addj("outside", paste(sh, n), 100, shape = sh, n = n)
for (s in list(c(300, .1), c(500, .1), c(1000, .05))) addj("few", paste(s, collapse = " "), 300, n = s[1], prev = s[2])
for (de in c("indep5", "mixed10")) for (nd in c(300L, 1000L)) addj("dev", paste(de, nd), 100, design = de, nd = nd, nv = 1000L)
one <- function(j) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_update.R"))
  set.seed(61000000L + j$r * 131L + sum(utf8ToInt(paste(j$exp, j$tag))) * 7L); a <- j$a; v <- NA_character_
  tryCatch({
    if (j$exp == "coef") { X <- gen_X(a$design, a$n); e <- as.numeric(A0[[a$design]] + X %*% BETA[[a$design]]); y <- rbinom(a$n, 1L, plogis(e + a$d * X[, a$j])); v <- ver(e, X, y) }
    if (j$exp == "ladder") { X <- gen_X("indep5", a$n); e <- as.numeric(-.5 + X %*% BETA$indep5); N4 <- ns(e, df = 4)
      lp <- switch(a$sc, flat = rep(-1, a$n), nearflat = -1 + .1 * e, steep_spline = as.numeric(-1 + N4 %*% c(4, -6, 6, -3)), extreme_transport = -3 + .3 * e); y <- rbinom(a$n, 1L, plogis(lp)) }
    if (j$exp == "outside") { X <- gen_X("indep5", a$n); e <- as.numeric(-.5 + X %*% BETA$indep5); z <- (e - mean(e)) / sd(e)
      y <- rbinom(a$n, 1L, plogis(e + switch(a$shape, z2 = .25 * (z^2 - 1), step = 1.5 * (z > 1), kink = 2 * pmax(z - 1, 0)))) }
    if (j$exp == "few") { X <- gen_X("indep5", a$n); lp <- as.numeric(X %*% BETA$indep5); aa <- uniroot(function(q) mean(plogis(q + lp)) - a$prev, c(-10, 5))$root; e <- aa + lp; y <- rbinom(a$n, 1L, plogis(e)) }
    if (j$exp == "dev") { b <- BETA[[a$design]]; aa <- A0[[a$design]]; Xd <- gen_X(a$design, a$nd); yd <- rbinom(a$nd, 1L, plogis(as.numeric(aa + Xd %*% b)))
      f <- suppressWarnings(glm(yd ~ Xd, family = binomial())); X <- gen_X(a$design, a$nv); tr <- as.numeric(aa + X %*% b); e <- as.numeric(cbind(1, X) %*% coef(f)); y <- rbinom(a$nv, 1L, plogis(tr)); v <- ver(e, X, y) }
    r <- localize_external_update(y, plogis(e), X, B = 199L)
    data.frame(exp = j$exp, tag = j$tag, t(setNames(G %in% r$named, G)), ver = v, err = "")
  }, error = function(x) data.frame(exp = j$exp, tag = j$tag, INTERCEPT = NA, SLOPE = NA, LINK = NA, COV = NA, ver = NA, err = conditionMessage(x)))
}
sz <- function(k) {
  suppressMessages(library(splines)); source(file.path(HERE, "R", "localize_groups.R")); source(file.path(HERE, "R", "localize_size.R"))
  b0 <- c(-.5, .5, .8, .6, .4, .3); cc <- SJ$c[k]; n <- SJ$n[k]; set.seed(100 + k); X <- matrix(rnorm(n * 5), n, 5); colnames(X) <- paste0("x", 1:5)
  e <- as.numeric(cbind(1, X) %*% b0); y <- rbinom(n, 1L, plogis(e + cc * X[, 1] * X[, 3])); s <- suppressWarnings(localize_size(y, plogis(e), X, R = 60L)); q <- s$parts[4, ]
  data.frame(c = cc, n = n, excess = q$excess, lo = q$lo, hi = q$hi)
}
SJ <- expand.grid(rep = 1:100, c = c(.2, .4, .8), n = c(1000L, 4000L))
truth <- sapply(c(.2, .4, .8), function(cc) { b0 <- c(-.5, .5, .8, .6, .4, .3); set.seed(5); n <- 4e5; X <- matrix(rnorm(n * 5), n, 5); e <- as.numeric(cbind(1, X) %*% b0)
  pi <- plogis(e + cc * X[, 1] * X[, 3]); N <- ns(e, df = 4); f <- suppressWarnings(glm(pi ~ N, family = quasibinomial())); 100 * mean(abs(pi - fitted(f))) })
cl <- parallel::makeCluster(as.integer(Sys.getenv("SG_WORKERS", "20"))); parallel::clusterExport(cl, c("HERE", "G", "gen_X", "BETA", "A0", "ver", "one", "sz", "SJ"))
d <- do.call(rbind, parallel::parLapplyLB(cl, jobs[sample(length(jobs))], one, chunk.size = 4L))
s <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(SJ)), sz, chunk.size = 4L)); parallel::stopCluster(cl)
print(table(d$err)); ok <- d[d$err == "", ]
res <- do.call(rbind, lapply(split(ok, paste(ok$exp, ok$tag)), function(g) data.frame(exp = g$exp[1], tag = g$tag[1], reps = nrow(g), t(colMeans(g[, G])), any = mean(rowSums(g[, G]) > 0),
  ctp_revision = if (all(is.na(g$ver))) NA else mean(g$ver == "revision"))))
utils::write.csv(res[order(res$exp, res$tag), ], file.path(HERE, "theory", "CHECK_AUDIT.csv"), row.names = FALSE); print(res[order(res$exp, res$tag), ], row.names = FALSE, digits = 3)
s$true <- truth[match(s$c, c(.2, .4, .8))]
z <- aggregate(cbind(true, excess, lo, hi) ~ c + n, s, mean); utils::write.csv(z, file.path(HERE, "theory", "CHECK_SIZE_TRUE.csv"), row.names = FALSE); print(z, digits = 3)
