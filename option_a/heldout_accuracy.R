# heldout_accuracy.R -- does following the verdict give better risks on real held-out patients? SUPPORT: each update is fitted on a random
# half of the validation patients and scored on the other half (20 halves; the first is the paper's split). Burn: models developed
# on half of the facilities and scored on the other half (the 20 facility splits of the paper), and by internal-external cross-validation.
# Also the planning numbers of Corollary A.4 for the slope test with the weights of the intercept update.
# Output heldout_support.csv, heldout_burn.csv, plan_updated.csv.
suppressMessages({ library(splines); library(glmnet); library(aplore3); library(nnet) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
s <- read.csv("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data/support2.csv", na.strings = c("NaN", "NA", ""))
d <- na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
  crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"), acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
VS <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co", "male")
FS <- list(M0 = y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male,
           M2 = y ~ age + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) + ns(wblc, 3) + num.co + male)
b <- aplore3::burn1000
bd <- data.frame(y = as.integer(b$death == "Dead"), age = b$age, tbsa = b$tbsa, male = as.numeric(b$gender == "Male"),
                 white = as.numeric(b$race == "White"), inh = as.numeric(b$inh_inj == "Yes"), flame = as.numeric(b$flame == "Yes"), facility = b$facility)
FB <- list(B0 = y ~ age + tbsa + male + white + inh + flame, B3 = y ~ ns(age, 3) + ns(tbsa, 3) + male + white + inh + flame + age:inh)
fac <- sort(unique(bd$facility))
rev_basis <- function(Xv) {
  cont <- which(apply(Xv, 2, function(u) length(unique(u)) > 4)); mu <- colMeans(Xv); sdv <- apply(Xv, 2, sd)
  kn <- lapply(cont, function(j) ns(Xv[, j], df = 3))
  function(Xn) { Zs <- sweep(sweep(Xn, 2, mu), 2, sdv, "/")
    sq <- do.call(cbind, lapply(cont, function(j) cbind(Zs[, j]^2, Zs[, j]^3)))
    pr <- do.call(cbind, combn(ncol(Xn), 2, function(ij) Zs[, ij[1]] * Zs[, ij[2]], simplify = FALSE))
    sp <- do.call(cbind, lapply(seq_along(cont), function(k) predict(kn[[k]], Xn[, cont[k]]))); cbind(Zs, sq, pr, sp) } }
auc <- function(y, p) { r <- rank(p); n1 <- sum(y); n0 <- length(y) - n1; (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
met <- function(y, p) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); c(ll = -100 * mean(y * log(p) + (1 - y) * log(1 - p)), brier = 100 * mean((y - p)^2), c = auc(y, p)) }
## every update fitted on half A and scored on half B
upd <- function(seed, p, yy, X) {
  e <- qlogis(p); n <- length(yy)
  set.seed(seed); A <- sort(sample(n, floor(n / 2))); B <- setdiff(seq_len(n), A)
  yA <- yy[A]; yB <- yy[B]; XA <- X[A, ]; XB <- X[B, ]
  fi <- glm(yA ~ offset(e[A]), family = binomial()); fr <- glm(yA ~ e[A], family = binomial())
  N <- ns(e[A], df = 3); ff <- glm(yA ~ N, family = binomial())
  rb <- rev_basis(XA); CA <- rb(XA); CB <- rb(XB)
  set.seed(seed + 1); cv <- cv.glmnet(CA, yA, family = "binomial", alpha = 0, offset = e[A], nfolds = 5)
  pr <- list(frozen = p[B], intercept = plogis(coef(fi)[1] + e[B]), recalibration = plogis(coef(fr)[1] + coef(fr)[2] * e[B]),
             flexible = as.numeric(plogis(cbind(1, predict(N, e[B])) %*% coef(ff))),
             revision = as.numeric(predict(cv, CB, newoffset = e[B], s = "lambda.min", type = "response")))
  do.call(rbind, lapply(names(pr), function(k) data.frame(seed = seed, update = k, t(met(yB, pr[[k]])))))
}
set.seed(20261006L); dev_r <- sort(sample(nrow(d), floor(nrow(d) / 2)))
v <- d[d$acute == 0, ]; a <- d[d$acute == 1, ]
SET <- list("similar patients" = list(dev = d[dev_r, ], val = d[-dev_r, ]), transported = list(dev = a, val = v))
seeds <- c(20261011L, 1:19)
J <- expand.grid(seed = seeds, m = names(FS), st = names(SET), stringsAsFactors = FALSE)
one <- function(i) { j <- J[i, ]; S <- SET[[j$st]]
  f <- suppressWarnings(glm(FS[[j$m]], data = S$dev, family = binomial())); p <- as.numeric(predict(f, S$val, type = "response"))
  cbind(setting = j$st, model = j$m, upd(j$seed, p, S$val$y, as.matrix(S$val[, VS]))) }
cl <- parallel::makeCluster(8L); parallel::clusterEvalQ(cl, suppressMessages({ library(splines); library(glmnet) })) -> .z
parallel::clusterExport(cl, c("SET", "FS", "VS", "J", "upd", "rev_basis", "auc", "met", "one"))
R <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(J)), one, chunk.size = 1L)); parallel::stopCluster(cl)
write.csv(R, file.path(HERE, "heldout_support.csv"), row.names = FALSE)
print(aggregate(cbind(ll, brier, c) ~ setting + model + update, R, mean), digits = 4)
for (st in names(SET)) for (m in names(FS)) { w <- function(u) R$ll[R$setting == st & R$model == m & R$update == u]
  cat(st, m, "revision better than flexible", sum(w("revision") < w("flexible")), "than intercept", sum(w("revision") < w("intercept")), "of 20; mean gain over flexible", round(mean(w("flexible") - w("revision")), 3), "\n") }
## burn: frozen B0 and B3 on held-out facilities
BR <- do.call(rbind, lapply(1:20, function(k) do.call(rbind, lapply(names(FB), function(m) {
  set.seed(20262000L + k); dev_f <- sample(fac, length(fac) %/% 2)
  f <- suppressWarnings(glm(FB[[m]], data = bd[bd$facility %in% dev_f, ], family = binomial())); vb <- bd[!bd$facility %in% dev_f, ]
  data.frame(design = "facility split", split = k, model = m, t(met(vb$y, as.numeric(predict(f, vb, type = "response"))))) }))))
IE <- do.call(rbind, lapply(names(FB), function(m) { p <- numeric(nrow(bd))
  for (h in fac) { i <- bd$facility == h; f <- suppressWarnings(glm(FB[[m]], data = bd[!i, ], family = binomial())); p[i] <- predict(f, bd[i, ], type = "response") }
  data.frame(design = "internal-external", split = 0, model = m, t(met(bd$y, p))) }))
BB <- rbind(BR, IE); write.csv(BB, file.path(HERE, "heldout_burn.csv"), row.names = FALSE)
print(aggregate(cbind(ll, brier, c) ~ design + model, BB, mean), digits = 4)
w0 <- BR$ll[BR$model == "B0"]; w3 <- BR$ll[BR$model == "B3"]; cat("B3 better in", sum(w3 < w0), "of 20 splits; mean gain", round(mean(w0 - w3), 3), "\n")
## planning: the slope test with frozen weights (Web Table) and with the weights of the intercept update
eps <- function(p, wt) { e <- qlogis(p); m <- sum(wt * e) / sum(wt); 2.8016 / sqrt(sum(wt * (e - m)^2)) }
Zs <- function(q, ref) scale(as.matrix(q[, VS]), center = colMeans(ref[, VS]), scale = apply(ref[, VS], 2, sd))
P <- do.call(rbind, lapply(names(SET), function(st) { S <- SET[[st]]; ps <- lapply(FS, function(fm) as.numeric(predict(suppressWarnings(glm(fm, data = S$dev, family = binomial())), S$val, type = "response")))
  set.seed(20261010L); nn <- nnet::nnet(Zs(S$dev, S$dev), S$dev$y, size = 5, decay = 0.01, maxit = 1000, entropy = TRUE, trace = FALSE); ps[["neural network"]] <- as.numeric(predict(nn, Zs(S$val, S$dev)))
  do.call(rbind, lapply(names(ps), function(m) { p <- ps[[m]]; e <- qlogis(p); pI <- fitted(glm(S$val$y ~ 1, offset = e, family = binomial()))
    data.frame(setting = st, model = m, eps_frozen = eps(p, p * (1 - p)), eps_updated = eps(p, pI * (1 - pI))) })) }))
write.csv(P, file.path(HERE, "plan_updated.csv"), row.names = FALSE); print(P, digits = 3)
