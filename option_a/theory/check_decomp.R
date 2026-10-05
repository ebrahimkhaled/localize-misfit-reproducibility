## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# check_decomp.R -- Proposition A.1 and Remark A.1, checked on the SAME data sets the simulations used.
# For each data set we know the true risks pi, so each member's noncentrality on that design is exact:
#   lambda = Q evaluated at the residual vector pi - p0  (U = B'(pi - p0), V = B' W B, Q = U' V^- U).
# A group whose largest member noncentrality is well above 0 really HAS a part of the misfit on the probability
# scale, so naming it is a correct detection, not a familywise error. We report, per cell and group, the mean
# largest noncentrality, the implied power of that member alone at .05, and the rate the procedure named the group.
# (1) PURITY_PROB cells (run_purity_prob.R, seeds 52000000 + id * 10000 + r), first 100 data sets per cell;
# (2) the external fixed-size pure cells of run_sim_groups.R (seeds 41000000 + id * 10000 + r), first 100 per cell;
# (3) Remark A.1: size of the spline approximation of functions of eta for the COV bases (n = 100,000).
suppressMessages(library(splines))
HERE <- file.path(.ROOT, "option_a")
source(file.path(HERE, "R", "localize_groups.R"))
b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")

## noncentrality of every member at the true residual; returns max per group and the member df
ncp_groups <- function(X, p0, pt) {
  eta0 <- qlogis(p0); w0 <- p0 * (1 - p0); gb <- .group_bases(eta0, X)
  one <- matrix(1, length(p0), 1); Z2 <- cbind(one, eta0); Z3 <- cbind(Z2, ns(eta0, df = 5))
  bases <- list(INTERCEPT = list(one), SLOPE = list(.wperp(matrix(eta0), one, w0)),
                LINK = lapply(gb$LINK, function(B) .wperp(B, Z2, w0)), COV = lapply(gb$COV, function(B) .wperp(B, Z3, w0)))
  r <- pt - p0
  sapply(G, function(g) {
    v <- sapply(bases[[g]], function(B) {
      B <- as.matrix(B); B <- B[, sqrt(colSums(B^2)) > 1e-8, drop = FALSE]
      V <- crossprod(B, B * w0); e <- eigen(V, symmetric = TRUE); k <- sum(e$values > 1e-9 * max(e$values))
      U <- crossprod(B, r); z <- crossprod(e$vectors[, 1:k, drop = FALSE], U)
      lam <- sum(z^2 / e$values[1:k])
      c(lam = lam, pow = pchisq(qchisq(.95, k), k, ncp = lam, lower.tail = FALSE))
    })
    v[, which.max(v["pow", ])]
  })
}

out <- list()
## (1) purity on the probability scale
CP <- data.frame(fam = rep(c("slope_p", "link_p", "cov_p"), each = 3), C = c(.6, 1.2, 2.4, .5, 1, 2, .5, 1, 2)); CP$id <- seq_len(nrow(CP))
rows_prob <- do.call(rbind, lapply(list.files(file.path(HERE, "rows_prob"), "\\.csv$", full.names = TRUE), read.csv))
for (q in seq_len(nrow(CP))) {
  A <- t(sapply(1:100, function(r) {
    n <- 1000L; set.seed(52000000L + CP$id[q] * 10000L + r)
    X <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- plogis(eta0); w0 <- p0 * (1 - p0)
    perp <- function(h, Z) { g <- as.numeric(.wperp(matrix(h), Z, w0)); g / sqrt(weighted.mean(g^2, w0)) }
    g <- switch(CP$fam[q], slope_p = perp(eta0, matrix(1, n, 1)), link_p = perp(eta0^2, cbind(1, eta0)),
                cov_p = perp(X[, 1] * X[, 3], cbind(1, eta0, ns(eta0, df = 5))))
    pt <- pmin(pmax(p0 + CP$C[q] * w0 * g, 0.001), 0.999)
    m <- ncp_groups(X, p0, pt); c(m["lam", ], m["pow", ])
  }))
  nm <- rows_prob[rows_prob$fam == CP$fam[q] & abs(rows_prob$C - CP$C[q]) < 1e-9 & rows_prob$rep <= 100, ]
  named <- colMeans(sapply(paste0("named_", G), function(v) as.logical(nm[[v]])))
  out[[length(out) + 1]] <- data.frame(study = "probability-scale", fam = CP$fam[q], C = CP$C[q], group = G,
                                       mean_ncp = colMeans(A[, 1:4]), member_power = colMeans(A[, 5:8]), named = named)
}
## (2) logit-scale pure cells of the main simulation (external, fixed size)
ext <- data.frame(fam = c("null", "intercept", "intercept", "slope", "slope", "slope_pure", "link", "link_pure", "link_pure",
                          "ushape", "thresh", "inter", "cov_pure", "cov_pure"), C = c(0, .2, .4, .8, .6, .6, .25, .25, .5, .8, 4, 1, 1, 2))
ins <- data.frame(fam = c("null", "link", "link_pure", "link_pure", "ushape", "thresh", "inter", "cov_pure", "cov_pure"), C = c(0, .25, .25, .5, .8, 4, 1, 1, 2))
CELLS <- rbind(merge(cbind(setting = "external", ext), data.frame(n = c(500L, 1000L))),
               merge(cbind(setting = "insample", ins), data.frame(n = c(500L, 1000L))))
CELLS$id <- seq_len(nrow(CELLS))
rows <- do.call(rbind, lapply(list.files(file.path(HERE, "rows"), "^external.*\\.csv$", full.names = TRUE), read.csv))
pick <- which(CELLS$setting == "external" & CELLS$n == 1000L & CELLS$fam %in% c("slope_pure", "link_pure", "cov_pure"))
for (q in pick) {
  cc <- CELLS[q, ]
  A <- t(sapply(1:100, function(r) {
    n <- cc$n; set.seed(41000000L + cc$id * 10000L + r)
    X <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- plogis(eta0); w0 <- p0 * (1 - p0)
    po <- function(h, Z) as.numeric(.wperp(matrix(h), Z, w0))
    d <- switch(cc$fam, slope_pure = (cc$C - 1) * po(eta0, matrix(1, n, 1)), link_pure = cc$C * po(eta0^2, cbind(1, eta0)),
                cov_pure = cc$C * po(X[, 1] * X[, 3], cbind(1, eta0, ns(eta0, df = 5))))
    m <- ncp_groups(X, p0, plogis(eta0 + d)); c(m["lam", ], m["pow", ])
  }))
  nm <- rows[rows$fam == cc$fam & abs(rows$C - cc$C) < 1e-9 & rows$n == cc$n & rows$rep <= 100, ]
  named <- colMeans(sapply(paste0("named_", G), function(v) as.logical(nm[[v]])))
  out[[length(out) + 1]] <- data.frame(study = "logit-scale", fam = cc$fam, C = cc$C, group = G,
                                       mean_ncp = colMeans(A[, 1:4]), member_power = colMeans(A[, 5:8]), named = named)
}
D <- do.call(rbind, out); rownames(D) <- NULL
write.csv(D, file.path(HERE, "theory", "DECOMP.csv"), row.names = FALSE)
print(D, digits = 3)

## (3) spline approximation of E_w(b | eta) for the COV bases
set.seed(7); n <- 1e5
X <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
eta0 <- as.numeric(cbind(1, X) %*% b0); p0 <- plogis(eta0); w0 <- p0 * (1 - p0)
gb <- .group_bases(eta0, X)
Z5 <- cbind(1, eta0, ns(eta0, df = 5)); Z40 <- cbind(1, eta0, ns(eta0, df = 40))
wn <- function(v) sqrt(colSums(v^2 * w0))
ap <- do.call(rbind, lapply(names(gb$COV), function(m) {
  B <- gb$COV[[m]]; b5 <- .wperp(B, Z5, w0); b40 <- .wperp(B, Z40, w0)
  rel <- wn(b5 - b40) / wn(b5)
  data.frame(member = m, columns = ncol(B), median_rel = median(rel), max_rel = max(rel))
}))
cat("\nRemark A.1: ||(P_G - P_G5) b|| / ||b off G5||, G approximated by ns(eta, 40), n = 1e5\n"); print(ap, digits = 3)
write.csv(ap, file.path(HERE, "theory", "SPLINE_APPROX.csv"), row.names = FALSE)
