## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# make_support_figure.R -- Figure 2: calibration of the transported SUPPORT models on half B of the chronic patients,
# frozen and after each update fitted on half A (the split-first update-and-retest of support_split_first.R; same
# data, models, split and seeds). Each curve is a loess smooth of the outcome on the predicted risk; the legend gives
# the parts named on B (support_split_first.csv). Black on white; line type by update.
suppressMessages(library(splines))
s <- utils::read.csv(file.path(.ROOT, "data/support2.csv"),
                     na.strings = c("NaN", "NA", ""))
d <- stats::na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                               crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                               acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
F0 <- y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + male
F2 <- y ~ age + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) + ns(wblc, 3) + num.co + male
v <- d[d$acute == 0, ]
set.seed(20261011L); A <- sort(sample(nrow(v), floor(nrow(v) / 2))); B <- setdiff(seq_len(nrow(v)), A)
sf <- utils::read.csv("../support_split_first.csv", stringsAsFactors = FALSE)
lab_named <- function(m, rung) {
  x <- sf$named[sf$model == m & sf$rung == rung]
  if (!length(x) || x == "none") "nothing named" else tolower(gsub(" \\+ ", ", ", x))
}
curves <- function(m) {
  fit <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = d[d$acute == 1, ], family = binomial()))
  p <- as.numeric(predict(fit, newdata = v, type = "response")); e <- qlogis(p)
  f1 <- glm(v$y[A] ~ e[A], family = binomial())
  f3 <- suppressWarnings(glm(if (m == "M0") F0 else F2, data = v[A, ], family = binomial()))
  presc <- if (m == "M0") { fi <- glm(v$y[A] ~ offset(e[A]), family = binomial()); plogis(coef(fi)[1] + e[B]) } else {
    ff <- glm(v$y[A] ~ ns(e[A], df = 3), family = binomial()); as.numeric(plogis(cbind(1, predict(ns(e[A], df = 3), e[B])) %*% coef(ff))) }
  list(list(p = p[B], lab = "frozen model", lty = 1, lwd = 1.6),
       list(p = presc, lab = paste0(if (m == "M0") "intercept update" else "flexible recalibration", " (prescribed): ",
                                    lab_named(m, if (m == "M0") "intercept update" else "flexible recalibration")), lty = 2, lwd = 1.3),
       list(p = plogis(coef(f1)[1] + coef(f1)[2] * e[B]), lab = paste0("logistic recalibration: ", lab_named(m, "logistic recalibration")),
            lty = 3, lwd = 1.5),
       list(p = as.numeric(predict(f3, newdata = v[B, ], type = "response")), lab = paste0("re-estimation: ", lab_named(m, "model revision")),
            lty = 4, lwd = 1.3))
}
grDevices::cairo_pdf("fig_support.pdf", width = 6.6, height = 3.4, family = "Helvetica")
op <- graphics::par(mfrow = c(1, 2), mar = c(3.4, 3.6, 1.8, 0.6), mgp = c(2.2, 0.6, 0), tcl = -0.3, las = 1, cex = 1,
                    cex.axis = 0.75, cex.lab = 0.8, cex.main = 0.8, font.main = 1)
for (m in c("M0", "M2")) {
  graphics::plot(NA, xlim = c(0, 0.65), ylim = c(0, 0.65), xlab = "Predicted risk", ylab = "Observed risk (loess)",
                 main = if (m == "M0") "Linear model M0" else "Spline model M2", bty = "l")
  graphics::abline(0, 1, col = "grey60", lwd = 0.8)
  cs <- curves(m)
  for (k in cs) {
    o <- order(k$p); lo <- stats::loess(v$y[B][o] ~ k$p[o], span = 0.75, degree = 1)
    xs <- seq(stats::quantile(k$p, .02), stats::quantile(k$p, .98), length.out = 100)
    graphics::lines(xs, pmin(pmax(stats::predict(lo, xs), 0), 1), lty = k$lty, lwd = k$lwd)
  }
  graphics::legend("topleft", legend = vapply(cs, `[[`, "", "lab"), lty = vapply(cs, `[[`, 0, "lty"),
                   lwd = vapply(cs, `[[`, 0, "lwd"), cex = 0.56, bty = "n", seg.len = 2.4)
}
graphics::par(op); grDevices::dev.off()
cat("fig_support.pdf written\n")
