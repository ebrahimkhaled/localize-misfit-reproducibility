# make_app_figures.R -- the two application figures, from ../apps_update.csv and ../size_apps.csv (no typed numbers):
# fig_verdicts.pdf: verdict and size of each part for SUPPORT; fig_compass.pdf and fig_compass_color.pdf: the display of
# plot() for the transported model M0, before and after the intercept update. Run from the paper folder.
source("../R/localize_plot.R")
a <- utils::read.csv("../apps_update.csv", stringsAsFactors = FALSE)
s <- utils::read.csv("../size_apps.csv", stringsAsFactors = FALSE)
G <- c("INTERCEPT", "SLOPE", "LINK", "COV")
MAP <- data.frame(a = c("random split", "transported", "transported, half A", "half B, intercept update", "half B, logistic recalibration", "half B, flexible recalibration"),
                  s = c("random split", "transported, all chronic", "transported, half A", "half B after intercept update", "half B after logistic recalibration", "half B after flexible recalibration"),
                  lab = c("Similar patients", "Transported", "Transported, half A", "Half B: intercept update", "Half B: logistic recal.", "Half B: flexible recal."),
                  stringsAsFactors = FALSE)
panel <- function(m, x0, labels) {
  graphics::text(x0 + 2, 6.95, paste("Model", m), font = 2, cex = 0.9)
  for (j in 1:4) graphics::text(x0 + j - .5, 6.45, tolower(G[j]), font = 3, cex = 0.72)
  for (i in seq_len(nrow(MAP))) {
    y <- nrow(MAP) - i
    r <- a[a$data == "SUPPORT" & a$model == m & a$setting == MAP$a[i], ]; stopifnot(nrow(r) == 1L)
    nm <- strsplit(r$named, " + ", fixed = TRUE)[[1]]
    if (labels) graphics::text(x0 - .15, y + .5, MAP$lab[i], adj = 1, cex = 0.78)
    for (j in 1:4) {
      q <- s[s$data == "SUPPORT" & s$model == m & s$setting == MAP$s[i] & s$part == G[j], ]; stopifnot(nrow(q) == 1L)
      on <- G[j] %in% nm
      graphics::rect(x0 + j - 1, y, x0 + j, y + 1, col = if (on) "grey15" else "grey93", border = "white", lwd = 2)
      graphics::text(x0 + j - .5, y + .5, sprintf("%.1f", max(q$excess, 0)), col = if (on) "white" else "grey25", cex = 0.78, font = if (on) 2 else 1)
    }
  }
}
figv <- function() {
  graphics::par(mar = c(.2, .2, .2, .2), family = "sans")
  graphics::plot.new(); graphics::plot.window(xlim = c(-3.3, 9.3), ylim = c(0, 7.2))
  panel("M0", 0, TRUE); panel("M2", 5, FALSE)
}
grDevices::cairo_pdf("fig_verdicts.pdf", width = 7.2, height = 3.0, family = "sans", pointsize = 11); figv(); grDevices::dev.off()
grDevices::png("fig_verdicts.png", width = 7.2, height = 3.0, units = "in", res = 170, pointsize = 11); figv(); grDevices::dev.off()

obj <- function(setting) {
  r <- a[a$data == "SUPPORT" & a$model == "M0" & a$setting == setting, ]
  p <- stats::setNames(as.numeric(r[G]), G)
  sets <- unlist(lapply(1:4, function(k) apply(utils::combn(G, k), 2, paste, collapse = "+")))
  full <- stats::setNames(rep(NA_real_, length(sets)), sets); full[G] <- p
  structure(list(named = G[stats::p.adjust(p, "holm") <= 0.05], intersection = full, adjusted = stats::p.adjust(p, "holm"),
                 naming = "holm", alpha = 0.05, n = r$n), class = "gof_localize")
}
draw <- function(colour) {
  graphics::layout(matrix(1:6, 2, 3, byrow = TRUE), widths = c(.16, .34, .5))
  graphics::par(mar = c(.2, .2, .2, .2), family = "sans", cex = 1)
  rows <- list(list(r = obj("transported, half A"), t = "Frozen model,\nhalf A", s = "named: intercept, cov\nnext step: revise\nthe model"),
               list(r = obj("half B, intercept update"), t = "After the intercept\nupdate, retest\non half B", s = "named: cov\nthe update did not\nrepair this part"))
  stopifnot(identical(rows[[1]]$r$named, c("INTERCEPT", "COV")), identical(rows[[2]]$r$named, "COV"))
  for (i in seq_along(rows)) {
    rw <- rows[[i]]
    graphics::plot.new(); graphics::text(.02, .62, rw$t, adj = 0, font = 2, cex = .78)
    graphics::text(.02, .3, rw$s, adj = 0, cex = .66, col = "grey25")
    .loc_compass(rw$r, colour, if (i == 1) "Misfit compass")
    .loc_lattice(rw$r, colour, if (i == 1) "Tests of each set of parts (Holm)")
  }
}
grDevices::cairo_pdf("fig_compass.pdf", width = 7.2, height = 6.2, family = "sans", pointsize = 11); draw(FALSE); grDevices::dev.off()
grDevices::cairo_pdf("fig_compass_color.pdf", width = 7.2, height = 6.2, family = "sans", pointsize = 11); draw(TRUE); grDevices::dev.off()
grDevices::png("fig_compass.png", width = 7.2, height = 6.2, units = "in", res = 170, pointsize = 11); draw(FALSE); grDevices::dev.off()
cat("figures written\n")
