## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# localize_plot.R -- the display of a localization verdict, for a result of localize_external() or localize_insample().
#   plot_localize(r)      the misfit compass (left) and the closed-testing lattice (right), side by side
#   localize_compass(r)   one spoke per part; distance from the centre = strength of evidence, -log10 p capped at
#                         p = 0.001; dashed ring = p 0.05; the parts named are drawn thick
#   localize_lattice(r)   every set of the parts tested, filled when its test rejects; a part is named when every
#                         set that contains it rejects (closure)
# colour = FALSE draws in greys only (journals that ask for black on white). Parts that cannot be tested in the
# setting (intercept and slope in-sample) are shown as grey spokes and left out of the lattice.

.LOC_PARTS <- c("INTERCEPT", "SLOPE", "LINK", "COV")
.loc_pal <- function(colour) if (colour) c(INTERCEPT = "#0072B2", SLOPE = "#009E73", LINK = "#9E3D8C", COV = "#D55E00") else
  c(INTERCEPT = "grey10", SLOPE = "grey10", LINK = "grey10", COV = "grey10")
.loc_lp <- function(p) pmin(-log10(pmax(p, 1e-3)), 3) / 3
.loc_pfmt <- function(p) ifelse(p < .001, "p < .001", paste0("p = ", sub("^0", "", sprintf("%.3f", p))))
.loc_split <- function(s) strsplit(s, "+", fixed = TRUE)[[1]]

localize_compass <- function(r, colour = TRUE, head = NULL, alpha = 0.05) {
  PC <- .loc_pal(colour); G <- .LOC_PARTS
  ang <- c(90, 0, 270, 180) * pi / 180                      # intercept top, slope right, link bottom, cov left
  graphics::plot(NA, xlim = c(-1.55, 1.55), ylim = c(-1.45, 1.6), asp = 1, axes = FALSE, xlab = "", ylab = "")
  if (!is.null(head)) graphics::text(0, 1.6, head, font = 2, cex = .8)
  tt <- seq(0, 2 * pi, length.out = 240)
  for (q in c(.01, .001)) graphics::lines(.loc_lp(q) * cos(tt), .loc_lp(q) * sin(tt), col = "grey80", lwd = .7)
  graphics::lines(.loc_lp(alpha) * cos(tt), .loc_lp(alpha) * sin(tt), lty = 2, lwd = 1)
  for (q in c(alpha, .01, .001)) graphics::text(.loc_lp(q) * cos(pi / 4) + .04, .loc_lp(q) * sin(pi / 4) + .04,
                                                sub("^0", "", format(q)), cex = .55, col = "grey35", adj = 0)
  tested <- G %in% names(r$intersection)
  graphics::segments(0, 0, 1.04 * cos(ang), 1.04 * sin(ang), col = ifelse(tested, "grey75", "grey88"), lty = ifelse(tested, 1, 3))
  p <- ifelse(tested, as.numeric(r$intersection[G]), NA); rad <- ifelse(tested, .loc_lp(p), 0)
  x <- rad * cos(ang); y <- rad * sin(ang); nm <- G %in% r$named
  if (sum(tested) >= 3) graphics::polygon(x[tested], y[tested], col = if (colour) "#F2F2F2" else "grey93", border = "grey45")
  else graphics::segments(0, 0, x[tested], y[tested], col = "grey45")
  graphics::segments(0, 0, x[nm], y[nm], lwd = 5, col = PC[G][nm], lend = 1)
  graphics::points(x[tested], y[tested], pch = 21, bg = ifelse(nm, PC[G], "white")[tested],
                   col = ifelse(nm, PC[G], "grey30")[tested], cex = 1.2)
  graphics::text(1.27 * cos(ang), 1.27 * sin(ang) + c(0.03, 0, -0.03, 0), tolower(G),
                 col = ifelse(tested, PC[G], "grey70"), font = ifelse(nm, 2, 1), cex = .8)
  graphics::text(1.27 * cos(ang), 1.27 * sin(ang) + c(-.15, -.15, -.17, -.15), ifelse(tested, .loc_pfmt(p), "not tested"),
                 cex = .55, col = ifelse(tested, "grey25", "grey65"))
  invisible(r)
}

localize_lattice <- function(r, colour = TRUE, head = NULL, alpha = 0.05) {
  PC <- .loc_pal(colour); REJ <- if (colour) "#E8EEF5" else "grey85"
  sets <- names(r$intersection); lev <- vapply(sets, function(s) length(.loc_split(s)), 1L); K <- max(lev)
  hh <- function(k) .13 * k + .17; hw <- .64                             # box half height grows with the names
  Y <- .55; for (k in seq_len(K - 1)) Y[k + 1] <- Y[k] + hh(k) + hh(k + 1) + .67
  xy <- do.call(rbind, lapply(seq_len(K), function(k) { s <- sets[lev == k]; nk <- length(s); W <- min(1.55, 8.6 / max(nk, 1))
    data.frame(s = s, k = k, x = 5 + (seq_len(nk) - (nk + 1) / 2) * W, y = Y[k]) }))
  top <- max(xy$y + hh(xy$k)); H <- 6.45                     # the height of the full four-part lattice, so that
  sh <- max(0, (H - top - .45) / 2); xy$y <- xy$y + sh; top <- top + sh   # boxes keep their size with fewer parts
  graphics::plot(NA, xlim = c(0.4, 9.6), ylim = c(0, max(H, top + .45)), axes = FALSE, xlab = "", ylab = "")
  if (!is.null(head)) graphics::text(1.0, top + .35, head, font = 2, cex = .8, adj = 0)
  for (i in seq_len(nrow(xy))) for (j in seq_len(nrow(xy)))
    if (xy$k[j] == xy$k[i] + 1 && all(.loc_split(xy$s[i]) %in% .loc_split(xy$s[j])))
      graphics::segments(xy$x[i], xy$y[i] + hh(xy$k[i]), xy$x[j], xy$y[j] - hh(xy$k[j]), col = "grey78", lwd = .8)
  for (i in seq_len(nrow(xy))) {
    s <- xy$s[i]; k <- xy$k[i]; rej <- r$intersection[[s]] <= alpha; named <- k == 1 && s %in% r$named
    graphics::rect(xy$x[i] - hw, xy$y[i] - hh(k), xy$x[i] + hw, xy$y[i] + hh(k),
                   col = if (named) PC[[s]] else if (rej) REJ else "white",
                   border = if (rej) "grey20" else "grey65", lwd = if (named) 3 else if (rej) 1.2 else .8)
    nm <- .loc_split(s)
    graphics::text(xy$x[i], xy$y[i] + hh(k) - .17 - (seq_along(nm) - 1) * .26, tolower(nm), cex = .56,
                   col = if (named) "white" else PC[nm], font = if (named) 2 else 1)
    graphics::text(xy$x[i], xy$y[i] - hh(k) + .13, .loc_pfmt(r$intersection[[s]]), cex = .46,
                   col = if (named) "white" else "grey20", font = if (rej) 2 else 1)
  }
  invisible(r)
}

plot_localize <- function(r, colour = TRUE, alpha = 0.05) {
  op <- graphics::par(mfrow = c(1, 2), mar = c(.2, .2, .2, .2)); on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, 1), widths = c(.4, .6))
  localize_compass(r, colour, "Misfit compass", alpha)
  localize_lattice(r, colour, "Closed testing: every set of parts", alpha)
  invisible(r)
}
