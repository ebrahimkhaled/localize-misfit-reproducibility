## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# make_compass_figure.R -- Figure 4 of the paper: the two panels that plot() draws for a verdict (../R/localize_plot.R), for the
# transported SUPPORT model M0 before and after its intercept update. Data: ../support_closure_full.rds
# (support_closure_full.R; its single-part p-values are checked against support_split_first.csv).
# Writes fig_compass.pdf in greys for the journal and fig_compass_color.pdf for the preprint and the reading copy.
source("../R/localize_plot.R")
R <- readRDS("../support_closure_full.rds")$M0
draw <- function(colour) {
  graphics::layout(matrix(1:6, 2, 3, byrow = TRUE), widths = c(.16, .34, .5))
  graphics::par(mar = c(.2, .2, .2, .2), family = "Helvetica", cex = 1)
  rows <- list(list(r = R[["verdict on A"]], t = "Frozen model,\nverdict on half A", s = "named: intercept\nnext step: update\nthe intercept"),
               list(r = R[["intercept update"]], t = "After the intercept\nupdate, retest\non half B", s = "named: cov\nnext step: revise\nthe model"))
  for (i in seq_along(rows)) {
    rw <- rows[[i]]
    graphics::plot.new(); graphics::text(.02, .62, rw$t, adj = 0, font = 2, cex = .78)
    graphics::text(.02, .3, rw$s, adj = 0, cex = .66, col = "grey25")
    .loc_compass(rw$r, colour, if (i == 1) "Misfit compass")
    .loc_lattice(rw$r, colour, if (i == 1) "Closed tests")
  }
}
grDevices::cairo_pdf("fig_compass.pdf", width = 7.2, height = 6.2, family = "Helvetica", pointsize = 11); draw(FALSE); grDevices::dev.off()
grDevices::cairo_pdf("fig_compass_color.pdf", width = 7.2, height = 6.2, family = "Helvetica", pointsize = 11); draw(TRUE); grDevices::dev.off()
if (nzchar(Sys.getenv("FIG_PNG"))) for (cl in c(FALSE, TRUE)) {
  grDevices::png(if (cl) "fig_compass_color.png" else "fig_compass.png", width = 7.2, height = 6.2, units = "in", res = 170,
                 family = "Helvetica", pointsize = 11); draw(cl); grDevices::dev.off() }
cat("fig_compass.pdf and fig_compass_color.pdf written\n")
