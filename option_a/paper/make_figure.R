## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# make_figure.R -- Figure 1 of the paper: the rate of naming a group with no misfit against n, for departures that
# shrink like 1/sqrt(n) (left) and departures of fixed size (right); external validation, Monte Carlo reference.
# Data: ../theory/LOCAL.csv. Black on white, line type and symbol by departure family, so it reads in greyscale.
L <- utils::read.csv("../theory/LOCAL.csv", na.strings = "")
L <- L[L$setting == "external", ]
fam <- data.frame(key = c("slope_pure", "link_pure", "cov_pure", "slope_mod", "link_mod", "cov_mod"),
                  lab = c("slope, large", "link, large", "covariate, large", "slope, moderate", "link, moderate", "covariate, moderate"),
                  lty = c(1, 2, 3, 1, 2, 3), pch = c(16, 17, 15, 1, 2, 0), stringsAsFactors = FALSE)
grDevices::cairo_pdf("fig_local.pdf", width = 6.6, height = 3.1, family = "Helvetica")
op <- graphics::par(mfrow = c(1, 2), mar = c(3.4, 3.6, 1.8, 0.6), mgp = c(2.2, 0.6, 0), tcl = -0.3, las = 1, cex = 1,
                    cex.axis = 0.75, cex.lab = 0.8, cex.main = 0.8, font.main = 1)
panel <- function(track, keys, main) {
  graphics::plot(NA, xlim = c(450, 36000), ylim = c(0, 1), log = "x", xaxt = "n", xlab = "Sample size n",
                 ylab = "Rate of naming a group with no misfit", main = main, bty = "l")
  graphics::axis(1, at = c(500, 1000, 2000, 4000, 8000, 16000, 32000), labels = c("500", "1k", "2k", "4k", "8k", "16k", "32k"))
  graphics::abline(h = 0.05, col = "grey55", lwd = 0.8)
  graphics::text(470, 0.075, "0.05", adj = c(0, 0), cex = 0.65, col = "grey35")
  for (k in keys) {
    f <- fam[fam$key == k, ]; x <- L[L$fam == k & L$track == track, ]; x <- x[order(x$n), ]
    if (!nrow(x)) next
    graphics::lines(x$n, x$wrong, lty = f$lty, lwd = 1.1)
    graphics::points(x$n, x$wrong, pch = f$pch, cex = 0.8)
  }
}
panel("local", fam$key, "Departure shrinking like 1/sqrt(n)")
panel("fixed", fam$key[1:3], "Departure of fixed size")
graphics::legend(2600, 0.62, legend = fam$lab, lty = fam$lty, pch = fam$pch, cex = 0.68, bty = "n", lwd = 1.1, seg.len = 2.2)
graphics::par(op); grDevices::dev.off()
cat("fig_local.pdf written\n")
