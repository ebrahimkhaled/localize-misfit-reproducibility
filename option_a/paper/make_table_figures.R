## archive root: the environment variable LOCALIZE_ARCHIVE_ROOT, or the working directory (see README)
.ROOT <- normalizePath(Sys.getenv("LOCALIZE_ARCHIVE_ROOT", "."), winslash = "/")
Sys.setenv(LOCALIZE_ARCHIVE_ROOT = .ROOT)
# make_table_figures.R -- Figures 1 and 3 of the paper: the simulation (n = 1,000) as a map of naming rates, and the
# SUPPORT verdicts as a ladder of p-values. Data come from ../SIM_GROUPS.csv, ../support_v2.csv and
# ../support_split_first.csv (no typed numbers). Black, white and greys only, so they read in print.
OA <- ".."   # run from the paper folder
PARTS <- c("INTERCEPT", "SLOPE", "LINK", "COV")
shade <- function(x) grDevices::grey(1 - 0.85 * pmin(pmax(x, 0), 1))   # 0 = white, 1 = near black
ink <- function(x) ifelse(x > 0.55, "white", "black")
both <- function(name, w, h, f) {                  # pointsize 13: the smallest label prints at about 7 pt
  grDevices::cairo_pdf(paste0(name, ".pdf"), width = w, height = h, family = "Helvetica", pointsize = 13); f(); grDevices::dev.off()
  if (nzchar(Sys.getenv("FIG_PNG")))
    { grDevices::png(paste0(name, ".png"), width = w, height = h, units = "in", res = 200, family = "Helvetica", pointsize = 13); f(); grDevices::dev.off() }
}

## ---- A. Table 2 as a map: share of data sets naming each part; the part each departure needs is framed
s <- utils::read.csv(file.path(OA, "SIM_GROUPS.csv"), na.strings = ""); s <- s[s$n == 1000, ]
FAM <- data.frame(fam = c("null", "intercept", "slope", "link", "ushape", "thresh", "inter"),
                  C = c(0, .4, .6, .25, .8, 4, 1),
                  lab = c("no misfit", "intercept shift", "slope 0.6", "curved link", "U-shape in x1", "threshold in x2", "interaction x1 x3"),
                  need = c(NA, "INTERCEPT", "SLOPE", "LINK", "COV", "COV", "COV"), stringsAsFactors = FALSE)
get <- function(setting, f, C, g) { r <- s[s$setting == setting & s$fam == f & abs(s$C - C) < 1e-9, ]; if (nrow(r)) r[[g]][1] else NA }
figA <- function() {
  graphics::par(mar = c(0.4, 0.4, 0.4, 0.4), family = "Helvetica")
  graphics::plot.new(); graphics::plot.window(xlim = c(-2.9, 10.6), ylim = c(-2.5, 8.6), asp = NA)
  nr <- nrow(FAM)
  blocks <- list(list(x0 = 0, set = "external", cols = PARTS, title = "External validation"),
                 list(x0 = 6.7, set = "insample", cols = c("LINK", "COV"), title = "In-sample"))
  for (b in blocks) {
    graphics::text(b$x0 + length(b$cols) / 2, 8.35, b$title, font = 2, cex = 0.8)
    for (j in seq_along(b$cols)) graphics::text(b$x0 + j - 0.5, 7.75, tolower(b$cols[j]), cex = 0.6, font = 3)
    graphics::text(b$x0 + length(b$cols) + 0.85, 7.75, "right action", cex = 0.62, font = 3)
    for (i in seq_len(nr)) {
      y <- nr - i + 0.5
      for (j in seq_along(b$cols)) {
        v <- get(b$set, FAM$fam[i], FAM$C[i], b$cols[j])
        if (is.na(v)) { graphics::rect(b$x0 + j - 1, y - .5, b$x0 + j, y + .5, col = "white", border = "grey80")
                        graphics::text(b$x0 + j - .5, y, "not visible", cex = .5, col = "grey45"); next }
        graphics::rect(b$x0 + j - 1, y - .5, b$x0 + j, y + .5, col = shade(v), border = "white", lwd = 1.5)
        graphics::text(b$x0 + j - .5, y, sprintf("%.2f", v), cex = 0.62, col = ink(v))
        if (!is.na(FAM$need[i]) && FAM$need[i] == b$cols[j])
          graphics::rect(b$x0 + j - 0.96, y - .46, b$x0 + j - .04, y + .46, border = "black", lwd = 2.2)
      }
      ra <- get(b$set, FAM$fam[i], FAM$C[i], "right_action")
      xb <- b$x0 + length(b$cols) + 0.2
      graphics::rect(xb, y - .28, xb + 1.4, y + .28, col = "grey93", border = NA)
      if (!is.na(ra) && !(b$set == "insample" && FAM$need[i] %in% c("INTERCEPT", "SLOPE"))) {
        graphics::rect(xb, y - .28, xb + 1.4 * ra, y + .28, col = "grey25", border = NA)
        graphics::text(xb + 1.4 * ra + 0.06, y, sprintf("%.2f", ra), adj = 0, cex = .55)
      }
    }
  }
  for (i in seq_len(nr)) graphics::text(-0.15, nr - i + 0.5, FAM$lab[i], adj = 1, cex = 0.68)
  graphics::text(-0.15, 7.75, "departure", adj = 1, cex = 0.62, font = 3)
  graphics::rect(0, -0.95, .45, -0.55, border = "black", lwd = 2.2)
  graphics::text(0.6, -0.75, "framed: the part the departure needs (dark framed cell = right verdict)", adj = 0, cex = .6)
  graphics::text(0, -1.5, "share of data sets naming the part:", adj = 0, cex = .6)
  graphics::text(4.05, -1.5, "0", adj = 1, cex = .6)
  for (k in 0:4) graphics::rect(4.15 + k * .4, -1.7, 4.55 + k * .4, -1.3, col = shade(k / 4), border = "grey70")
  graphics::text(6.25, -1.5, "1", adj = 0, cex = .6)
  graphics::text(0, -2.2, "right action: the highest part named is the part the departure needs (no misfit: nothing named)", adj = 0, cex = .6)
}
both("fig_sim", 7.2, 4.7, figA)

## ---- B. Table 3 as a ladder of p-values along the SUPPORT story (M0 and M2 side by side)
v <- utils::read.csv(file.path(OA, "support_v2.csv"), na.strings = "NA", stringsAsFactors = FALSE)
sf <- utils::read.csv(file.path(OA, "support_split_first.csv"), na.strings = "NA", stringsAsFactors = FALSE)
sf$analysis <- ifelse(sf$stage == "verdict on A", "split: verdict on A", paste0("split: ", sf$rung))
v <- rbind(v[, c("analysis", "model", "named", PARTS)], sf[, c("analysis", "model", "named", PARTS)])
ROWS <- data.frame(a = c("in-sample", "in-sample, de-aliased", "external random, montecarlo", "external transported, montecarlo",
                         "split: verdict on A", "split: intercept update", "split: logistic recalibration",
                         "split: flexible recalibration", "split: model revision"),
                   lab = c("all patients", "de-aliased link group", "random half", "acute to chronic", "verdict on half A",
                           "after intercept update", "after logistic recalibration", "after flexible recalibration", "after re-estimation"),
                   grp = c(1, 1, 2, 3, 4, 5, 5, 5, 5), stringsAsFactors = FALSE)
GRP <- c("In-sample", "Random split", "Transported", "Split first", "Retest on half B")
figB <- function() {
  graphics::par(mar = c(0.4, 0.4, 0.4, 0.4), family = "Helvetica")
  graphics::plot.new(); graphics::plot.window(xlim = c(-5.2, 10.2), ylim = c(-1.3, 12.2))
  nr <- nrow(ROWS); gap <- cumsum(c(0, diff(ROWS$grp) != 0)) * 0.35
  yy <- (nr - seq_len(nr)) + 0.5 + (max(gap) - gap)
  for (m in 1:2) {
    x0 <- c(0, 5.6)[m]; mod <- c("M0", "M2")[m]
    graphics::text(x0 + 2, max(yy) + 1.45, c("M0: linear model", "M2: spline model")[m], font = 2, cex = 0.8)
    for (j in 1:4) graphics::text(x0 + j - .5, max(yy) + 0.8, tolower(PARTS[j]), cex = .6, font = 3)
    for (i in seq_len(nr)) {
      r <- v[v$analysis == ROWS$a[i] & v$model == mod, ]
      named <- if (nrow(r)) strsplit(r$named[1], " + ", fixed = TRUE)[[1]] else character(0)
      for (j in 1:4) {
        p <- if (nrow(r)) r[[PARTS[j]]][1] else NA
        xl <- x0 + j - 1; y <- yy[i]
        if (is.na(p)) { graphics::rect(xl, y - .45, xl + 1, y + .45, col = "white", border = "grey85")
                        graphics::text(xl + .5, y, "not visible", cex = .5, col = "grey45"); next }
        isn <- PARTS[j] %in% named
        graphics::rect(xl, y - .45, xl + 1, y + .45, col = if (isn) "grey15" else shade(min(-log10(p), 3) / 3 * .45),
                       border = "white", lwd = 1.5)
        graphics::text(xl + .5, y, if (p < .001) "<.001" else sub("^0", "", sprintf("%.3f", p)), cex = .55,
                       col = if (isn) "white" else "black", font = if (isn) 2 else 1)
      }
    }
  }
  for (i in seq_len(nr)) graphics::text(-0.15, yy[i], ROWS$lab[i], adj = 1, cex = .64)
  for (g in unique(ROWS$grp)) {
    ys <- yy[ROWS$grp == g]
    graphics::segments(-3.55, min(ys) - .4, -3.55, max(ys) + .4, lwd = 1.2)
    graphics::text(-3.7, mean(ys), GRP[g], adj = 1, cex = .64, font = 2)
  }
  graphics::rect(0, -1.15, .5, -0.75, col = "grey15", border = NA); graphics::text(.6, -.95, "part named", adj = 0, cex = .58)
  graphics::rect(2.6, -1.15, 3.1, -0.75, col = shade(.45), border = NA)
  graphics::text(3.2, -.95, "small p, not named (closure needs every larger set to reject)", adj = 0, cex = .58)
}
both("fig_pvalues", 7.2, 4.9, figB)

cat("fig_sim.pdf and fig_pvalues.pdf written\n")
