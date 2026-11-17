# One large figure: one time axis, then the thesis, Gaussian and flat-top truths at the same 7 times.
# Exploratory: the two smooth truths are alternatives, not the truth used in the thesis.
# Run from the project root, the folder that holds the .here file:
#   Rscript code/extra/ground_truth_comparison/ground_truth_comparison.R
source("code/00_config.R")

amplitude <- 1.5
onset <- 160
offset <- 340
mid <- (onset + offset) / 2
radius_taper <- 0.70     # end of the Gaussian taper
sigma_space <- 0.50
r_thesis <- sigma_space * sqrt(2 * log(2))
lim <- 1.1               # half-width of the scalp grid
col_in <- "#D62728"
col_out <- "#1F77B4"
col_spline <- "#247BA0"
col_thesis <- "#E8590C"
palette <- hcl.colors(100, "viridis")

# --- temporal functions ---------------------------------------------------
bump_bspline <- function(z) {
    a <- abs(z)
    out <- numeric(length(a))
    i1 <- a < 1
    out[i1] <- (4 - 6 * a[i1]^2 + 3 * a[i1]^3) / 4
    i2 <- a >= 1 & a < 2
    out[i2] <- (2 - a[i2])^3 / 4
    out
}
temporal_spline <- function(t) bump_bspline(4 * (t - mid) / (offset - onset))

source(project_path("code", "R", "superpopulation_truth.R"))
wf <- population_effect_waveform(
    bump_amplitude    = amplitude,
    onset_support_ms  = seq(140, 180, by = 2),
    downsample_factor = 1L,
    template_path     = project_path("resources", "brms_meeg", "code", "erp_template.R"))
temporal_thesis <- function(t) {
    fit <- approx(
        x      = wf$time * 1000,                  # template time, converted from s to ms
        y      = wf$effect_waveform / amplitude,  # waveform scaled to peak units
        xout   = t,                               # times at which we want the value
        yleft  = 0,                               # value before the first template time
        yright = 0)                               # value after the last template time
    fit$y
}

# --- spatial functions ----------------------------------------------------
gauss <- function(d) exp(-d^2 / (2 * sigma_space^2))
thesis_cut <- function(d) gauss(d) * (d <= r_thesis)
gaussian_taper <- function(d) {
    u <- pmin(1, pmax(0, (d - taper_start) / (radius_taper - taper_start)))
    gauss(d) * (1 - u)^3 * (1 + 3 * u + 6 * u^2)
}
near_flat_top <- function(d) pmax(0, 1 - (d / radius_flat)^5)^3

# --- electrodes -----------------------------------------------------------
elec <- read.csv(project_path("resources", "electrodes", "electrode_layout_61ch.csv"))
elec$d <- sqrt(elec$xproj^2 + elec$yproj^2)
radius_flat <- sort(elec$d)[10]   # near-flat-top support radius: the 10th smallest distance (C3 and C4)
nine <- elec[elec$d <= radius_taper, ]
taper_start <- max(nine$d)   # start of the taper: the distance of the outer ring (FC1, FC2, CP1, CP2)
others <- elec[elec$d > radius_taper & abs(elec$xproj) <= lim & abs(elec$yproj) <= lim, ]
stopifnot(nrow(nine) == 9, sum(elec$d < radius_flat) == 9)

# --- common time points (taken from the spline) ---------------------------
at_level <- function(p) uniroot(function(t) temporal_spline(t) - p, c(onset, mid))$root
t25 <- at_level(0.25)
t75 <- at_level(0.75)
snap_times <- c(onset, t25, t75, mid, 2 * mid - t75, 2 * mid - t25, offset)
snap_labels <- c("Onset", "25% of peak", "75% of peak", "Peak", "75% of peak", "25% of peak", "Offset")

# --- plotting helpers -----------------------------------------------------
xy <- seq(-lim, lim, length.out = 301)
dgrid <- outer(xy, xy, function(x, y) sqrt(x^2 + y^2))

# Top panel: both time courses at Cz with the seven time points
draw_timecourse <- function() {
    par(mar = c(4.5, 5, 3, 2), las = 1)
    plot(NA, xlim = c(100, 400), ylim = c(0, 1.7), xlab = "Time (ms)", ylab = "Effect at Cz (µV)",
         main = "Time courses at Cz and the seven common time points")
    abline(v = c(onset, offset), lty = 3, col = "#8795A1")
    abline(v = snap_times, lty = 2, col = "grey55")
    tt <- seq(0, 500, by = 2)
    lines(tt, amplitude * temporal_spline(tt), lwd = 3.5, col = col_spline)
    lines(tt, amplitude * temporal_thesis(tt), lwd = 3.5, col = col_thesis)
    points(snap_times, amplitude * temporal_spline(snap_times), pch = 21, bg = col_spline, cex = 1.5)
    points(snap_times, amplitude * temporal_thesis(snap_times), pch = 21, bg = col_thesis, cex = 1.5)
    text(snap_times, 1.68, seq_along(snap_times), font = 2, cex = 1.1)
    legend("topleft", bty = "n", lwd = 3.5, col = c(col_thesis, col_spline), cex = 1.0,
           legend = c("Thesis: mean of shifted ERP templates (onsets 140 to 180 ms)",
                      "Smooth B-spline (used with the two new spatial profiles)"))
}

# A boxed value at the left or right edge of the panel, joined to (cx, cy) by a line
draw_value_label <- function(value, cx, cy, side, col) {
    u <- par("usr")
    lab <- sprintf("%.2f", value)
    lx <- if (side == "L") u[1] + 0.07 * diff(u[1:2]) else u[2] - 0.07 * diff(u[1:2])
    segments(cx, cy, lx, cy, col = col, lwd = 1)
    rect(lx - strwidth(lab) * 0.6, cy - strheight(lab) * 0.9,
         lx + strwidth(lab) * 0.6, cy + strheight(lab) * 0.9,
         col = "white", border = col)
    text(lx, cy, lab, col = col, font = 2, cex = 0.8)
}

# One 3D scalp panel at time tm
draw_panel <- function(profile, temporal, tm, k, row_title, labels) {
    h <- amplitude * temporal(tm)
    z <- h * matrix(profile(dgrid), nrow = length(xy))
    nr <- nrow(z)
    nc <- ncol(z)
    zf <- (z[-nr, -nc] + z[-1, -nc] + z[-nr, -1] + z[-1, -1]) / 4
    col <- palette[pmin(100, pmax(1, round(zf / amplitude * 99) + 1))]
    col[zf < 0.03] <- "#FAFBFC"
    title <- if (nzchar(labels[k])) sprintf("%d. %s", k, labels[k]) else sprintf("%d.", k)
    view <- persp(xy, xy, z, zlim = c(0, amplitude), theta = 35, phi = 30, expand = 0.7,
                  col = col, border = NA, shade = 0.2, box = FALSE, axes = FALSE,
                  main = sprintf("%s\n%.0f ms (%.2f µV)", title, tm, max(z)), cex.main = 1.0)
    if (k == 1) mtext(row_title, side = 3, line = 4.6, adj = 0, font = 2, cex = 1.2)

    lines(trans3d(c(-lim, lim, lim, -lim, -lim), c(-lim, -lim, lim, lim, -lim), 0, view),
          col = "grey60", lwd = 1.0)
    points(trans3d(others$xproj, others$yproj, 0, view), pch = 19, cex = 0.35, col = "grey35")

    height <- h * profile(nine$d)
    pts <- trans3d(nine$xproj, nine$yproj, height + 0.02, view)
    for (grp in list(list(c("FCZ", "C2", "CPZ", "C1"), col_in, "L"),
                     list(c("FC2", "CP2", "CP1", "FC1"), col_out, "R"))) {
        idx <- match(grp[[1]], nine$channel)
        polygon(pts$x[idx], pts$y[idx], col = adjustcolor(grp[[2]], 0.35),
                border = grp[[2]], lwd = 1.2)
        value <- mean(height[idx])
        if (value >= 0.005)
            draw_value_label(value, mean(pts$x[idx]), mean(pts$y[idx]), grp[[3]], grp[[2]])
    }
    points(pts, pch = 21, bg = "white", col = "black", cex = 0.9)
    if (tm == mid) {
        dx <- 0.004 * diff(par("usr")[1:2])
        for (o in list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1)))
            text(pts$x + o[1] * dx, pts$y + o[2] * dx, nine$channel, pos = 3,
                 cex = 0.6, font = 2, offset = 0.4, col = "white")
        text(pts$x, pts$y, nine$channel, pos = 3, cex = 0.6, font = 2, offset = 0.4)
    }
}

# --- the three rows of scalp panels ----------------------------------------
rows <- list(
    list(thesis_cut, temporal_thesis,
         sprintf("Thesis ground truth: Gaussian cut at the half-maximum radius (r = %.3f)", r_thesis),
         c("Window start", "", "", "Peak", "", "", "Window end")),
    list(gaussian_taper, temporal_spline,
         sprintf("Gaussian with a smooth taper from %.3f to %.2f", taper_start, radius_taper),
         snap_labels),
    list(near_flat_top, temporal_spline,
         sprintf("Near-flat-top bump (radius %.4f, the 10th smallest distance)", radius_flat),
         snap_labels))

# --- main: draw the figure -------------------------------------------------
out_dir <- fig_dir(project_path("code", "extra", "ground_truth_comparison"))
png(file.path(out_dir, "ground_truth_comparison.png"), width = 5600, height = 3900, res = 240, bg = "white")
layout(rbind(rep(1, 7), matrix(2:22, nrow = 3, byrow = TRUE)), heights = c(0.8, 1, 1, 1))
par(oma = c(2, 0, 3, 0))
draw_timecourse()
par(mar = c(0.5, 0.5, 6.5, 0.5))
for (rw in rows)
    for (k in 1:7) draw_panel(rw[[1]], rw[[2]], snap_times[k], k, rw[[3]], rw[[4]])
mtext("Ground truth: the thesis and two smooth alternatives at the same seven time points",
      outer = TRUE, side = 3, font = 2, cex = 1.3)
mtext("Dots: the nine active electrodes. Red plane joins the inner four, blue plane the outer four; labels give the mean effect of each plane in µV.",
      outer = TRUE, side = 1, cex = 0.95)
dev.off()
