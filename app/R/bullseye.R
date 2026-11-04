# Bullseye plot: one nested hexagon per column; the four rings are the pSI
# thresholds (outermost 0.05, centre 0.0001); each ring is scaled to the
# number of genes specific at that threshold (radius proportional to the square
# root of the list size); ring fill is the Benjamini-Hochberg corrected p-value
# on the original tools' yellow-to-dark-red scale from 0.1 down to 0, grey when
# above 0.1.
#
# Base graphics only: the exported site runs R in the browser, where every
# package is part of the download (ggplot2 and its dependencies were 14 MB).
# This is a port of the lab's original plotting code (draw.bullseyes() and
# add.heatmap(), kept verbatim in dev/original_bullseye_functions.R):
# plotrix::draw.circle(nv = 6) becomes polygon() on six vertices, the
# gridBase legend becomes rasterImage() and text().

BULLSEYE_GREY <- "#D3D3D3"      # R's "lightgrey": fill of a ring that is not significant
BULLSEYE_CUTOFF <- 0.1          # coloured at or below this BH-corrected p-value
BULLSEYE_P_MIN <- 1e-5          # smaller values are clamped to this (the darkest red)
BULLSEYE_BAND <- 0.10           # height of the legend strip under the frame, as a fraction of the frame width

# The original ramp: 128 colours from yellow2 through orange and red to darkred,
# interpolated in Lab space with the bias towards the yellow end, shown left to
# right on the legend. The original code labelled its bins with the reversed
# palette, so BULLSEYE_PALETTE[1] (darkred) is the smallest p-value.
BULLSEYE_RAMP <- grDevices::colorRampPalette(c("yellow2", "#FF6700", "#D11500", "darkred"),
                                             bias = 0.15, space = "Lab")(128)
BULLSEYE_PALETTE <- rev(BULLSEYE_RAMP)
# 128 equal-width bins on the raw p scale between 1e-5 and 0.1, the two outer
# breaks pushed out by 0.1% of the range: exactly what the original's
# cut(c(min_pval, sigcutoff, p), breaks = 128) computed.
BULLSEYE_BREAKS <- local({
  b <- seq(BULLSEYE_P_MIN, BULLSEYE_CUTOFF, length.out = length(BULLSEYE_PALETTE) + 1)
  d <- (BULLSEYE_CUTOFF - BULLSEYE_P_MIN) / 1000
  b[1] <- b[1] - d
  b[length(b)] <- b[length(b)] + d
  b
})

# Colour for a BH-corrected p-value on the bullseye scale: dark red at (or
# below) 1e-5 through red and orange to yellow at 0.1, NA above 0.1 (drawn
# grey in the plot, left unshaded in the table). `lighten` blends towards
# white so the same hue can sit behind black text in the results table
# (0 = full plot colour, 1 = white).
bullseye_colour <- function(bh, lighten = 0) {
  out <- rep(NA_character_, length(bh))
  ok <- !is.na(bh) & bh <= BULLSEYE_CUTOFF
  if (any(ok)) {
    bin <- as.integer(cut(pmax(bh[ok], BULLSEYE_P_MIN), BULLSEYE_BREAKS))
    col <- BULLSEYE_PALETTE[bin]
    if (lighten > 0) {
      rgb <- t(grDevices::col2rgb(col)) * (1 - lighten) + 255 * lighten
      col <- grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
    }
    out[ok] <- col
  }
  out
}

# Height of the plot (frame plus legend strip) relative to its width.
bullseye_aspect <- function(lay) {
  (lay$frame[["height"]] + BULLSEYE_BAND * lay$frame[["width"]]) / lay$frame[["width"]]
}

# Vertices at 0, 60, ..., 300 degrees (a vertex points right, flat top and
# bottom), the hexagon plotrix::draw.circle(nv = 6) drew.
HEX_ANGLES <- (0:5) * pi / 3
hex_xy <- function(cx, cy, r) list(x = cx + r * cos(HEX_ANGLES), y = cy + r * sin(HEX_ANGLES))

# The white hatch lines across a significant outer ring, as the original
# draw.bullseyes() drew them: one line through the centre from the lower-left
# vertex to the upper-right one, then parallels on either side every 7 of the
# original's plot units, each ending 4 units short of the vertices. `unit` is
# the size of one of those units in frame units: the original radius was
# sqrt(n / pi) * circle_size and ours is k * sqrt(n), and the TSEA tool's own
# image (hatch period 4.5 px for a largest ring of 55 px) puts circle_size at 3.
hatch_segments <- function(cx, cy, r, unit) {
  shrink <- 4 * unit
  spacer <- 7 * unit
  c60 <- cos(pi / 3); s60 <- sin(pi / 3); s30 <- sin(pi / 6); c30 <- cos(pi / 6)
  x1 <- cx - c60 * (r - shrink); y1 <- cy - s60 * (r - shrink)
  x2 <- cx + c60 * (r - shrink); y2 <- cy + s60 * (r - shrink)
  d <- seq_len(floor((r - unit) / spacer)) * spacer      # side length = r
  n <- length(d)
  data.frame(x0 = c(x1, x1 - d * s30, x1 + d),
             y0 = c(y1, y1 + d * c30, rep(y1, n)),
             x1 = c(x2, x2 - d, x2 + d * s30),
             y1 = c(y2, rep(y2, n), y2 - d * c30))
}

# Clip line segments to a convex polygon given counter-clockwise (Cyrus-Beck).
# Segments entirely outside are dropped.
clip_segments <- function(seg, poly) {
  n <- length(poly$x)
  ex <- poly$x[c(2:n, 1)] - poly$x
  ey <- poly$y[c(2:n, 1)] - poly$y
  nx <- -ey; ny <- ex                                       # inward normals
  dx <- seg$x1 - seg$x0; dy <- seg$y1 - seg$y0
  t0 <- rep(0, nrow(seg)); t1 <- rep(1, nrow(seg)); keep <- rep(TRUE, nrow(seg))
  for (i in seq_len(n)) {
    num <- nx[i] * (seg$x0 - poly$x[i]) + ny[i] * (seg$y0 - poly$y[i])   # >= 0 inside this edge
    den <- nx[i] * dx + ny[i] * dy
    t <- -num / den
    t0 <- ifelse(den > 0, pmax(t0, t), t0)
    t1 <- ifelse(den < 0, pmin(t1, t), t1)
    keep <- keep & !(den == 0 & num < 0)
  }
  keep <- keep & t0 <= t1
  data.frame(x0 = seg$x0 + t0 * dx, y0 = seg$y0 + t0 * dy,
             x1 = seg$x0 + t1 * dx, y1 = seg$y0 + t1 * dy)[keep, , drop = FALSE]
}

# lwd for a line `w` frame units wide (lwd is in 1/96 inch; the frame's size
# in inches follows from the plot region), so that borders stay proportional
# to the plot whatever the device size.
lwd_units <- function(w) {
  usr <- graphics::par("usr"); pin <- graphics::par("pin")
  w * pin[1] / (usr[2] - usr[1]) * 96
}

# tab: long result table from sea_fisher(); lay: output of layout_nodes().
# Draws on the current device: the frame (lay$frame units, origin bottom-left)
# with the legend strip under it; the device should have the aspect ratio
# bullseye_aspect(lay), otherwise the drawing is centred with blank margins.
bullseye_plot <- function(tab, lay, base_size = 14, tree_dir = ".") {
  nodes <- lay$nodes
  fw <- lay$frame[["width"]]; fh <- lay$frame[["height"]]
  band <- BULLSEYE_BAND * fw
  op <- graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", xpd = NA,
                      ps = base_size, cex = 1, lheight = 0.85)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, fw), ylim = c(-band, fh), asp = 1)

  if (!is.null(lay$tree)) {
    f <- file.path(tree_dir, lay$tree)
    if (file.exists(f)) graphics::rasterImage(png::readPNG(f), 0, 0, fw, fh, interpolate = TRUE)
  }

  # Rings, outermost first, one column at a time in alphabetical order as the
  # original did (where two hexagons touch, the later one lies on top). The
  # TSEA tool drew the outermost (0.05) ring white with a light-grey border
  # unless significant, and hatched it when it was (`outer_faint`); the CSEA
  # tool drew it grey with a black border like the inner rings.
  lwd_border <- lwd_units(0.002 * fw)
  hatch_unit <- lay$k * sqrt(pi) / 3
  tab <- tab[tab$code %in% nodes$code, ]
  for (code in sort(unique(tab$code), method = "radix")) {
    i <- match(code, nodes$code)
    rows <- tab[tab$code == code, ]
    rows <- rows[order(-rows$threshold), ]
    for (j in seq_len(nrow(rows))) {
      r <- lay$k * sqrt(rows$n_specific[j])
      h <- hex_xy(nodes$x[i], nodes$y[i], r)
      col <- bullseye_colour(rows$bh[j])
      if (j == 1 && isTRUE(lay$outer_faint)) {
        graphics::polygon(h$x, h$y, col = if (is.na(col)) "white" else col, border = BULLSEYE_GREY, lwd = lwd_border)
        if (!is.na(col)) {
          seg <- clip_segments(hatch_segments(nodes$x[i], nodes$y[i], r, hatch_unit), h)
          graphics::segments(seg$x0, seg$y0, seg$x1, seg$y1, col = "white", lwd = lwd_border * 2 / 3)
        }
      } else {
        graphics::polygon(h$x, h$y, col = if (is.na(col)) BULLSEYE_GREY else col, border = "black", lwd = lwd_border)
      }
    }
  }

  # Labels: under the hexagon, or above it with the sublabel under it.
  t05 <- tab[tab$threshold == max(tab$threshold), ]
  r05 <- lay$k * sqrt(t05$n_specific[match(nodes$code, t05$code)])
  top <- nodes$y + r05 * sqrt(3) / 2
  bottom <- nodes$y - r05 * sqrt(3) / 2
  pad <- 0.012 * fh
  has <- !is.na(nodes$label)
  sub <- has & !is.na(nodes$sublabel)
  if (isTRUE(lay$label_above)) {
    if (any(has)) graphics::text(nodes$x[has], top[has] + pad, nodes$label[has], adj = c(0.5, 0), cex = 0.74)
    if (any(sub)) graphics::text(nodes$x[sub], bottom[sub] - pad, nodes$sublabel[sub], adj = c(0.5, 1), cex = 0.74)
  } else if (any(has)) {
    graphics::text(nodes$x[has], bottom[has] - pad, nodes$label[has], adj = c(0.5, 1), cex = 0.74)
  }
  if (!is.null(lay$headers)) {
    graphics::text(lay$headers$x, lay$headers$y, lay$headers$text, cex = 0.95)
  }

  # Legend in the strip under the frame, laid out as the original
  # add.heatmap(): the ramp from yellow (0.1) to dark red (0) on a bar half
  # the frame wide, tick marks on its upper edge, the tick labels and the title
  # "p-values" above it.
  bar_h <- 0.2 * band
  bar_top <- -0.6 * band
  graphics::rasterImage(grDevices::as.raster(matrix(BULLSEYE_RAMP, nrow = 1)),
                        0.25 * fw, bar_top - bar_h, 0.75 * fw, bar_top, interpolate = TRUE)
  at <- 0.25 * fw + (0:4) / 4 * 0.5 * fw
  graphics::segments(at, bar_top - 0.45 * bar_h, at, bar_top + 0.3 * bar_h, lwd = lwd_border * 2 / 3)
  graphics::text(at, bar_top + 0.6 * bar_h, as.character(BULLSEYE_CUTOFF * (4:0) / 4), adj = c(0.5, 0), cex = 0.7)
  graphics::text(0.5 * fw, bar_top + 1.5 * bar_h, "p-values", adj = c(0.5, 0), cex = 0.8)
  invisible(NULL)
}
