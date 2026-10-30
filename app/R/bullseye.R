# Bullseye plot: one nested hexagon per column; the four rings are the pSI
# thresholds (outermost 0.05, centre 0.0001); each ring is scaled to the
# number of genes specific at that threshold (radius proportional to the square
# root of the list size); ring fill is the Benjamini-Hochberg corrected p-value
# on a yellow-to-red scale from 0.1 down to 0, grey when above 0.1.

BULLSEYE_GREY <- "#D3D3D3"
BULLSEYE_COLOURS <- c("#8B0000", "#FF0000", "#FFA500", "#FFFF00")   # BH 0 -> 0.1
BULLSEYE_CUTOFF <- 0.1

# Colour for a BH-corrected p-value on the bullseye scale: dark red at 0 through
# red and orange to yellow at 0.1, NA above 0.1 (drawn grey in the plot, left
# unshaded in the table). `lighten` blends towards white so the same hue can sit
# behind black text in the results table (0 = full plot colour, 1 = white).
bullseye_colour <- function(bh, lighten = 0) {
  ramp <- grDevices::colorRamp(BULLSEYE_COLOURS)
  out <- rep(NA_character_, length(bh))
  ok <- !is.na(bh) & bh <= BULLSEYE_CUTOFF
  if (any(ok)) {
    rgb <- ramp(pmin(pmax(bh[ok], 0), BULLSEYE_CUTOFF) / BULLSEYE_CUTOFF)
    rgb <- rgb * (1 - lighten) + 255 * lighten
    out[ok] <- grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
  }
  out
}

hex_polygon <- function(cx, cy, r) {
  a <- seq(0, 2 * pi, length.out = 7)[-7]          # flat-top hexagon
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
}

# tab: long result table from sea_fisher(); lay: output of layout_nodes().
bullseye_polygons <- function(tab, lay) {
  nodes <- lay$nodes
  tab <- tab[tab$code %in% nodes$code, ]
  tab <- tab[order(-tab$threshold), ]              # outer rings first, drawn underneath
  i <- match(tab$code, nodes$code)
  r <- lay$k * sqrt(tab$n_specific)
  polys <- lapply(seq_len(nrow(tab)), function(j) {
    p <- hex_polygon(nodes$x[i[j]], nodes$y[i[j]], r[j])
    p$id <- paste(tab$code[j], tab$threshold[j])
    p$bh <- tab$bh[j]
    p$threshold <- tab$threshold[j]
    p
  })
  out <- do.call(rbind, polys)
  # ggplot2 draws polygon groups in factor-level order, so fix the levels to
  # the row order (outer ring first) or the 0.05 ring would cover the others
  out$id <- factor(out$id, levels = unique(out$id))
  out$fill <- ifelse(out$bh <= BULLSEYE_CUTOFF, out$bh, NA_real_)
  # the TSEA tool draws the outermost (0.05) ring as a faint outline unless it
  # is significant; the other panels outline every ring in black
  out$faint <- isTRUE(lay$outer_faint) & out$threshold == 0.05 & is.na(out$fill)
  out
}

bullseye_plot <- function(tab, lay, base_size = 14, tree_dir = ".") {
  polys <- bullseye_polygons(tab, lay)
  nodes <- lay$nodes
  t05 <- tab[tab$threshold == 0.05, ]
  r05 <- lay$k * sqrt(t05$n_specific[match(nodes$code, t05$code)])
  top <- nodes$y + r05 * sqrt(3) / 2
  bottom <- nodes$y - r05 * sqrt(3) / 2
  fw <- lay$frame[["width"]]; fh <- lay$frame[["height"]]
  pad <- 0.012 * fh

  p <- ggplot2::ggplot()
  if (!is.null(lay$tree)) {
    f <- file.path(tree_dir, lay$tree)
    if (file.exists(f)) {
      p <- p + ggplot2::annotation_raster(png::readPNG(f), xmin = 0, xmax = fw, ymin = 0, ymax = fh)
    }
  }
  p <- p +
    ggplot2::geom_polygon(data = polys[polys$faint, ], ggplot2::aes(x, y, group = id),
                          fill = "white", colour = "#DCDCDC", linewidth = 0.6) +
    ggplot2::geom_polygon(data = polys[!polys$faint, ], ggplot2::aes(x, y, group = id, fill = fill),
                          colour = "black", linewidth = 0.4) +
    ggplot2::scale_fill_gradientn(
      colours = BULLSEYE_COLOURS, limits = c(0, BULLSEYE_CUTOFF), na.value = BULLSEYE_GREY,
      breaks = c(0, 0.025, 0.05, 0.075, 0.1), name = "BH-corrected p-value",
      guide = ggplot2::guide_colourbar(reverse = TRUE, direction = "horizontal", title.position = "top",
                                       title.hjust = 0.5, barwidth = grid::unit(0.45, "npc"),
                                       barheight = grid::unit(0.6, "lines"))) +
    ggplot2::coord_fixed(xlim = c(0, fw), ylim = c(0, fh), expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(4, 4, 4, 4),
                   legend.title = ggplot2::element_text(size = base_size * 0.8),
                   legend.text = ggplot2::element_text(size = base_size * 0.7))

  lab <- nodes[!is.na(nodes$label), ]
  if (nrow(lab)) {
    if (isTRUE(lay$label_above)) {
      lab$ty <- top[!is.na(nodes$label)] + pad
      p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x, ty, label = label), vjust = 0,
                                  size = base_size / ggplot2::.pt * 0.74, lineheight = 0.85)
      sub <- lab[!is.na(lab$sublabel), ]
      if (nrow(sub)) {
        sub$sy <- bottom[match(sub$code, nodes$code)] - pad
        p <- p + ggplot2::geom_text(data = sub, ggplot2::aes(x, sy, label = sublabel), vjust = 1,
                                    size = base_size / ggplot2::.pt * 0.74)
      }
    } else {
      lab$ty <- bottom[!is.na(nodes$label)] - pad
      p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x, ty, label = label), vjust = 1,
                                  size = base_size / ggplot2::.pt * 0.74)
    }
  }
  if (!is.null(lay$headers)) {
    p <- p + ggplot2::geom_text(data = lay$headers, ggplot2::aes(x, y, label = text),
                                size = base_size / ggplot2::.pt * 0.95, lineheight = 0.85)
  }
  p
}
