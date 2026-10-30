# Recover the hand-drawn bullseye layouts of the old tools from their saved PNGs.
#
# For a reference image this script
#   1. finds the hexagon fill regions (grey 211,211,211) as connected components,
#   2. keeps those that contain a smaller, co-centred component (the nested
#      threshold rings) and are hexagon-shaped (bbox aspect ~ 2/sqrt(3)),
#      which separates hexagons from the decorative tree drawn in the same grey,
#   3. prints the hexagon centres and outer widths (pixels, origin top-left),
#      to be labelled by hand in app/R/layouts.R,
#   4. writes the remaining grey pixels (the tree) as a transparent PNG that the
#      app draws underneath the vector hexagons.
#
# Usage: Rscript dev/extract_layout.R image.png out_tree.png csea|tsea
args <- commandArgs(trailingOnly = TRUE)
img <- png::readPNG(args[[1]])
r <- round(img[, , 1] * 255); g <- round(img[, , 2] * 255); b <- round(img[, , 3] * 255)
grey <- (r == g & g == b)
mask <- grey & r == 211
cat("image:", ncol(r), "x", nrow(r), "px;", sum(mask), "fill-grey pixels\n")

# connected components (4-neighbour) by iterative minimum-label propagation
lab <- matrix(seq_along(mask), nrow = nrow(mask)); lab[!mask] <- NA
repeat {
  up <- rbind(NA, lab[-nrow(lab), ]); dn <- rbind(lab[-1, ], NA)
  lf <- cbind(NA, lab[, -ncol(lab)]); rt <- cbind(lab[, -1], NA)
  new <- pmin(lab, up, dn, lf, rt, na.rm = TRUE); new[!mask] <- NA
  if (identical(new, lab)) break
  lab <- new
}
ids <- sort(unique(as.vector(lab[!is.na(lab)])))
comp <- do.call(rbind, lapply(ids, function(i) {
  w <- which(lab == i, arr.ind = TRUE)
  data.frame(id = i, n = nrow(w), x0 = min(w[, 2]), x1 = max(w[, 2]), y0 = min(w[, 1]), y1 = max(w[, 1]))
}))
comp$cx <- (comp$x0 + comp$x1) / 2; comp$cy <- (comp$y0 + comp$y1) / 2
comp$w <- comp$x1 - comp$x0 + 1; comp$h <- comp$y1 - comp$y0 + 1
comp <- comp[comp$n >= 20, ]

is_hex <- vapply(seq_len(nrow(comp)), function(i) {
  a <- comp[i, ]
  if (abs(a$w / a$h - 2 / sqrt(3)) > 0.15) return(FALSE)
  inner <- comp[-i, ]
  any(abs(inner$cx - a$cx) <= 3 & abs(inner$cy - a$cy) <= 3 &
      inner$w < a$w & inner$x0 > a$x0 & inner$x1 < a$x1)
}, logical(1))
hex <- comp[is_hex, ]
hex <- hex[order(-hex$w), ]
keep <- rep(TRUE, nrow(hex))
for (i in seq_len(nrow(hex))) for (j in seq_len(i - 1)) {
  if (keep[j] && abs(hex$cx[i] - hex$cx[j]) <= 3 && abs(hex$cy[i] - hex$cy[j]) <= 3) keep[i] <- FALSE
}
hex <- hex[keep, ]
hex <- hex[order(hex$cy, hex$cx), ]
rownames(hex) <- NULL
cat(nrow(hex), "hexagons (cx, cy, outer width):\n")
for (i in seq_len(nrow(hex))) cat(sprintf("  %2d  x=%6.1f  y=%6.1f  w=%3d\n", i, hex$cx[i], hex$cy[i], hex$w[i]))

# tree = grey pixels outside every hexagon (0.05 ring plus a margin), plus 1px
# of anti-aliased neighbours. Hexagon centres/sizes come from app/R/layouts.R
# and the panel data, so the mask is independent of the detection above.
source("app/R/layouts.R"); source("app/R/stats.R")
panel <- args[[3]]
dat <- readRDS(file.path("app/data", c(csea = "csea_mouse", tsea = "tsea_tissue")[[panel]], fsep = "/") |> paste0(".rds"))
n05 <- colSums(!is.na(dat$psi) & dat$psi < 0.05)
lay <- layout_nodes(panel, dat$columns, n05)
yy <- matrix(rep(seq_len(nrow(r)), ncol(r)), nrow = nrow(r)); xx <- matrix(rep(seq_len(ncol(r)), each = nrow(r)), nrow = nrow(r))
inside <- matrix(FALSE, nrow(r), ncol(r))
for (i in seq_len(nrow(lay$nodes))) {
  cx <- lay$nodes$x[i]; cy <- lay$frame[["height"]] - lay$nodes$y[i]; R <- lay$k * sqrt(n05[[lay$nodes$code[i]]]) + 4
  dx <- abs(xx - cx); dy <- abs(yy - cy)
  inside <- inside | (dy <= sqrt(3) / 2 * R & sqrt(3) * dx + dy <= sqrt(3) * R)
}
tree <- mask & !inside
dil <- tree
dil[-1, ] <- dil[-1, ] | tree[-nrow(tree), ]; dil[-nrow(dil), ] <- dil[-nrow(dil), ] | tree[-1, ]
dil[, -1] <- dil[, -1] | tree[, -ncol(tree)]; dil[, -ncol(dil)] <- dil[, -ncol(dil)] | tree[, -1]
alpha <- pmin(1, pmax(0, (255 - r) / (255 - 211)))
alpha[!(dil & grey & r >= 205 & !inside)] <- 0
out <- array(0, c(nrow(r), ncol(r), 4))
out[, , 1:3] <- 211 / 255
out[, , 4] <- alpha
png::writePNG(out, args[[2]])
cat("tree written to", args[[2]], "(", round(file.size(args[[2]]) / 1024), "KB,", sum(alpha > 0), "px )\n")
