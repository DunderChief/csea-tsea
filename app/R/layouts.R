# Bullseye layouts.
#
# The CSEA (mouse cell types) and TSEA (tissues) bullseyes of the original
# tools place the hexagons on a hand-drawn tree. The hexagon centres below were
# measured from the tools' own output images (dev/extract_layout.R); the tree
# itself is shipped as a transparent PNG drawn underneath the vector hexagons.
# The adult-brain and developmental panels use a row and a region x period grid.
#
# Coordinates are pixels of the original images, origin top-left (y grows
# downwards); layout_nodes() flips y so that ggplot draws them the right way up.

CSEA_FRAME <- c(width = 1500, height = 1125)
CSEA_NODES <- data.frame(
  code = c("Ctx.Ntsr", "Cpu.D1", "Ctx.glt25d2", "Ctx.cort", "RetC.Cones", "Cpu.D2",
           "Cb.grp", "Cb.Lypd6", "RetR.Rods", "Ctx.etv1_ts88", "Ctx.cnp", "Cb.Neurod1",
           "Ctx.pdgfrajd340", "Cb.pcp2", "Cb.Grm2", "Epi.Chat", "Cb.cnp", "Ctx.Aldh1L1",
           "BF.Chat", "Hyp.Hcrt", "Spc.Chat", "Cb.Aldh1L1", "Cpu.Chat", "BS.chat",
           "BS.Slc6a4", "Cb.septin4", "Ctx.pnoc"),
  x = c(856, 1189.5, 712.5, 999, 500, 1314, 345, 194, 626.5, 807.5, 1369, 177.5, 554.5,
        858.5, 188, 427.5, 1367, 1048, 300.5, 641.5, 913.5, 1320.5, 283.5, 778, 590,
        1121.5, 451),
  y = c(201, 206.5, 207, 239.5, 287, 293.5, 298, 324, 363.5, 407.5, 474.5, 489, 563,
        609, 671, 671, 671, 699, 715, 751, 803.5, 863.5, 879.5, 886.5, 901, 902, 902.5),
  stringsAsFactors = FALSE
)
CSEA_R_MAX <- 87.5   # outer (pSI < 0.05) radius of the largest cell type, px

TSEA_FRAME <- c(width = 1000, height = 750)
TSEA_NODES <- data.frame(
  code = c("Vagina", "Skin", "Prostate", "Esophagus", "Pituitary", "Stomach", "Nerve",
           "Adrenal.Gland", "Brain", "Liver", "Testis", "Lung", "Thyroid", "Pancreas",
           "Heart", "Colon", "Muscle", "Adipose.Tissue", "Ovary", "Breast",
           "Fallopian.Tube", "Kidney", "Blood.Vessel", "Uterus", "Blood"),
  x = c(329.5, 432.5, 544.5, 259.5, 657, 716.5, 396, 500.5, 339.5, 625.5, 767, 349.5,
        463.5, 751, 602, 689.5, 519.5, 235.5, 410.5, 234, 654.5, 362, 439.5, 653.5, 544),
  y = c(110.5, 114, 140, 167, 168.5, 229, 238, 276.5, 284, 296.5, 315.5, 403, 424, 437,
        438.5, 465, 469, 471.5, 480.5, 551, 547.5, 603, 632, 634, 640),
  stringsAsFactors = FALSE
)
TSEA_R_MAX <- 55

# Node positions, label text and hexagon scale for one panel.
#
# columns: the panel's column table (code, label, group, ...).
# n05:     named vector, number of genes with pSI < 0.05 per column (hexagon
#          sizes are proportional to sqrt(list size), so the largest list in
#          the panel gets radius r_max).
#
# Returns list(nodes = data.frame(code, x, y, label, sublabel), frame, k, tree)
# where k converts sqrt(number of genes) to a radius in frame units.
layout_nodes <- function(layout, columns, n05) {
  if (layout == "csea") {
    nodes <- CSEA_NODES
    nodes$y <- CSEA_FRAME[["height"]] - nodes$y
    i <- match(nodes$code, columns$code)
    nodes$label <- vapply(columns$label[i], function(s) paste(strwrap(s, 16), collapse = "\n"), "")
    nodes$sublabel <- columns$group[i]
    return(list(nodes = nodes, frame = CSEA_FRAME, k = CSEA_R_MAX / sqrt(max(n05[nodes$code])),
                tree = "www/csea_tree.png", label_above = TRUE))
  }
  if (layout == "tsea") {
    nodes <- TSEA_NODES
    nodes$y <- TSEA_FRAME[["height"]] - nodes$y
    i <- match(nodes$code, columns$code)
    nodes$label <- columns$label[i]
    nodes$sublabel <- NA_character_
    return(list(nodes = nodes, frame = TSEA_FRAME, k = TSEA_R_MAX / sqrt(max(n05[nodes$code])),
                tree = "www/tsea_tree.png", label_above = FALSE, outer_faint = TRUE))
  }
  if (layout == "row") {
    n <- nrow(columns); sp <- 190
    nodes <- data.frame(code = columns$code, x = sp * (seq_len(n) - 0.5) + 180, y = 250,
                        label = columns$label, sublabel = columns$group, stringsAsFactors = FALSE)
    return(list(nodes = nodes, frame = c(width = sp * n + 360, height = 420),
                k = 92 / sqrt(max(n05)), tree = NULL, label_above = TRUE))
  }
  if (layout == "grid") {
    regions <- unique(columns$region)
    periods <- unique(columns$period[order(columns$period_order)])
    cs <- 167; rs <- 94
    col <- match(columns$region, regions); row <- match(columns$period, periods)
    nodes <- data.frame(code = columns$code, x = 250 + cs * (col - 1), y = 150 + rs * (row - 1),
                        label = NA_character_, sublabel = NA_character_, stringsAsFactors = FALSE)
    frame <- c(width = 250 + cs * (length(regions) - 1) + 150, height = 150 + rs * (length(periods) - 1) + 80)
    nodes$y <- frame[["height"]] - nodes$y
    headers <- data.frame(
      x = c(250 + cs * (seq_along(regions) - 1), rep(95, length(periods))),
      y = frame[["height"]] - c(rep(75, length(regions)), 150 + rs * (seq_along(periods) - 1)),
      text = c(regions, vapply(periods, function(s) paste(strwrap(s, 12), collapse = "\n"), "")),
      stringsAsFactors = FALSE)
    return(list(nodes = nodes, frame = frame, k = 56 / sqrt(max(n05)), tree = NULL,
                label_above = TRUE, headers = headers))
  }
  stop("unknown layout: ", layout)
}
