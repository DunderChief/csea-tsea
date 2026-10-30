# Cross-check sea_fisher() against pSI::fisher.iteration() / candidate.overlap()
# on every panel, for the fixture list and the five bundled example lists.
# Usage: Rscript dev/compare_with_pSI.R [extra_lib_path]
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(pSI))
source("app/R/stats.R")

lists <- c(list(fixture20 = readLines("tests/fixtures/genes_20.txt")),
           setNames(lapply(readRDS("app/data/examples.rds"), `[[`, "genes"),
                    vapply(readRDS("app/data/examples.rds"), `[[`, "", "id")))
panels <- c("csea_mouse", "sea_adult", "sea_development", "tsea_tissue")
worst <- 0; n_cells <- 0; ok_overlap <- TRUE
for (panel in panels) {
  psi <- readRDS(file.path("app/data", paste0(panel, ".rds")))$psi
  df <- as.data.frame(psi)
  for (nm in names(lists)) {
    g <- lists[[nm]]
    ours <- sea_fisher(psi, g)$table
    ref <- fisher.iteration(df, g, background = "data.set", p.adjust = TRUE)
    ref_nom <- fisher.iteration(df, g, background = "data.set", p.adjust = FALSE)
    for (k in seq_along(SEA_THRESHOLDS)) {
      o <- ours[ours$threshold == SEA_THRESHOLDS[k], ]
      stopifnot(identical(o$code, rownames(ref)))
      worst <- max(worst, abs(o$bh - ref[, k]), abs(o$p - ref_nom[, k]))
      n_cells <- n_cells + nrow(o)
    }
    ov <- candidate.overlap(df, g)
    for (k in seq_along(SEA_THRESHOLDS)) {
      lst <- ov[[paste0("pSi_", format(SEA_THRESHOLDS[k], scientific = FALSE))]]
      o <- ours[ours$threshold == SEA_THRESHOLDS[k], ]
      for (j in seq_len(ncol(psi))) {
        theirs <- sort(as.character(stats::na.omit(lst[[j]])))
        mine <- sort(strsplit(o$overlap[j], ";")[[1]])
        if (!identical(theirs, mine)) { ok_overlap <- FALSE; message("overlap differs: ", panel, " ", nm, " ", o$code[j], " ", SEA_THRESHOLDS[k]) }
      }
    }
    cat(sprintf("%-16s %-10s max |diff| so far = %.2e\n", panel, nm, worst))
  }
}
cat(sprintf("\n%d cells compared; largest absolute difference in p or BH: %.2e; overlaps identical: %s\n",
            n_cells, worst, ok_overlap))
