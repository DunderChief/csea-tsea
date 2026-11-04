# Render the four bullseyes for a gene list to PNG files, for eyeballing.
# Usage: Rscript dev/render_examples.R out_dir [list: fixture|autdb|rdnv|hcrt|retina]
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]; which <- if (length(args) > 1) args[[2]] else "fixture"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (f in list.files("app/R", full.names = TRUE)) source(f)
genes <- if (which == "fixture") readLines("tests/fixtures/genes_20.txt") else {
  ex <- readRDS("app/data/examples.rds"); ex[[match(which, vapply(ex, `[[`, "", "id"))]]$genes }
specs <- list(csea = c("csea_mouse", "csea"), adult = c("sea_adult", "row"),
              dev = c("sea_development", "grid"), tsea = c("tsea_tissue", "tsea"))
for (nm in names(specs)) {
  dat <- readRDS(file.path("app/data", paste0(specs[[nm]][1], ".rds")))
  res <- sea_fisher(dat$psi, genes)
  n05 <- setNames(res$table$n_specific[res$table$threshold == 0.05], res$table$code[res$table$threshold == 0.05])
  lay <- layout_nodes(specs[[nm]][2], dat$columns, n05)
  w <- 1100; h <- round(w * bullseye_aspect(lay))
  png(file.path(out, paste0(nm, "_", which, ".png")), width = w, height = h, res = 96)
  bullseye_plot(res$table, lay); dev.off()
  cat(nm, ": ", w, "x", h, " min BH = ", signif(min(res$table$bh), 3), "\n", sep = "")
}
