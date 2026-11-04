# Regression fixtures: outputs saved from the old CGI tools (doughertytools.wustl.edu)
# on 2026-10-22/23 for the 20-gene list in tests/fixtures/genes_20.txt.
# CSEA/SEA pages print p-values to three decimals; the TSEA page prints full
# doubles from a different (Python) Fisher implementation, so its values of
# 0.99999999999x correspond to R's exact 1.

root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
genes20 <- readLines(file.path(root, "tests/fixtures/genes_20.txt"))
meta <- read.csv(file.path(root, "tests/fixtures/meta.csv"), stringsAsFactors = FALSE)

panels <- c("csea_mouse", "csea_mouse_flag", "sea_adult", "sea_development", "tsea_tissue")
# fixture -> (data file, gene list, species selector of the old form)
fixture_spec <- list(
  csea_mouse      = list(data = "csea_mouse", genes = "genes_20.txt", species = "human"),
  csea_mouse_flag = list(data = "csea_mouse", genes = "genes_10_mouseflag.txt", species = "mouse"),
  sea_adult       = list(data = "sea_adult", genes = "genes_20.txt", species = "human"),
  sea_development = list(data = "sea_development", genes = "genes_20.txt", species = "human"),
  tsea_tissue     = list(data = "tsea_tissue", genes = "genes_20.txt", species = "human")
)

# The old CGI used slightly different column names for two of the panels.
old_name <- function(panel, code) {
  if (panel %in% c("csea_mouse", "csea_mouse_flag")) {
    map <- c(BS.chat = "BS.Chat", Cb.septin4 = "Cb.Septin4", Cb.pcp2 = "Cb.Pcp2", Cb.grp = "Cb.Grp",
             Cb.Aldh1L1 = "Cb.Fthfd", Cb.cnp = "Cb.Cnp", Ctx.pdgfrajd340 = "Ctx.Pdgfrajd340",
             Ctx.etv1_ts88 = "Ctx.Etv1_ts88", Ctx.pnoc = "Ctx.Pnoc", Ctx.glt25d2 = "Ctx.Glt25d2",
             Ctx.Aldh1L1 = "Ctx.Fthfd", Ctx.cort = "Ctx.Cort", Ctx.cnp = "Ctx.Cnp")
    return(ifelse(code %in% names(map), map[code], code))
  }
  if (panel == "sea_development") {
    code <- sub("Neonatal...Early.Infancy", "Neotal.Early.Infancy", code, fixed = TRUE)
    code <- sub("Middle...Late.Childhood", "Middle.Late.Childhood", code, fixed = TRUE)
    return(code)
  }
  if (panel == "tsea_tissue") return(gsub(".", " ", code, fixed = TRUE))
  code
}

for (panel in panels) {
  test_that(paste(panel, "reproduces the saved CGI output"), {
    spec <- fixture_spec[[panel]]
    dat <- readRDS(file.path(root, "app/data", paste0(spec$data, ".rds")))
    fx <- read.csv(file.path(root, "tests/fixtures", paste0(panel, ".csv")), stringsAsFactors = FALSE)
    fx$overlap[is.na(fx$overlap)] <- ""
    genes <- readLines(file.path(root, "tests/fixtures", spec$genes))
    # the old form's Human/Mouse selector chose pSI's background: the whole table
    # for a same-species list, the genes with a homolog for the other species
    res <- sea_fisher(dat$psi, genes, total = dat$universe[[spec$species]])
    tab <- res$table
    tab$row <- old_name(panel, tab$code)

    m <- meta[meta$panel == panel, ]
    expect_equal(res$n_genes, m$n_input)
    expect_equal(length(res$found), m$n_found)
    expect_setequal(res$missing, strsplit(m$missing, ";")[[1]])

    expect_setequal(unique(tab$row), unique(fx$row))
    key <- paste(fx$row, fx$threshold)
    idx <- match(key, paste(tab$row, tab$threshold))
    expect_false(anyNA(idx))

    ours <- lapply(strsplit(tab$overlap[idx], ";"), sort)
    theirs <- lapply(strsplit(fx$overlap, ";"), sort)
    expect_equal(ours, theirs)

    if (panel != "tsea_tissue") {
      tol <- 5e-4 + 1e-9   # the CSEA/SEA pages print three decimals
      expect_true(all(abs(tab$p[idx] - fx$p) <= tol),
                  info = paste("nominal p mismatch at", key[abs(tab$p[idx] - fx$p) > tol]))
      expect_true(all(abs(tab$bh[idx] - fx$bh) <= tol),
                  info = paste("BH mismatch at", key[abs(tab$bh[idx] - fx$bh) > tol]))
    } else {
      # The old TSEA CGI (a separate Python script) counted one gene too many in
      # every tissue's specific list except Testis. Its values are reproduced
      # exactly by that off-by-one, and our values from the published table are
      # within 1% (relative) of them. Both facts are asserted so that the
      # documented difference stays exactly as documented.
      rel <- abs(tab$p[idx] - fx$p) / pmax(fx$p, 1e-12)
      expect_true(all(rel <= 0.01), info = paste("TSEA p differs by >1% at", key[rel > 0.01]))
      in_list <- rownames(dat$psi) %in% toupper(genes)
      p_old <- vapply(seq_along(idx), function(i) {
        z <- tab[idx[i], ]
        spec <- !is.na(dat$psi[, z$code]) & dat$psi[, z$code] < z$threshold
        a <- sum(spec & in_list); n_spec <- sum(spec) + if (z$code == "Testis") 0 else 1
        phyper(a - 1, n_spec, nrow(dat$psi) - n_spec, sum(in_list), lower.tail = FALSE)
      }, numeric(1))
      expect_true(all(abs(p_old - fx$p) <= 1e-9),
                  info = paste("off-by-one model fails at", key[abs(p_old - fx$p) > 1e-9]))
    }
  })
}

test_that("the CSEA human-list background is the historic pSI constant", {
  dat <- readRDS(file.path(root, "app/data/csea_mouse.rds"))
  expect_equal(dat$universe$human, 13421)
  expect_equal(dat$universe$mouse, nrow(dat$psi))
  # with the full table instead, the old page is NOT reproduced (guards against silently dropping the background)
  fx <- read.csv(file.path(root, "tests/fixtures/csea_mouse.csv"), stringsAsFactors = FALSE)
  full <- sea_fisher(dat$psi, genes20)$table
  expect_gt(abs(full$p[full$code == "Ctx.Ntsr" & full$threshold == 0.05] - fx$p[fx$row == "Ctx.Ntsr" & fx$threshold == 0.05]), 0.005)
})

test_that("parse_gene_list handles separators, case and duplicates", {
  g <- parse_gene_list("Shank3, scn2a\nSCN2A;CHD8\tchd8  PTEN\n\n")
  expect_equal(g$genes, c("SHANK3", "SCN2A", "CHD8", "PTEN"))
  expect_equal(g$n_entered, 6)
  expect_equal(g$n_duplicates, 2)
  expect_equal(parse_gene_list("")$genes, character(0))
  expect_equal(parse_gene_list(NULL)$genes, character(0))
})

test_that("sea_fisher is a one-sided Fisher test on the right 2x2 table", {
  psi <- matrix(c(0.001, 0.2, NA, 0.04, 0.5, 0.0005, NA, NA, 0.02, 0.3), ncol = 2,
                dimnames = list(c("A", "B", "C", "D", "E"), c("X", "Y")))
  res <- sea_fisher(psi, c("a", "b", "zzz"))
  expect_equal(res$missing, "ZZZ")
  expect_equal(res$found, c("A", "B"))
  r <- res$table[res$table$code == "X" & res$table$threshold == 0.05, ]
  # specific in X at 0.05: A, D (NA excluded); list: A, B -> a=1, b=1, c=1, d=2
  expect_equal(r$n_specific, 2); expect_equal(r$n_overlap, 1); expect_equal(r$overlap, "A")
  expect_equal(r$p, fisher.test(matrix(c(1, 1, 1, 2), 2), alternative = "greater")$p.value)
  expect_equal(r$p, phyper(0, 2, 3, 2, lower.tail = FALSE))
})

test_that("sea_wide formats one cell per threshold", {
  dat <- readRDS(file.path(root, "app/data/sea_adult.rds"))
  w <- sea_wide(sea_fisher(dat$psi, genes20)$table)
  expect_equal(names(w), c("code", "0.05", "0.01", "0.001", "0.0001"))
  expect_equal(w[w$code == "Cortex.Young.Adulthood", "0.0001"], "0.007 (0.040)")
})

test_that("table shading uses the bullseye colour scale", {
  expect_true(is.na(bullseye_colour(0.2)))
  expect_true(is.na(bullseye_colour(NA)))
  expect_equal(bullseye_colour(0), "#8B0000")        # darkred at 0
  expect_equal(bullseye_colour(1e-5), "#8B0000")     # values below 1e-5 are clamped: the same darkest red
  expect_equal(bullseye_colour(0.1), "#EDED00")      # yellow2 (through Lab space) at the cutoff
  expect_equal(bullseye_colour(c(0.2, 0.05))[1], NA_character_)
  # the original tools' binning: cut(c(1e-5, 0.1, p), breaks = 128) labelled with the reversed 128-colour ramp
  ramp <- grDevices::colorRampPalette(c("yellow2", "#FF6700", "#D11500", "darkred"), bias = 0.15, space = "Lab")(128)
  p <- c(2e-5, 0.0007, 0.03, 0.049, 0.0999)
  expect_equal(bullseye_colour(p), as.character(cut(c(1e-5, 0.1, p), breaks = 128, labels = rev(ramp)))[-(1:2)])
  light <- grDevices::col2rgb(bullseye_colour(0, lighten = 0.55))
  expect_true(all(light > grDevices::col2rgb("#8B0000")))   # blended towards white
  expect_match(bullseye_colour(0.03, lighten = 0.55), "^#[0-9A-F]{6}$")
})
