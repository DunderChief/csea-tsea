# The CSEA/TSEA layouts generated from the lab's original node coordinates.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))

test_that("CSEA layout tables are complete and consistent", {
  expect_equal(nrow(CSEA_NODES), 27)
  expect_true(all(CSEA_NODES$code %in% CSEA_TREE_NODES$node))
  expect_true(all(c(CSEA_TREE_EDGES$from, CSEA_TREE_EDGES$to) %in% CSEA_TREE_NODES$node))
  expect_equal(nrow(CSEA_TREE_EDGES), nrow(CSEA_TREE_NODES) - 1)   # a tree
  expect_false(anyNA(CSEA_NODES$label))
})

test_that("Layer 6 and Layer 5b labels are on the right cell types (the 2017 fix)", {
  lab <- setNames(CSEA_NODES$label, CSEA_NODES$code)
  expect_match(lab[["Ctx.Ntsr"]], "^Layer 6")
  expect_match(lab[["Ctx.glt25d2"]], "^Layer 5b")
})

test_that("TSEA layout tables are complete and consistent", {
  dat <- readRDS(file.path(root, "app/data/tsea_tissue.rds"))
  expect_true(all(dat$columns$code %in% TSEA_TREE_NODES$node))
  expect_true(all(c(TSEA_TREE_EDGES$from, TSEA_TREE_EDGES$to) %in% TSEA_TREE_NODES$node))
  expect_equal(nrow(TSEA_TREE_EDGES), nrow(TSEA_TREE_NODES) - 1)
  lay <- layout_nodes("tsea", dat$columns, colSums(!is.na(dat$psi) & dat$psi < 0.05))
  expect_equal(nrow(lay$nodes), 25)
  expect_equal(nrow(lay$segments), 48)
})
