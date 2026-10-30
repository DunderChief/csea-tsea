# Whole-app test: every tab must analyse its own dataset (a lazy-evaluation bug
# once made all four tabs use the last panel's data).
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
genes20 <- readLines(file.path(root, "tests/fixtures/genes_20.txt"))

test_that("each tab uses its own panel data, layout and downloads", {
  skip_if_not_installed("shiny"); skip_if_not_installed("bslib")
  withr::local_dir(file.path(root, "app"))
  app <- shiny::shinyAppDir(".")
  # the summary line uses the original pages' dataset phrases (adult and dev share one),
  # so the row count of the table tells the two human panels apart
  expected <- c(csea = "brain cell type expression dataset", adult = "brain region and development expression dataset",
                dev = "brain region and development expression dataset", tsea = "tissue expression dataset")
  n_rows <- c(csea = 35, adult = 6, dev = 60, tsea = 25)
  shiny::testServer(app, {
    for (id in names(expected)) {
      set_inputs <- function(...) { args <- list(...); names(args) <- paste0(id, "-", names(args)); do.call(session$setInputs, args) }
      set_inputs(species = "human", genes = paste(genes20, collapse = "\n"))
      set_inputs(run = 1)
      summary_html <- as.character(output[[paste0(id, "-summary")]]$html)
      expect_match(summary_html, expected[[id]], fixed = TRUE, info = id)
      table_html <- as.character(output[[paste0(id, "-table")]]$html)
      # header row, plus two block-heading rows on the CSEA tab (cell types / whole regions)
      n_extra <- if (id == "csea") 3 else 1
      expect_equal(lengths(regmatches(table_html, gregexpr("<tr", table_html, fixed = TRUE))), n_rows[[id]] + n_extra, info = id)
      csv <- read.csv(output[[paste0(id, "-csv")]], stringsAsFactors = FALSE)
      expect_equal(nrow(csv), n_rows[[id]] * 4, info = id)
      expect_true("on_bullseye" %in% names(csv), info = id)
      if (id == "csea") {
        expect_equal(sum(csv$on_bullseye), 27 * 4)
        expect_true(all(diff(as.integer(csv$on_bullseye)) <= 0))   # plotted rows first, then whole regions
        expect_setequal(unique(csv$code[!csv$on_bullseye]), c("Hyp", "BF", "BS", "Cb", "Cpu", "Ctx", "Epi", "Spc"))
      } else expect_true(all(csv$on_bullseye), info = id)
      expect_true(file.exists(output[[paste0(id, "-csv")]]), info = paste(id, "csv"))
      expect_true(file.exists(output[[paste0(id, "-png")]]), info = paste(id, "png"))
      expect_gt(file.size(output[[paste0(id, "-png")]]), 10000)
    }
  })
})
