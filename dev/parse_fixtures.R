# Parse the saved output pages of the old CGI tools into CSV regression fixtures
# (tests/fixtures/*.csv, already committed; rerun only if the saved pages change).
# The source pages are copies of the old tools' result pages saved on
# 2026-10-22/23 for the gene lists in tests/fixtures/genes_*.txt; they are not
# part of this repository.
# Usage: Rscript dev/parse_fixtures.R <directory with the saved result pages>

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript dev/parse_fixtures.R <directory with the saved result pages>")
ref_dir <- args[[1]]

pages <- c(
  csea_mouse      = "result_seaCSEA_human.html",
  sea_adult       = "result_seaAdultBrain_human.html",
  sea_development = "result_seaBrainRegion_human.html",
  tsea_tissue     = "result_TissueFisherTest_human.html",
  csea_mouse_flag = "result_seaCSEA_mouseflag.html"   # first 10 genes, "Mouse" flag (tests/fixtures/genes_10_mouseflag.txt)
)

num_re <- "[0-9]+\\.[0-9]+(e-[0-9]+)?"

parse_page <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  txt <- gsub("null device\\s+1", "", txt)     # stray R output embedded in one CSEA cell
  txt <- gsub("\\s+", " ", txt)

  # header line: "Of the inputed N genes, K genes exist ..., M genes do not, which are: A,B <br>"
  n_input <- as.integer(sub(".*Of the inputed ([0-9]+) genes.*", "\\1", txt))
  n_found <- as.integer(sub(".*<br> ([0-9]+) genes exist in.*", "\\1", txt))
  missing <- sub(".*which are: ?(?:<br> )?(.*?) ?<br>.*", "\\1", txt, perl = TRUE)
  missing <- trimws(strsplit(missing, ",")[[1]])

  tbl <- regmatches(txt, regexpr("<table>.*?</table>", txt))
  rows <- regmatches(tbl, gregexpr("<tr>.*?</tr>", tbl))[[1]]
  rows <- rows[grepl("&#40", rows, fixed = TRUE)]        # drop header rows (one page has two)
  thresholds <- c(0.05, 0.01, 0.001, 1e-04)
  res <- do.call(rbind, lapply(rows, function(r) {
    cells <- strsplit(r, "<td>", fixed = TRUE)[[1]][-1]
    name <- trimws(sub("</td>.*", "", cells[[1]]))
    vals <- lapply(cells[-1], function(cell) {
      cell <- gsub("href=#[0-9_]+", "", cell)
      as.numeric(regmatches(cell, gregexpr(num_re, cell))[[1]])
    })
    stopifnot(length(vals) == 4, all(lengths(vals) == 2))
    data.frame(row = name, threshold = thresholds,
               p = vapply(vals, `[`, 0, 1), bh = vapply(vals, `[`, 0, 2),
               stringsAsFactors = FALSE)
  }))

  # overlap sections, one per threshold
  secs <- strsplit(txt, "pSI threshold </B>", fixed = TRUE)[[1]][-1]
  ov <- do.call(rbind, lapply(secs, function(s) {
    th <- as.numeric(regmatches(s, regexpr(num_re, s)))
    m <- gregexpr("<a name = [0-9_]+ id = [0-9_]+> (.*?) (?:</a> )?: (.*?) <br>", s, perl = TRUE)
    hits <- regmatches(s, m)[[1]]
    data.frame(
      row = sub("<a name = [0-9_]+ id = [0-9_]+> (.*?) (?:</a> )?: (.*?) <br>", "\\1", hits, perl = TRUE),
      threshold = th,
      overlap = sub("<a name = [0-9_]+ id = [0-9_]+> (.*?) (?:</a> )?: (.*?) <br>", "\\2", hits, perl = TRUE),
      stringsAsFactors = FALSE)
  }))
  ov$overlap[ov$overlap == "none"] <- ""
  ov$overlap <- gsub(", ", ";", ov$overlap)
  res <- merge(res, ov, by = c("row", "threshold"), sort = FALSE)
  res <- res[order(match(res$row, unique(sub("</td>.*", "", sub("^<tr><td>", "", rows)))), -res$threshold), ]
  rownames(res) <- NULL
  list(table = res, meta = data.frame(n_input = n_input, n_found = n_found,
                                      missing = paste(missing, collapse = ";")))
}

dir.create("tests/fixtures", showWarnings = FALSE, recursive = TRUE)
meta <- list()
for (panel in names(pages)) {
  out <- parse_page(file.path(ref_dir, pages[[panel]]))
  write.csv(out$table, file.path("tests/fixtures", paste0(panel, ".csv")), row.names = FALSE)
  meta[[panel]] <- cbind(panel = panel, out$meta)
  message(panel, ": ", nrow(out$table) / 4, " rows, ", out$meta$n_found, "/", out$meta$n_input,
          " genes found, missing: ", out$meta$missing)
}
write.csv(do.call(rbind, meta), "tests/fixtures/meta.csv", row.names = FALSE)
