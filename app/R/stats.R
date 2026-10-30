# Specific Expression Analysis statistics.
#
# This is the computation behind the CSEA / TSEA web tools, written against the
# behaviour of pSI::fisher.iteration() (pSI 1.1, Xu et al. 2014) so that the
# p-values are identical:
#
#   * for each column (cell type / region / tissue) and each pSI threshold t,
#     the "specific" gene set is the rows with pSI < t (NA = not specific);
#   * the candidate list is matched to the table by upper-cased symbol, with
#     duplicates counted once; candidates absent from the table are reported
#     and otherwise ignored;
#   * a one-sided (greater) Fisher's exact test on the 2x2 table
#       [in list & specific, in list & not specific;
#        not in list & specific, not in list & not specific]
#     with the whole table (all rows) as the gene universe, or, when the
#     candidate list comes from the other species, the number of genes in the
#     table with a clear homolog (pSI's background = "mouse.human" /
#     "human.mouse"; only the universe size changes, the specific lists and
#     the matching do not);
#   * Benjamini-Hochberg adjustment across columns, separately per threshold.
#
# The functions here are plain R (no package dependencies) so they run
# unchanged under webR. Once pSI 1.2 is on CRAN and webR, sea_fisher() can
# delegate the p-values to pSI::fisher.iteration(); the regression tests in
# tests/ compare the two.

SEA_THRESHOLDS <- c(0.05, 0.01, 0.001, 0.0001)

# Split free text into gene symbols. Accepts one symbol per line, or symbols
# separated by commas, semicolons, tabs or spaces. Returns the symbols as
# entered plus the upper-cased, de-duplicated list used for matching.
parse_gene_list <- function(text) {
  text <- if (is.null(text) || is.na(text)) "" else text
  raw <- unlist(strsplit(text, "[[:space:],;]+"))
  raw <- raw[nzchar(raw)]
  up <- toupper(raw)
  dup <- duplicated(up)
  list(
    as_entered = raw,
    genes = up[!dup],
    n_entered = length(raw),
    n_duplicates = sum(dup)
  )
}

# One-sided Fisher's exact test for every column x threshold.
#
# psi:   numeric matrix, genes (upper-case symbols as row names) x columns.
# genes: character vector of candidate gene symbols (any case, may contain
#        duplicates).
# total: size of the gene universe; defaults to all genes in the table
#        (pSI background = "data.set").
#
# Returns a list:
#   table    data.frame with one row per column x threshold: code, threshold,
#            n_specific (genes with pSI < threshold), n_overlap, p (nominal),
#            bh (Benjamini-Hochberg adjusted across columns at that
#            threshold), overlap (character, ";"-separated overlapping genes)
#   found    candidate genes present in the table (upper case)
#   missing  candidate genes absent from the table (upper case)
#   n_genes  number of distinct candidate genes submitted
#   total    size of the gene universe (rows of psi)
sea_fisher <- function(psi, genes, total = nrow(psi), thresholds = SEA_THRESHOLDS) {
  stopifnot(is.matrix(psi), is.numeric(psi), !is.null(rownames(psi)), total >= 1, total <= nrow(psi))
  genes <- unique(toupper(trimws(genes)))
  genes <- genes[nzchar(genes)]
  universe <- rownames(psi)
  in_list <- universe %in% genes
  found <- universe[in_list]
  missing <- genes[!genes %in% universe]
  n_list <- sum(in_list)

  rows <- lapply(thresholds, function(th) {
    specific <- !is.na(psi) & psi < th             # logical matrix
    n_specific <- colSums(specific)
    a <- colSums(specific & in_list)               # in list and specific
    b <- n_list - a                                # in list, not specific
    c_ <- n_specific - a                           # specific, not in list
    d <- total - a - b - c_                        # neither
    p <- vapply(seq_len(ncol(psi)), function(j) {
      stats::fisher.test(matrix(c(a[j], c_[j], b[j], d[j]), nrow = 2),
                         alternative = "greater")$p.value
    }, numeric(1))
    overlap <- vapply(seq_len(ncol(psi)), function(j) {
      paste(universe[specific[, j] & in_list], collapse = ";")
    }, character(1))
    data.frame(
      code = colnames(psi), threshold = th,
      n_specific = unname(n_specific), n_overlap = unname(a),
      p = p, bh = stats::p.adjust(p, method = "BH"),
      overlap = overlap, stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(match(tab$code, colnames(psi)), -tab$threshold), ]
  rownames(tab) <- NULL
  list(table = tab, found = found, missing = missing, n_genes = length(genes), total = total)
}

# Wide view of the results, one row per column, the four thresholds as columns,
# each cell "p (bh)" formatted like the original tools. Used for the on-screen
# table; the CSV download uses the long table from sea_fisher().
sea_wide <- function(tab, digits = 3) {
  fmt <- function(x) ifelse(x < 10^-digits & x > 0, formatC(x, digits = 2, format = "e"),
                            formatC(x, digits = digits, format = "f"))
  cell <- paste0(fmt(tab$p), " (", fmt(tab$bh), ")")
  codes <- unique(tab$code)
  ths <- sort(unique(tab$threshold), decreasing = TRUE)
  out <- data.frame(code = codes, stringsAsFactors = FALSE)
  for (th in ths) {
    idx <- match(codes, tab$code[tab$threshold == th])
    out[[format(th, scientific = FALSE, drop0trailing = TRUE)]] <- cell[tab$threshold == th][idx]
  }
  out
}
