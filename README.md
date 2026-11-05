# CSEA / TSEA: Specific Expression Analysis in the browser

**Live site: <https://dunderchief.github.io/csea-tsea/>**

A rebuild of the Dougherty lab's **CSEA** (Cell-type Specific Expression Analysis) and **TSEA** (Tissue Specific
Expression Analysis) web tools (formerly `doughertytools.wustl.edu/CSEAtool.html` and `TSEAtool.html`) as a
Shiny app compiled with [Shinylive](https://github.com/posit-dev/r-shinylive) and served as static files from
GitHub Pages. R runs in the browser through WebAssembly: no server, no accounts, no bills, and gene lists never
leave the user's machine.

The four analyses, their inputs, the statistics (one-sided Fisher's exact test per column and pSI threshold,
Benjamini-Hochberg corrected) and the bullseye plots are those of the original tools; see `app/about.md` for the
exact definition, the data provenance and the documented differences.

| Panel | Data | Rows |
|---|---|---|
| CSEA: brain cell types | `pSI.data` mouse bacTRAP cell types (Xu et al. 2014) | 35 |
| SEA: adult brain regions | `pSI.data` human BrainSpan, young adulthood | 6 |
| SEA: brain regions across development | `pSI.data` human BrainSpan, 6 regions x 10 periods | 60 |
| TSEA: tissues | Table S3 of Wells et al. 2015 (GTEx, 25 tissues) | 25 |

## Repository layout

```
app/                 the Shiny app (this directory is what Shinylive exports)
  app.R              UI (bslib navbar, one tab per analysis + About) and server
  about.md           About page: citations, data provenance, differences from the old tools
  R/stats.R          parse_gene_list(), sea_fisher(): the statistics (plain R, matches pSI::fisher.iteration)
  R/bullseye.R       bullseye plot (base graphics: nested hexagons, the original tools' colour scale and legend)
  R/layouts.R        layouts: the CSEA/TSEA trees (R/layout_tables.R, generated from the lab's original node
                     coordinates by dev/layouts_from_originals.R), row and grid layouts
  R/panels.R         the four panel definitions and explanatory text
  R/mod_analysis.R   Shiny module: gene list -> summary, table, bullseye, overlapping genes, downloads
  data/*.rds         the pSI tables (xz-compressed, ~1.5 MB in total) and the example gene lists
data-raw/            prepare_data.R rebuilds app/data from the published tables (see data-raw/README.md)
tests/               regression fixtures parsed from the old tools' output pages, and testthat tests
dev/                 one-off scripts: parse the fixture pages, compare with the pSI package, regenerate the layout
                     tables from original_layouts/, render example plots, export/serve the site
                     original_bullseye_functions.R is the lab's original plotting code, kept for reference
.github/workflows/   GitHub Actions: run the tests, export the site, publish to GitHub Pages
```

## Run locally

```r
shiny::runApp("app")                       # ordinary Shiny, for development
```

```bash
Rscript tests/run_tests.R                  # regression tests (fixture p-values, overlaps, parsing)
Rscript dev/compare_with_pSI.R             # cross-check against pSI::fisher.iteration() (needs pSI installed)
Rscript dev/export_site.R serve            # export the Shinylive site to ./site and serve it on :8008
```

Requirements: R >= 4.2 with shiny, bslib, png, testthat, and shinylive (>= 0.3) for the export. The
`archive` dependency of shinylive needs libarchive; on Ubuntu `sudo apt install libarchive-dev`, or install the
prebuilt binary from Posit Package Manager.

## Deploy

Pushing to `main` runs `.github/workflows/deploy.yaml`, which runs the tests, exports the site and publishes it
to GitHub Pages (Settings > Pages > Source must be "GitHub Actions"). The site is then at
<https://dunderchief.github.io/csea-tsea/>.

## Verification status

- `sea_fisher()` reproduces `pSI::fisher.iteration()` and `pSI::candidate.overlap()` exactly on all four panels for
  the regression list and the five bundled example lists (3,024 cells, largest difference 0).
- The CSEA and both SEA panels reproduce the saved output pages of the old CGI tools to the three decimals they
  print, including the lists of overlapping genes, for both settings of the Human/Mouse selector (a human list
  against the mouse table uses pSI's homolog-restricted background of 13,421 genes, the value the old tool used).
- The old TSEA page differs by an off-by-one in its specific-gene counts (see `tests/testthat/test-stats.R`);
  the tests assert both that our values match the published table and that adding one to each count
  reproduces the old page exactly.
- `dev/webr_validate.mjs` runs the data loading, statistics and plotting under webR (R compiled to WebAssembly)
  in Node and reproduces the native numbers; per panel the data load in under 0.1 s, the tests take 0.1-0.3 s and
  a bullseye renders in under a second.

## Citation

Xu X, Wells AB, O'Brien DR, Nehorai A, Dougherty JD. *J Neurosci* 2014;34(4):1420-31 (CSEA).
Wells A, Kopp N, Xu X, O'Brien DR, et al. *Nucleic Acids Res* 2015;43(22):10804-20 (TSEA).
Dougherty JD, Schmidt EF, Nakajima M, Heintz N. *Nucleic Acids Res* 2010;38(13):4218-30 (specificity index).

License: GPL (>= 2), like the pSI package the statistics come from.
