# Export the Shinylive site to ./site and (optionally) serve it locally.
#   Rscript dev/export_site.R          # export only
#   Rscript dev/export_site.R serve    # export, then serve on http://localhost:8008
args <- commandArgs(trailingOnly = TRUE)
shinylive::export("app", "site", quiet = FALSE)
cat("\nsite size:", round(sum(file.size(list.files("site", recursive = TRUE, full.names = TRUE))) / 1e6), "MB\n")
if ("serve" %in% args) httpuv::runStaticServer("site", port = 8008)
