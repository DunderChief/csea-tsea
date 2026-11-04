# Export the Shinylive site to ./site and (optionally) serve it locally.
#   Rscript dev/export_site.R          # export only (this is what .github/workflows/deploy.yaml runs)
#   Rscript dev/export_site.R serve    # export, then serve on http://localhost:8008
args <- commandArgs(trailingOnly = TRUE)

# The page Shinylive generates shows only a spinner while R and its packages download (about 45 MB on a
# first visit, a minute or two). Give the page a title and a line of explanation; the message is removed as
# soon as the app's iframe appears.
loading_html <- '
<div id="sea-loading" style="position: fixed; left: 0; right: 0; top: 58%; padding: 0 1rem; text-align: center;
     font-family: system-ui, -apple-system, sans-serif; color: #555; font-size: 0.95rem; line-height: 1.4;">
  <p style="margin: 0 0 0.3rem;"><strong>Loading CSEA / TSEA</strong></p>
  <p style="margin: 0;">R is starting in your browser. The first visit downloads about 45 MB and can take a
  minute or two; later visits take a few seconds.</p>
</div>
<script>
  // Shinylive creates the iframe at once; the app is there when its navbar exists inside it.
  var seaLoadingTimer = setInterval(function () {
    var frame = document.querySelector("#root iframe"), ready = false;
    try { ready = !!(frame && frame.contentDocument && frame.contentDocument.querySelector(".navbar")); }
    catch (e) { ready = true; }
    if (ready) { var el = document.getElementById("sea-loading"); if (el) el.remove(); clearInterval(seaLoadingTimer); }
  }, 300);
</script>'

shinylive::export("app", "site", quiet = FALSE,
  template_params = list(title = "CSEA / TSEA: Specific Expression Analysis",
                         include_after_body = loading_html))

# webR fetches R.wasm and the base-library image first, then the package binaries one after another, then the
# translations file. Preload hints let the browser download those later files from the first moment, in parallel
# with the runtime (measured against GitHub Pages: parallel downloads take about half the time of sequential
# ones; in a normal profile the worker's later request reuses the preload, even while it is still in flight).
# The runtime files are left out: webR asks for them within a second anyway, and a private/incognito window,
# whose memory-only cache rejects entries of more than a few MB, would download them twice.
webr_dir <- file.path("site", "shinylive", "webr")
candidates <- c(list.files(file.path(webr_dir, "packages"), pattern = "\\.tgz$", recursive = TRUE, full.names = TRUE),
                file.path(webr_dir, "vfs/usr/lib/R/library/translations.data.gz"))
candidates <- candidates[file.exists(candidates)]
big <- candidates[file.size(candidates) > 1e6]
big <- big[order(-file.size(big))]
preload <- sprintf('    <link rel="preload" href="./%s" as="fetch" crossorigin="anonymous" />', sub("^site/", "", big))
index <- file.path("site", "index.html")
html <- readLines(index, warn = FALSE)
i <- grep("</head>", html, fixed = TRUE)[1]
writeLines(c(html[seq_len(i - 1)], preload, html[i:length(html)]), index)
cat("preload hints for", length(big), "files:", paste(basename(big), collapse = ", "), "\n")

cat("site size:", round(sum(file.size(list.files("site", recursive = TRUE, full.names = TRUE))) / 1e6), "MB\n")
if ("serve" %in% args) httpuv::runStaticServer("site", port = 8008)
