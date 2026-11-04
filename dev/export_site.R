# Export the Shinylive site to ./site and (optionally) serve it locally.
#   Rscript dev/export_site.R          # export only (this is what .github/workflows/deploy.yaml runs)
#   Rscript dev/export_site.R serve    # export, then serve on http://localhost:8008
args <- commandArgs(trailingOnly = TRUE)

# The page Shinylive generates shows only a spinner while R and its packages download (about 30 MB on a
# first visit, a minute or two). Give the page a title and a line of explanation; the message is removed as
# soon as the app's iframe appears.
loading_html <- '
<div id="sea-loading" style="position: fixed; left: 0; right: 0; top: 58%; padding: 0 1rem; text-align: center;
     font-family: system-ui, -apple-system, sans-serif; color: #555; font-size: 0.95rem; line-height: 1.4;">
  <p style="margin: 0 0 0.3rem;"><strong>Loading CSEA / TSEA</strong></p>
  <p style="margin: 0;">R is starting in your browser. The first visit downloads about 30 MB and can take a
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

# Not done, measured and rejected (2026-11-04): <link rel="preload"> hints for the package binaries. The browser
# fetches everything over one HTTP/2 connection, so the hints only reorder the queue (and held back app.json,
# which is on the critical path); first-visit time was unchanged on GitHub Pages and under a throttled network.
# Load time is the download size; the lever that worked was fewer bytes: the bullseye was ported from ggplot2 to
# base graphics, which took ggplot2 and its dependencies (14 MB) out of the site.

cat("site size:", round(sum(file.size(list.files("site", recursive = TRUE, full.names = TRUE))) / 1e6), "MB\n")
if ("serve" %in% args) httpuv::runStaticServer("site", port = 8008)
