# CSEA / TSEA: Specific Expression Analysis in the browser.
# Shiny app compiled with Shinylive; everything runs client-side.

library(shiny)
library(bslib)
library(ggplot2)

if (!exists("sea_fisher", mode = "function")) {
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f, local = FALSE)
}

APP_VERSION <- "1.2"
APP_REPO <- "https://github.com/DunderChief/csea-tsea"
GOATCOUNTER_CODE <- ""   # e.g. "csea-tsea" once a goatcounter.com site exists; leave empty for no analytics

data_cache <- new.env()
get_panel_data <- function(panel) {
  force(panel)   # evaluate now: the server loop below reuses the variable name
  function() {
    if (is.null(data_cache[[panel$id]])) data_cache[[panel$id]] <- readRDS(file.path("data", panel$file))
    data_cache[[panel$id]]
  }
}
EXAMPLES <- readRDS("data/examples.rds")

about_md <- paste(readLines("about.md", warn = FALSE), collapse = "\n")
about_md <- gsub("{{version}}", APP_VERSION, about_md, fixed = TRUE)
about_md <- gsub("{{repo}}", APP_REPO, about_md, fixed = TRUE)
about_parts <- strsplit(about_md, "<!-- TABLE1 -->", fixed = TRUE)[[1]]
table1_html <- tags$table(
  class = "table table-sm table-striped w-auto", id = "table1",
  tags$caption(class = "caption-top", tags$strong("Table 1"), " of Xu et al. 2014 (J Neurosci 34:1420-31), as linked from the original CSEA page. ",
               # footnote and references, verbatim
               "*This line was generated to target Layer 5a Corticostriatal interneurons, but gene expression analysis revealed",
               "substantial contamination with lymphoid cell data. Additional Information about these mouse lines and cell types",
               "can be found in the following references: Doyle, J.P. et al. Application of a translational profiling approach",
               "for the comparative analysis of CNS cell types. Cell 135, 749-62 (2008). Dougherty, J.D. et al. Candidate Pathways",
               "for Promoting Differentiation or Quiescence of Oligodendrocyte Progenitor-like Cells in Glioma. Cancer research 72,",
               "4856-68 (2012). ",
               tags$em("Note for this version: Cb.Fthfd and Ctx.Fthfd are Cb.Aldh1L1 and Ctx.Aldh1L1 in the data (same driver gene, renamed).")),
  tags$thead(tags$tr(tags$th("Cell types and regions assayed"), tags$th("Abbreviations used"))),
  tags$tbody(lapply(seq_len(nrow(XU2014_TABLE1)), function(i)
    tags$tr(tags$td(XU2014_TABLE1$cell[i]), tags$td(tags$code(XU2014_TABLE1$code[i]))))))

ui <- page_navbar(
  id = "nav",
  # an ordinary flowing page: with a fillable page the About tab's text overflowed
  # the viewport-high content area and the footer was drawn across it
  fillable = FALSE,
  title = tags$span(tags$strong("CSEA / TSEA"), tags$small(class = "text-muted ms-2", "Specific Expression Analysis")),
  window_title = "CSEA / TSEA: Specific Expression Analysis",
  theme = bs_theme(version = 5, preset = "shiny", primary = "#8B0000"),
  header = tags$head(
    tags$style(HTML("
      .sea-table td, .sea-table th { font-size: 0.9rem; }
      .sea-cell:hover { outline: 2px solid #8B0000; }
      .navbar-brand small { font-weight: normal; font-size: 0.8rem; }
      /* the page scrolls as a whole (fillable = FALSE); keep the form in view while reading a long table */
      .bslib-sidebar-layout > .sidebar { position: sticky; top: 0; align-self: start; max-height: 100vh; overflow-y: auto; }
      .selectize-dropdown { max-width: calc(100vw - 2rem); }
      .selectize-dropdown .option { white-space: normal; }
    ")),
    tags$script(HTML("
      // switch to the About tab and scroll to an anchor (used by in-text links such as 'Table 1')
      function seaGoToAbout(anchor) {
        var link = document.querySelector('.navbar a.nav-link[data-value=\"about\"]');
        if (link && window.bootstrap) { bootstrap.Tab.getOrCreateInstance(link).show(); }
        setTimeout(function () { var el = document.getElementById(anchor); if (el) el.scrollIntoView({behavior: 'smooth'}); }, 150);
        return false;
      }
      // Shinylive runs the app inside an iframe, so a plain link would load its target in that frame
      // (github.com and the journal sites refuse to be framed). Open every external link, and the
      // sample PDF, in a new tab instead; in-page anchors are left alone.
      document.addEventListener('click', function (e) {
        var a = e.target.closest ? e.target.closest('a[href]') : null;
        if (!a || a.getAttribute('target')) return;
        var href = a.getAttribute('href');
        if (/^(https?:)?\\/\\//i.test(href) || /\\.pdf$/i.test(href)) {
          a.setAttribute('target', '_blank');
          a.setAttribute('rel', 'noopener');
        }
      });
    ")),
    if (nzchar(GOATCOUNTER_CODE)) tags$script(
      async = NA, `data-goatcounter` = sprintf("https://%s.goatcounter.com/count", GOATCOUNTER_CODE),
      src = "//gc.zgo.at/count.js")
  ),
  nav_panel(PANELS$csea$nav,  value = "csea",  analysis_ui("csea",  PANELS$csea)),
  nav_panel(PANELS$adult$nav, value = "adult", analysis_ui("adult", PANELS$adult)),
  nav_panel(PANELS$dev$nav,   value = "dev",   analysis_ui("dev",   PANELS$dev)),
  nav_panel(PANELS$tsea$nav,  value = "tsea",  analysis_ui("tsea",  PANELS$tsea)),
  nav_panel("About", value = "about",
            div(class = "container-md py-3", style = "max-width: 56rem;", id = "about-top",
                markdown(about_parts[1]), table1_html, markdown(about_parts[2]))),
  nav_spacer(),
  nav_item(tags$a(href = APP_REPO, target = "_blank", rel = "noopener", class = "nav-link", icon("github"), "Source")),
  footer = div(class = "text-center text-muted small py-3",
               sprintf("CSEA / TSEA v%s", APP_VERSION), " · Dougherty lab, Washington University in St. Louis · ",
               "computed entirely in your browser (R via WebAssembly)")
)

server <- function(input, output, session) {
  lapply(PANELS, function(p) analysis_server(p$id, p, get_panel_data(p), EXAMPLES))
}

shinyApp(ui, server)
