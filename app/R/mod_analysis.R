# Shiny module: one analysis panel (gene list in, table + bullseye + overlaps out).

# downloadButton() without the HTML `download` attribute. Chromium does not
# route links that carry that attribute through service workers (Chromium
# bug 468227), and under Shinylive the download is served by a service
# worker, so in Chrome the click reached the static host instead and failed
# with "File wasn't available on site". Without the attribute the link opens a
# navigation the service worker handles; the Content-Disposition header still
# makes the browser save it. This is the workaround in Shiny's own
# downloadButton documentation; it is harmless for a normal Shiny server.
download_button <- function(...) {
  tag <- shiny::downloadButton(...)
  tag$attribs$download <- NULL
  tag
}

analysis_ui <- function(id, panel) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    fillable = FALSE,
    sidebar = bslib::sidebar(
      width = 360, open = TRUE,
      shiny::tags$p(shiny::tags$strong(panel$prompt)),
      shiny::textAreaInput(ns("genes"), NULL, rows = 12, width = "100%"),
      shiny::radioButtons(ns("species"), "Candidate gene list from", choices = c("Human" = "human", "Mouse" = "mouse"),
                          selected = "human", inline = TRUE),
      shiny::div(
        class = "d-flex gap-2",
        shiny::actionButton(ns("run"), "Run", class = "btn-primary", icon = shiny::icon("play")),
        shiny::actionButton(ns("clear"), "Clear", class = "btn-outline-secondary")
      ),
      # dropdown attached to <body>: the sticky sidebar scrolls its own overflow and would clip it
      shiny::selectizeInput(ns("example"), "Load an example list", choices = c("Choose an example..." = ""),
                            options = list(
                              dropdownParent = "body",
                              # Shiny's selectize gives the dropdown width: max-content, which overruns a
                              # narrow screen; cap it at the control's width (long labels wrap)
                              onDropdownOpen = I("function($dropdown) { $dropdown.css('max-width', this.$control.outerWidth() + 'px'); }")))
    ),
    shiny::tags$h4(class = "text-center mt-2", page_heading(panel)),
    shiny::tags$p(style = "text-align: justify;", page_intro(panel)),
    shiny::uiOutput(ns("results")),
    shiny::tags$hr(),
    page_notes(panel)
  )
}

analysis_server <- function(id, panel, get_data, examples) {
  force(panel); force(get_data)   # R evaluates arguments lazily; bind this panel now, not at first use
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    shiny::updateSelectizeInput(session, "example",
      choices = c("Choose an example..." = "",
                  stats::setNames(vapply(examples, `[[`, "", "id"), vapply(examples, `[[`, "", "name"))))

    shiny::observeEvent(input$example, {
      if (nzchar(input$example)) {
        ex <- examples[[match(input$example, vapply(examples, `[[`, "", "id"))]]
        shiny::updateTextAreaInput(session, "genes", value = paste(ex$genes, collapse = "\n"))
      }
    })
    shiny::observeEvent(input$clear, {
      shiny::updateTextAreaInput(session, "genes", value = "")
      shiny::updateSelectizeInput(session, "example", selected = "")
    })

    result <- shiny::eventReactive(input$run, {
      parsed <- parse_gene_list(input$genes)
      shiny::validate(shiny::need(length(parsed$genes) > 0, "Enter at least one gene symbol, then press Run."))
      dat <- get_data()
      total <- dat$universe[[input$species]]
      res <- sea_fisher(dat$psi, parsed$genes, total = total)
      t05 <- res$table[res$table$threshold == 0.05, ]
      n05 <- stats::setNames(t05$n_specific, t05$code)
      list(res = res, parsed = parsed, dat = dat, species = input$species,
           lay = layout_nodes(panel$layout, dat$columns, n05))
    })

    output$results <- shiny::renderUI({
      shiny::req(input$run)
      shiny::tagList(
        shiny::uiOutput(ns("summary")),
        bslib::card(
          bslib::card_header(shiny::tags$span("Bullseye plot", class = "fw-bold"),
                             download_button(ns("png"), "Download PNG", class = "btn-sm float-end")),
          shiny::imageOutput(ns("plot"), height = "auto"),
          shiny::tags$p(class = "small text-muted mb-0 mt-2", BULLSEYE_CAPTION)
        ),
        bslib::card(
          bslib::card_header(shiny::tags$span("Results table", class = "fw-bold"),
                             download_button(ns("csv"), "Download CSV", class = "btn-sm float-end")),
          shiny::tags$p(class = "mb-2",
                        shiny::tags$strong(panel$interpretation[1]), shiny::tags$br(),
                        panel$interpretation[2], shiny::tags$br(),
                        panel$interpretation[3], shiny::tags$br(),
                        panel$interpretation[4]),
          shiny::uiOutput(ns("table")),
          shiny::uiOutput(ns("cell"))
        ),
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel("Overlapped genes at every pSI threshold", shiny::uiOutput(ns("overlaps")))
        )
      )
    })

    output$summary <- shiny::renderUI({
      r <- result()
      p <- r$parsed; res <- r$res
      n_missing <- length(res$missing)
      shiny::tags$div(
        class = "alert alert-light border",
        # the original result page's opening line (its "inputed" corrected to "inputted")
        shiny::tags$p(class = "mb-1",
          sprintf("Of the inputted %d genes,", p$n_entered), shiny::tags$br(),
          sprintf("%d genes exist in the %s, %d genes do not%s", length(res$found), panel$dataset_phrase,
                  n_missing, if (n_missing) ", which are:" else ".")),
        if (n_missing) shiny::tags$p(class = "mb-1 font-monospace small", paste(res$missing, collapse = ", ")),
        if (p$n_duplicates) shiny::tags$p(class = "mb-0 small text-muted",
          sprintf("(%d distinct symbols after removing %d duplicate%s)", res$n_genes, p$n_duplicates,
                  if (p$n_duplicates == 1) "" else "s")),
        shiny::tags$p(class = "mb-0 small text-muted",
          sprintf("Background for the Fisher test: %s genes (%s list; %s).", format(res$total, big.mark = ","),
                  if (r$species == "human") "human" else "mouse",
                  if (res$total == nrow(r$dat$psi)) "all genes in the dataset" else "genes with a clear homolog in the other species"))
      )
    })

    # Columns drawn on the bullseye first, then the rest (the mouse table has
    # eight whole-region profiles that the plot has never shown). Within each
    # block the original dataset order is kept.
    plot_order <- function(cols, codes) {
      in_plot <- if (!is.null(cols$in_plot)) cols$in_plot[match(codes, cols$code)] else rep(TRUE, length(codes))
      in_plot[is.na(in_plot)] <- TRUE
      ord <- order(!in_plot, seq_along(codes))
      list(codes = codes[ord], in_plot = in_plot[ord])
    }
    BLOCK_LABELS <- c("Cell types (drawn on the bullseye)", "Whole regions (not on the bullseye)")

    # Readable name for a column: the cell type with its region, the region with
    # its period, or the tissue. Shown as "Readable (code)" in the table.
    row_name <- function(cols, code) {
      i <- match(code, cols$code)
      if (!is.null(cols$period)) paste(cols$region[i], cols$period[i], sep = ", ")
      else if (!is.na(cols$group[i]) && !grepl(tolower(cols$group[i]), tolower(cols$label[i]), fixed = TRUE)) paste(cols$label[i], cols$group[i], sep = ", ")
      else cols$label[i]
    }
    row_label <- function(cols, code) shiny::tagList(
      row_name(cols, code), " ", shiny::tags$code(class = "text-muted", paste0("(", code, ")")))

    output$table <- shiny::renderUI({
      r <- result()
      tab <- r$res$table; cols <- r$dat$columns
      fmt <- function(x) ifelse(x < 0.001 & x > 0, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "f"))
      ths <- SEA_THRESHOLDS
      po <- plot_order(cols, unique(tab$code))
      make_row <- function(code) {
        cells <- lapply(ths, function(th) {
          z <- tab[tab$code == code & tab$threshold == th, ]
          # same colour scale as the bullseye rings, lightened to sit behind text
          bg <- bullseye_colour(z$bh, lighten = 0.55)
          shiny::tags$td(
            class = "text-nowrap sea-cell",
            style = paste0("cursor: pointer;", if (!is.na(bg)) paste0(" background-color: ", bg, ";") else ""),
            title = sprintf("%d of %d specific genes overlap; click for the list", z$n_overlap, z$n_specific),
            onclick = sprintf("Shiny.setInputValue('%s', {code: '%s', threshold: %s}, {priority: 'event'})",
                              ns("cell_click"), code, format(th, scientific = FALSE)),
            sprintf("%s (%s)", fmt(z$p), fmt(z$bh)))
        })
        shiny::tags$tr(shiny::tags$td(row_label(cols, code)), cells)
      }
      block_row <- function(label) shiny::tags$tr(class = "table-light",
        shiny::tags$th(colspan = length(ths) + 1, class = "small text-uppercase text-muted fw-semibold", label))
      rows <- if (all(po$in_plot)) lapply(po$codes, make_row) else list(
        block_row(BLOCK_LABELS[1]), lapply(po$codes[po$in_plot], make_row),
        block_row(BLOCK_LABELS[2]), lapply(po$codes[!po$in_plot], make_row))
      shiny::tags$div(class = "table-responsive",
        shiny::tags$table(class = "table table-sm table-hover align-middle sea-table",
          shiny::tags$thead(shiny::tags$tr(
            shiny::tags$th(panel$table_corner),
            lapply(ths, function(th) shiny::tags$th(format(th, scientific = FALSE))))),
          shiny::tags$tbody(rows)))
    })

    output$cell <- shiny::renderUI({
      click <- input$cell_click
      r <- result()
      if (is.null(click)) return(shiny::tags$p(class = "small text-muted", "Click a p-value to see the overlapping genes."))
      z <- r$res$table[r$res$table$code == click$code & abs(r$res$table$threshold - click$threshold) < 1e-12, ]
      if (!nrow(z)) return(NULL)
      genes <- strsplit(z$overlap, ";")[[1]]
      shiny::tags$div(class = "alert alert-secondary mt-2 mb-0",
        shiny::tags$strong(sprintf("%s (%s) at pSI < %s: ", row_name(r$dat$columns, click$code), click$code, format(click$threshold, scientific = FALSE))),
        sprintf("%d of your genes among the %d specific genes (p = %s, BH = %s). ", z$n_overlap, z$n_specific, fmt3(z$p), fmt3(z$bh)),
        if (length(genes)) shiny::tags$span(class = "font-monospace", paste(genes, collapse = ", ")) else "No overlap.")
    })

    output$overlaps <- shiny::renderUI({
      r <- result(); tab <- r$res$table; cols <- r$dat$columns
      po <- plot_order(cols, unique(tab$code))
      dl <- function(z) shiny::tags$dl(class = "row small mb-0",
        lapply(seq_len(nrow(z)), function(i) list(
          shiny::tags$dt(class = "col-sm-4 fw-normal", row_label(cols, z$code[i])),
          shiny::tags$dd(class = "col-sm-8 font-monospace", if (nzchar(z$overlap[i])) gsub(";", ", ", z$overlap[i]) else "none"))))
      lapply(SEA_THRESHOLDS, function(th) {
        z <- tab[tab$threshold == th, ]
        z <- z[match(po$codes, z$code), ]
        shiny::tags$div(
          shiny::tags$h6(class = "mt-3", sprintf("The overlapped genes for every %s ( pSI threshold %s ) are:",
                                                 panel$overlap_noun, format(th, scientific = FALSE))),
          if (all(po$in_plot)) dl(z) else list(
            shiny::tags$p(class = "small text-uppercase text-muted fw-semibold mb-1", BLOCK_LABELS[1]), dl(z[po$in_plot, ]),
            shiny::tags$p(class = "small text-uppercase text-muted fw-semibold mb-1 mt-2", BLOCK_LABELS[2]), dl(z[!po$in_plot, ])))
      })
    })

    draw_plot <- function() {
      r <- result()
      bullseye_plot(r$res$table, r$lay, base_size = 13)
    }
    # The plot is drawn once per run as a fixed-size PNG and scaled by CSS.
    # renderPlot() would redraw it at the container width times the screen's
    # pixel ratio on every resize, which under WebAssembly means rasterising
    # and encoding several megapixels each time; in Shinylive that never
    # completed.
    plot_file <- shiny::reactive({
      f <- tempfile(fileext = ".png")
      w <- 1200
      grDevices::png(f, width = w, height = round(w * bullseye_aspect(result()$lay)), res = 100)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_plot()
      f
    })
    output$plot <- shiny::renderImage({
      list(src = plot_file(), contentType = "image/png", alt = "Bullseye plot",
           style = "max-width: 100%; height: auto;")
    }, deleteFile = FALSE)

    output$csv <- shiny::downloadHandler(
      filename = function() paste0(panel$id, "_results.csv"),
      content = function(file) {
        r <- result(); tab <- r$res$table; cols <- r$dat$columns
        po <- plot_order(cols, unique(tab$code))
        tab <- tab[order(match(tab$code, po$codes), -tab$threshold), ]
        i <- match(tab$code, cols$code)
        out <- data.frame(code = tab$code, name = cols$label[i], group = cols$group[i],
                          on_bullseye = po$in_plot[match(tab$code, po$codes)],
                          pSI_threshold = tab$threshold, specific_genes = tab$n_specific,
                          overlapping_genes = tab$n_overlap, p_value = tab$p, bh_adjusted_p = tab$bh,
                          overlap = tab$overlap, stringsAsFactors = FALSE)
        utils::write.csv(out, file, row.names = FALSE)
      })
    output$png <- shiny::downloadHandler(
      filename = function() paste0(panel$id, "_bullseye.png"),
      content = function(file) {
        grDevices::png(file, width = 12, height = 12 * bullseye_aspect(result()$lay), units = "in", res = 150, bg = "white")
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_plot()
      })
  })
}

fmt3 <- function(x) ifelse(x < 0.001 & x > 0, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "f"))
