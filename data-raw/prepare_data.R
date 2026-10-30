# Build the compressed pSI tables the app ships with (app/data/*.rds).
#
# Inputs (raw, not committed; see data-raw/README.md for provenance):
#   pSI.data_1.0.tar.gz   mouse bacTRAP cell types + human BrainSpan tables
#   TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt   GTEx tissue pSI table (Wells 2015)
#   pSI_1.1.tar.gz        only for the bundled candidate.genes example lists
#
# Usage (from the repo root):
#   Rscript data-raw/prepare_data.R [input_dir]
# input_dir is a directory holding the three files below; it defaults to
# $CSEA_RAW_DATA, then to data-raw/cache (git-ignored). Files that are missing
# there are downloaded from the lab website into data-raw/cache.

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[[1]] else Sys.getenv("CSEA_RAW_DATA", unset = "data-raw/cache")

files <- c(
  psi_data = "pSI.data_1.0.tar.gz",
  psi_pkg  = "pSI_1.1.tar.gz",
  table_s3 = "TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt"
)
urls <- c(
  psi_data = "https://sites.wustl.edu/doughertylab/files/2014/01/pSI.data_1.0.tar_.gz",
  psi_pkg  = "https://sites.wustl.edu/doughertylab/files/2014/01/pSI_1.1.tar_.gz",
  table_s3 = "https://sites.wustl.edu/doughertylab/files/2013/10/TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt"
)
paths <- file.path(input_dir, files)
names(paths) <- names(files)
for (k in names(files)) {
  if (!file.exists(paths[[k]])) {
    dir.create("data-raw/cache", showWarnings = FALSE, recursive = TRUE)
    paths[[k]] <- file.path("data-raw/cache", files[[k]])
    if (!file.exists(paths[[k]])) {
      message("Downloading ", urls[[k]])
      utils::download.file(urls[[k]], paths[[k]], mode = "wb", quiet = TRUE)
    }
  }
}

out_dir <- "app/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- cross-species backgrounds (pSI::fisher.iteration, background = "mouse.human" / "human.mouse") ----
# When the candidate list comes from the other species, pSI restricts the gene
# universe of the Fisher test to the genes in the table that have a clear
# homolog, defined as "symbol present in the UCSC kgXref table of the other
# assembly" (hg19 for a mouse table, mm10 for a human table). The old CSEA web
# tool used 13,421 for the mouse table with a human list: that is the constant
# left in pSI 1.1's source and it reproduces the saved output pages exactly;
# today's hg19 kgXref gives a smaller number, so the historic value is kept
# for the CSEA panel and the current tables are used where no historic value
# exists (human tables with a mouse list).
cache_dir <- "data-raw/cache"; dir.create(cache_dir, showWarnings = FALSE)
kgxref_symbols <- function(assembly) {
  f <- file.path(cache_dir, paste0("kgXref_", assembly, ".txt.gz"))
  if (!file.exists(f)) {
    url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/", assembly, "/database/kgXref.txt.gz")
    message("Downloading ", url); utils::download.file(url, f, mode = "wb", quiet = TRUE)
  }
  k <- utils::read.table(gzfile(f), sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  unique(toupper(k[, 5]))
}
hg19_symbols <- kgxref_symbols("hg19")
mm10_symbols <- kgxref_symbols("mm10")
kgxref_note <- sprintf("UCSC kgXref tables downloaded %s", format(Sys.Date()))
CSEA_MOUSE_HUMAN_TOTAL <- 13421L

# ---- pSI.data: mouse cell types, human adult regions, human development ----
tmp <- tempfile("psidata_"); dir.create(tmp)
utils::untar(paths[["psi_data"]], exdir = tmp)
e <- new.env()
load(file.path(tmp, "pSI.data", "data", "mouse.rda"), envir = e)
load(file.path(tmp, "pSI.data", "data", "human.rda"), envir = e)

as_psi_matrix <- function(df) {
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  rn <- toupper(rownames(df))
  if (anyDuplicated(rn)) stop("duplicated gene symbols after upper-casing")
  rownames(m) <- rn
  m
}

labels <- utils::read.csv("data-raw/cell_type_labels.csv", stringsAsFactors = FALSE)
mouse <- as_psi_matrix(e$mouse$psi.out)
stopifnot(setequal(colnames(mouse), labels$code))
labels <- labels[match(colnames(mouse), labels$code), ]   # keep the data's column order

write_panel <- function(name, psi, columns, source, universe) {
  obj <- list(
    psi = psi,                       # genes (upper-case symbols) x columns, pSI values, NA = not specific
    columns = columns,               # data.frame: code, label, group (+ extra display columns)
    universe = universe,             # named list: gene-universe size per candidate-list species
    source = source,                 # provenance text for the About page
    built = format(Sys.Date())
  )
  f <- file.path(out_dir, paste0(name, ".rds"))
  saveRDS(obj, f, compress = "xz")
  message(sprintf("%-18s %5d genes x %2d columns  -> %s (%.2f MB)", name, nrow(psi), ncol(psi), f,
                  file.size(f) / 1e6))
}

write_panel(
  "csea_mouse", mouse,
  data.frame(code = labels$code, label = labels$cell_type, group = labels$region,
             in_plot = labels$in_plot, stringsAsFactors = FALSE),
  "pSI.data 1.0, mouse$psi.out: bacTRAP cell-type profiles as of October 2013 (Xu et al. 2014).",
  list(mouse = nrow(mouse), human = CSEA_MOUSE_HUMAN_TOTAL,
       note = sprintf("human = %d mouse genes with a human homolog, the value used by pSI 1.1 and the original web tool; the current hg19 kgXref would give %d",
                      CSEA_MOUSE_HUMAN_TOTAL, sum(rownames(mouse) %in% hg19_symbols)))
)

human_cols <- function(cn) {
  region <- sub("\\..*$", "", cn)
  period <- sub("^[^.]+\\.", "", cn)
  period_label <- c(
    "Early.Fetal" = "Early Fetal", "Early.Mid.Fetal" = "Early Mid Fetal",
    "Late.Mid.Fetal" = "Late Mid Fetal", "Late.Fetal" = "Late Fetal",
    "Neonatal...Early.Infancy" = "Neonatal / Early Infancy", "Late.Infancy" = "Late Infancy",
    "Early.Childhood" = "Early Childhood", "Middle...Late.Childhood" = "Middle / Late Childhood",
    "Adolescence" = "Adolescence", "Young.Adulthood" = "Young Adulthood"
  )
  stopifnot(all(period %in% names(period_label)))
  data.frame(code = cn, label = region, group = unname(period_label[period]),
             region = region, period = unname(period_label[period]),
             period_order = match(period, names(period_label)),
             stringsAsFactors = FALSE)
}

adult <- as_psi_matrix(e$human$young.adulthood$psi.out)
human_universe <- function(psi) list(
  human = nrow(psi), mouse = sum(rownames(psi) %in% mm10_symbols),
  note = sprintf("mouse = genes with a mouse homolog (symbol in mm10 kgXref; %s)", kgxref_note))
write_panel("sea_adult", adult, human_cols(colnames(adult)),
            "pSI.data 1.0, human$young.adulthood$psi.out: BrainSpan RNA-seq, young adulthood, six brain regions.",
            human_universe(adult))

dev <- as_psi_matrix(e$human$developmental.periods$psi.out)
write_panel("sea_development", dev, human_cols(colnames(dev)),
            "pSI.data 1.0, human$developmental.periods$psi.out: BrainSpan RNA-seq, six brain regions x ten developmental periods.",
            human_universe(dev))

# ---- TSEA: GTEx tissues (Wells et al. 2015, Table S3) ----
s3 <- utils::read.table(paths[["table_s3"]], header = TRUE, sep = "\t", row.names = 1,
                        check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
tissue <- as_psi_matrix(s3)
write_panel("tsea_tissue", tissue,
            data.frame(code = colnames(tissue), label = gsub(".", " ", colnames(tissue), fixed = TRUE),
                       group = gsub(".", " ", colnames(tissue), fixed = TRUE), stringsAsFactors = FALSE),
            "Table S3 of Wells et al. 2015 (TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt): GTEx pilot RNA-seq (2013-01-31 release) averaged to 25 tissues, RefSeq release 60 protein-coding genes.",
            human_universe(tissue))

# ---- example gene lists bundled with pSI 1.1 ----
tmp2 <- tempfile("psi11_"); dir.create(tmp2)
utils::untar(paths[["psi_pkg"]], exdir = tmp2)
e2 <- new.env()
load(file.path(tmp2, "pSI", "data", "candidate.genes.rda"), envir = e2)
examples <- list(
  list(id = "autdb", name = "AutDB autism candidate genes (328)", genes = e2$candidate.genes$AutDB),
  list(id = "rdnv", name = "Protein-disrupting rare de novo variants, ASD probands (122)", genes = e2$candidate.genes$protein.disrupting.rdnv),
  list(id = "silent", name = "Silent rare de novo variants, unaffected siblings (122)", genes = e2$candidate.genes$silent.rdnv),
  list(id = "hcrt", name = "Hypocretin neuron ablation genes (63)", genes = e2$candidate.genes$hcrt.genes),
  list(id = "retina", name = "RetNet retinopathy genes (120)", genes = e2$candidate.genes$retinopathy.genes)
)
saveRDS(examples, file.path(out_dir, "examples.rds"), compress = "xz")
message("examples.rds written (", length(examples), " lists)")
