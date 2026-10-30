# The four analyses.
#
# Every string marked VERBATIM below is the wording of the original pages
# (http://doughertytools.wustl.edu/CSEAtool.html and TSEAtool.html, and their
# result pages, as saved on 2026-10-22) reproduced unchanged, except that three
# typos were corrected at David's request ("desease", "18,056 rof", "inputed")
# and link targets that no longer exist were re-pointed; those changes are
# listed on the About page. Text written for this version lives in about.md.

# Links as HTML strings, and paragraphs assembled with para(), so that no
# whitespace is inserted around inline links (the original punctuation,
# "(Table 1)", "(Dougherty et al.)", must sit tight against the link text).
ext_link <- function(text, href) sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', href, text)
about_link <- function(text, anchor) sprintf('<a href="#" onclick="return seaGoToAbout(\'%s\');">%s</a>', anchor, text)
para <- function(..., class = NULL) shiny::tags$p(class = class, shiny::HTML(paste0(...)))

# ---- "Specific Expression Analysis (SEA)" page: CSEA, adult regions, development ----

SEA_HEADING <- "Specific Expression Analysis (SEA)"   # VERBATIM

SEA_INTRO <- paste(  # VERBATIM
  "Recent advances have substantially increased the number of genes associated with complex genetic disorders",
  "of the CNS such as autism and schizophrenia. It is now clear that there are hundreds of distinct loci",
  "contributing to these disorders, underscoring a remarkable genetic heterogeneity. However, the commonality",
  "of symptoms across patients indicates there must be a functional convergence downstream of these loci upon",
  "the circuits mediating the affected behaviors. One possible mechanism for this convergence would be the",
  "selective expression of at least a subset of these genes in the cell types that comprise the circuit. Here,",
  "we have provided a simple web server implementing our approach of Cell-type Specific Expression Analysis",
  "(CSEA), for identifying candidate cell populations likely to be disrupted across a set of patients with",
  "diverse genetic lesions. This leverages a mouse transcriptomic profiling dataset to identify sets of",
  "transcripts specifically expressed in particular mouse cells. Though data from specific cell types in the",
  "brain is not available in humans, we have also analyzed human data from the Brainspan collection to provide",
  "a parallel analysis of transcripts enriched in particular human brain regions, and/or developmental windows.")

sea_notes <- function() shiny::tagList(  # VERBATIM paragraphs; links as on the original page
  para(
    "This server can accept as an input lists of gene symbols (as assigned by HUGO), and returns an enrichment ",
    "analysis of their expression across the currently available cell types and regions (",
    about_link("Table 1", "table1"),
    "). Candidate genes lists that overlap with lists of transcripts enriched in a particular cell type or region ",
    "are identified by Fisher's exact test with Benjamini-Hochberg correction."),
  para(
    "User input candidate gene lists may be from any source including, candidate gene lists from genetic ",
    "analyses of psychiatric disorder, lists of genes with altered expression from RNA-seq or Microarray analysis ",
    "from postmortem patient tissue, animal models, tumor samples, etc. To properly calculate expected overlap ",
    "for Fisher's Test, user should indicate whether candidate list is from human or mouse data."),
  para(
    "For each cell type (mouse) or brain region (human), transcripts that are specifically expressed or ",
    "enriched have been identified as described (",
    ext_link("Dougherty et al", "https://pubmed.ncbi.nlm.nih.gov/20308160/"),
    ".) at Specificity Index thresholds (pSI) of varying stringency (e.g., pSI smaller than .01 identifies a ",
    "larger number of relatively enriched transcripts, while those at pSI smaller than .0001 will be just the ",
    "subset that is relatively specific). Results will appear for all thresholds. Sample results and explanation ",
    "are ", ext_link("here", "CSEA_Sample_File.pdf"), "."),
  para(
    "Human data are from BrainSpan: Atlas of the Developing Human Brain [Internet]. Funded by ARRA Awards ",
    "1RC2MH089921-01, 1RC2MH090047-01, and 1RC2MH089929-01. © 2011. Available from ",
    ext_link("here", "https://www.brainspan.org/"), "."),
  para(class = "small text-muted mb-1",
    "Version 1.1: Updated 7/31/17 to address a bug fix in which layer 6 neurons and layer 5b neurons were ",
    "swapped in the CSEA plot output.")
)

# ---- "Tissue Specific Expression Analysis (TSEA)" page ----

TSEA_HEADING <- "Tissue Specific Expression Analysis (TSEA)"   # VERBATIM

TSEA_INTRO <- paste(  # VERBATIM
  "Every cell in the human body harbors a nearly identical copy of the genome, yet each tissue only expresses",
  "the subset of genes that are required to fulfill its particular functions. Human diseases are often quite",
  "clearly caused by dysfunction or disruption of particular tissues, yet human disease genetics typically",
  "treat all regions of the genome equally in regard to their prior probability of contributing to disease",
  "risk. If it could be established that genes with highly enriched or specific expression in a particular",
  "tissue were more likely to contribute to diseases of that tissue, then tissue specific expression could be",
  "reasonably used to weight tests of statistical association to improve our ability to identify risk loci.",
  "Here, we provide a simple web tool to test this 'selective expression' hypothesis by testing whether",
  "currently identified disease genes are over-represented by enriched expression in the disease-relevant",
  "tissue. We call our tool Tissue Specific Expression Analysis (TSEA).")

tsea_notes <- function() shiny::tagList(  # VERBATIM paragraphs; links as on the original page
  para(
    "This server can accept an input list of gene symbols (as assigned by HUGO), and return an enrichment ",
    "analysis of their expressions across 25 tissues. Candidate gene lists that overlap with lists of ",
    "transcripts enriched in a particular tissue are identified by Fisher's Exact test with Benjamini-Hochberg ",
    "correction."),
  para(
    "User gene lists may be from any source including, lists of genes from annotated GWAS studies, Gene ",
    "Ontological (GO) analyses, genetic analyses of disease/traits related tissues, lists of genes with altered ",
    "expression from RNA-seq or Microarray analysis from postmortem patient tissue, animal models, tumor ",
    "samples, etc."),
  para(
    "The gene expression data used in TSEA was collected using published RNA-Seq data from the Genotype-Tissue ",
    "Expression project (",
    ext_link("GTEx Analysis Pilot Data", "https://pubmed.ncbi.nlm.nih.gov/?term=The+Genotype-Tissue+Expression+(GTEx)+project"),
    " 2013-01-31, summarized to genes). The GTEx project was supported by the Common Fund of the Office of the ",
    "Director of the National Institutes of Health (commonfund.nih.gov/GTEx). The raw GTEx data consists of ",
    "1,839 samples derived from 189 post-mortem subjects. This included samples from 45 different tissues, with ",
    "some tissues offering multiple 'sub-tissue' types (e.g., multiple brain dissections). To analyze the data ",
    "at tissue level, RPKM values for the sub-tissue types were averaged resulting in 25 'whole-tissue' types. ",
    "Transcripts were further filtered to include well-annotated protein-coding genes designated by RefSeq ",
    "(release 60) gene annotations. After filtering, 18,056 of the original 52,576 transcripts remained."),
  para(
    "For each tissue, transcripts from the processed GTEx transcripts that are specifically expressed or ",
    "enriched have been identified by using our pSI R package function specificity.index (download ",
    ext_link("here", "https://github.com/DunderChief/pSI"),
    ") to calculate the Specificity Index thresholds (pSI) of varying stringency (",
    ext_link("Dougherty et al", "https://pubmed.ncbi.nlm.nih.gov/20308160/"), " and ",
    ext_link("CSEA publication", "https://pubmed.ncbi.nlm.nih.gov/24453331/"),
    "). For example, pSI smaller than .01 identifies a larger number of relatively enriched transcripts, while ",
    "those at pSI smaller than .0001 will be just the subset that is relatively specific. Results will appear ",
    "for all thresholds. Sample results and explanation are (",
    if (file.exists("www/TSEA_Sample_File.pdf")) ext_link("here", "TSEA_Sample_File.pdf") else "here",
    ")."),
)

# ---- Table 1 of Xu et al. 2014 (J Neurosci 34:1420), as linked from the original page; VERBATIM ----

XU2014_TABLE1 <- data.frame(
  cell = c("Cone cells", "Rod cells", "Hypothalamus", "Hypocretinergic neurons of hypothalamus",
           "Basal forebrain", "Cholinergic neurons of basal forebrain", "Brain stem",
           "Cholinergic neurons of brain stem", "Serotonergic neurons of brain stem", "Cerebellum",
           "Bergman glia and mature oligodendrocytes of cerebellum", "Purkinje cells of cerebellum",
           "Granule cells of cerebellum", "Stellate and basket cells of cerebellum",
           "Unipolar brush cells and Bergman glia of cerebellum", "Golgi neurons of cerebellum",
           "Astroglia of cerebellum", "Myelinating Oligodendrocyte of cerebellum",
           "Striatum(Caudate and Putamen)", "Drd2+ medium spiny neurons of striatum",
           "Drd1+ medium spiny neurons of striatum", "Cholinergic neurons of striatum", "Cortex",
           "Oligodendrocyte Progenitor Cells of cortex", "Immune Cells*", "Pnoc+ neurons of cortex",
           "Ntsr+ neurons of cortex", "Glt25d2 neurons of cortex", "Astrocytes of cortex",
           "Cort+ neurons of cortex", "Myelinating Oligodendrocyte of cortex", "Habenula (epithalamus)",
           "Cholinergic neurons of habenula", "Spinal cord", "Cholinergic neurons of spinal cord"),
  code = c("RetC.Cones", "RetR.Rods", "Hyp", "Hyp.Hcrt", "BF", "BF.Chat", "BS", "BS.Chat", "Bs.Slc6a4", "Cb",
           "Cb.Septin4", "Cb.Pcp2", "Cb.Neurod1", "Cb.Lypd6", "Cb.Grp", "Cb.Grm2", "Cb.Fthfd", "Cb.Cnp", "Cpu",
           "Cpu.D2", "Cpu.D1", "Cpu.Chat", "Ctx", "Ctx.Pdgfrajd340", "Ctx.Etv1_ts88", "Ctx.Pnoc", "Ctx.Ntsr",
           "Ctx.Glt25d2", "Ctx.Fthfd", "Ctx.Cort", "Ctx.Cnp", "Epi", "Epi.Chat", "Spc", "Spc.Chat"),
  stringsAsFactors = FALSE)

# ---- Sample-file annotation of the bullseye (CSEA_Sample_File.pdf); VERBATIM ----
BULLSEYE_CAPTION <- paste(
  "Graphical Version of the table (Bullseye plot), with B-H P-values plotted by color.",
  "Varying stringencies for enrichment (pSI): are represented by the size of the hexagons going from least",
  "specific lists (outer hexagons) to most specific lists (center). Hexagons scaled to size of gene lists.")

# ---- panel definitions ----
# prompt, dataset_phrase, table_corner, interpretation and overlap_noun are the
# original pages' wording (the result pages' interpretation text had two
# copy-paste slips, "brain cell types" on the development page and "adult brain
# regions and development" on the CSEA page; the nouns are corrected here).

PANELS <- list(
  csea = list(
    id = "csea", file = "csea_mouse.rds", layout = "csea", page = "sea",
    nav = "CSEA: brain cell types",
    prompt = "Enter candidate gene symbols for specific expression analysis across cell types (CSEA):",
    dataset_phrase = "brain cell type expression dataset",
    table_corner = "Brain Cell Types and P-Values",
    interpretation = c(
      "Interpretation of CSEA results table:",
      "The first row lists the pSI thresholds, and the first column lists the brain cell types.",
      "The values in the tables are the Fisher's Exact p-values, followed by the Benjamini-Hochberg (BH) corrected values in parentheses.",
      "Click the linked p-values to see the overlapped genes for the corresponding brain cell types and thresholds."),
    overlap_noun = "cell type"
  ),
  adult = list(
    id = "adult", file = "sea_adult.rds", layout = "row", page = "sea",
    nav = "SEA: adult brain regions",
    prompt = "Enter candidate gene symbols for SEA of adult brain regions:",
    dataset_phrase = "brain region and development expression dataset",
    table_corner = "Brain Regions and Development and P-Values",
    interpretation = c(
      "Interpretation of SEA brain region and development results table:",
      "The first row lists the pSI thresholds, and the first column lists the adult brain regions and development.",
      "The values in the tables are the Fisher's Exact p-values, followed by the Benjamini-Hochberg (BH) corrected values in parentheses.",
      "Click the linked p-values to see the overlapped genes for the corresponding adult brain regions and development thresholds."),
    overlap_noun = "brain region and development"
  ),
  dev = list(
    id = "dev", file = "sea_development.rds", layout = "grid", page = "sea",
    nav = "SEA: brain regions across development",
    prompt = "Enter candidate gene symbols for SEA across brain regions and development:",
    dataset_phrase = "brain region and development expression dataset",
    table_corner = "Brain Regions and Development and P-Values",
    interpretation = c(
      "Interpretation of SEA brain region and development results table:",
      "The first row lists the pSI thresholds, and the first column lists the brain regions and development.",
      "The values in the tables are the Fisher's Exact p-values, followed by the Benjamini-Hochberg (BH) corrected values in parentheses.",
      "Click the linked p-values to see the overlapped genes for the corresponding brain regions and development thresholds."),
    overlap_noun = "brain region and development"
  ),
  tsea = list(
    id = "tsea", file = "tsea_tissue.rds", layout = "tsea", page = "tsea",
    nav = "TSEA: tissues",
    prompt = "Enter candidate gene symbols for specific expression analysis of tissues:",
    dataset_phrase = "tissue expression dataset",
    table_corner = "Tissues and P-Values",
    interpretation = c(
      "Interpretation of TSEA results table:",
      "The first row lists the pSI thresholds, and the first column lists the tissues.",
      "The values in the tables are the Fisher's Exact p-values, followed by the Benjamini-Hochberg (BH) corrected values in parentheses.",
      "The overlapped genes for the corresponding tissues and thresholds are listed at the end."),
    overlap_noun = "tissue"
  )
)

page_heading <- function(panel) if (panel$page == "sea") SEA_HEADING else TSEA_HEADING
page_intro <- function(panel) if (panel$page == "sea") SEA_INTRO else TSEA_INTRO
page_notes <- function(panel) if (panel$page == "sea") sea_notes() else tsea_notes()
