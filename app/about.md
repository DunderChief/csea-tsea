## About

**CSEA / TSEA** version {{version}} is a rebuild of the Dougherty lab's Cell-type Specific Expression
Analysis (CSEA) and Tissue Specific Expression Analysis (TSEA) web tools, formerly served from
`doughertytools.wustl.edu`. The four analyses, their inputs, the statistics and the bullseye plots are the
same, and the text on the analysis pages is the original pages' text, reproduced verbatim. The difference is
that this version is a static web page: the R code runs in your browser through WebAssembly, so there is no
server, nothing to maintain, and your gene lists never leave your computer.

### What is computed

For each cell type, region or tissue, the genes with a specificity index p-value (pSI) below each of four
thresholds (0.05, 0.01, 0.001, 0.0001) form a "specific" gene list. Your candidate list is matched to the
dataset by upper-cased gene symbol (duplicates are counted once; symbols not in the dataset are reported and
left out). For every column and threshold a one-sided Fisher's exact test asks whether the overlap between
your list and the specific list is larger than expected by chance. The background of the test is every gene
in the dataset when the candidate list comes from the same species as the dataset, and the genes with a clear
homolog in the other species when it does not (the "Candidate gene list from" selector; the results summary
states the number used). P-values are Benjamini-Hochberg corrected across the columns, separately at each
threshold. This is exactly what `fisher.iteration()` in the
[pSI R package](https://github.com/DunderChief/pSI) computes.

### How to cite

If you use these tools, please cite the paper for the analysis you ran:

- **CSEA (brain cell types)**: Xu X, Wells AB, O'Brien DR, Nehorai A, Dougherty JD. Cell type-specific expression
  analysis to identify putative cellular mechanisms for neurogenetic disorders. *J Neurosci.* 2014;34(4):1420-31.
  [doi:10.1523/JNEUROSCI.4488-13.2014](https://doi.org/10.1523/JNEUROSCI.4488-13.2014), PMID 24453331.
- **SEA (human brain regions, adult and developmental)**: also Xu et al. 2014, above. Human data are from
  [BrainSpan: Atlas of the Developing Human Brain](https://www.brainspan.org/), funded by ARRA Awards
  1RC2MH089921-01, 1RC2MH090047-01 and 1RC2MH089929-01.
- **TSEA (tissues)**: Wells A, Kopp N, Xu X, O'Brien DR, Yang W, Nehorai A, Adair-Kirk TL, Kopan R, Dougherty JD.
  The anatomical distribution of genetic associations. *Nucleic Acids Res.* 2015;43(22):10804-20.
  [doi:10.1093/nar/gkv1262](https://doi.org/10.1093/nar/gkv1262), PMID 26586807.
- **The specificity index (pSI)**: Dougherty JD, Schmidt EF, Nakajima M, Heintz N. Analytical approaches to RNA
  profiling data for the identification of genes enriched in specific cells. *Nucleic Acids Res.*
  2010;38(13):4218-30. [doi:10.1093/nar/gkq130](https://doi.org/10.1093/nar/gkq130), PMID 20308160.

### Data and provenance

- **Brain cell types (CSEA)**: `pSI.data` 1.0, `mouse$psi.out`: 16,907 genes x 35 columns (27 bacTRAP
  cell types plus 8 whole-region profiles), the bacTRAP dataset as of October 2013 used by Xu et al. 2014.
  Column abbreviations are those of Table 1 of that paper (below); `Cb.Aldh1L1` and `Ctx.Aldh1L1` appear there
  as `Cb.Fthfd` and `Ctx.Fthfd` (same driver gene, renamed).
- **Adult brain regions and development (SEA)**: `pSI.data` 1.0, `human$young.adulthood$psi.out` (16,866 genes
  x 6 regions) and `human$developmental.periods$psi.out` (16,866 genes x 6 regions x 10 periods), derived from
  BrainSpan RNA-seq.
- **Tissues (TSEA)**: Table S3 of Wells et al. 2015 (`TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt`):
  18,056 RefSeq protein-coding genes x 25 tissues, from GTEx pilot RNA-seq (2013-01-31 release) with
  sub-tissues averaged to whole tissues. The GTEx project is supported by the Common Fund of the Office of the
  Director of the NIH.
- The example gene lists are those bundled with the pSI package; see `?candidate.genes` in the package for
  their sources.

The tables are shipped with this page as compressed R data files built by `data-raw/prepare_data.R` in the
[source repository]({{repo}}), where the raw inputs and their download URLs are documented.

<!-- TABLE1 -->

### Privacy

The analysis runs in your browser (R compiled to WebAssembly via [Shinylive](https://github.com/posit-dev/r-shinylive)).
Gene lists are never transmitted. The page is served as static files from GitHub Pages.

### Source code

Source, issue tracker and build instructions: [{{repo}}]({{repo}}).
