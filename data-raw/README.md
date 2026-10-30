# Raw data (not committed)

`prepare_data.R` builds the compressed tables in `app/data/` from three files. They are not stored in this
repository; the script looks for them in the directory given on the command line (or `$CSEA_RAW_DATA`, or
`data-raw/cache/`), and downloads any that are missing from the lab website into `data-raw/cache/` (git-ignored):

| File | Source | What the app uses |
|---|---|---|
| `pSI.data_1.0.tar.gz` (10 MB) | https://sites.wustl.edu/doughertylab/files/2014/01/pSI.data_1.0.tar_.gz | `mouse$psi.out`, `human$young.adulthood$psi.out`, `human$developmental.periods$psi.out` |
| `TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt` (1.9 MB) | https://sites.wustl.edu/doughertylab/files/2013/10/TableS3_NAR_Dougherty_Tissue_gene_pSI_v3-1.txt | the whole table (TSEA) |
| `pSI_1.1.tar.gz` (0.8 MB) | https://sites.wustl.edu/doughertylab/files/2014/01/pSI_1.1.tar_.gz | `candidate.genes` example lists |

`cell_type_labels.csv` maps the 35 mouse column codes to the cell-type names and regions printed on the original
CSEA bullseye (Xu et al. 2014, Table 1; `Cb.Fthfd`/`Ctx.Fthfd` there are `Cb.Aldh1L1`/`Ctx.Aldh1L1` in the data).
`in_plot` is FALSE for the eight whole-region columns, which the bullseye omits but the table keeps.

Output: `app/data/{csea_mouse,sea_adult,sea_development,tsea_tissue}.rds`, each a list with `psi` (numeric matrix,
upper-cased gene symbols x columns, NA = not specific), `columns` (code, label, group, ...), `source`, `built`;
plus `app/data/examples.rds`. About 1.5 MB in total with xz compression.

A copy of `TableS3_..._v3.txt` (the pre-fix version with Excel-mangled SEPT/MARCH/DEC gene names, otherwise
identical) is on the Internet Archive:
http://web.archive.org/web/20160322230737/http://genetics.wustl.edu/jdlab/files/2015/10/TableS3_NAR_Dougherty_Tissue_gene_pSI_v3.txt
