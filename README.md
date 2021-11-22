
# G34-gliomas

This repository contains the code & data for the bulk analysis included
the G34R/V HGG manuscript ([Chen, Deshmukh, Jessa, Hadjadj, et al, Cell, 2020](https://doi.org/10.1016/j.cell.2020.11.012)),
for the analysis that was performed by our lab.

Contents:

* [Directory organization](https://github.com/fungenomics/G34-gliomas#directory-organization)
* [Notes for reproducibility](https://github.com/fungenomics/G34-gliomas#notes-for-reproducibility)
* [GitHub / version control](https://github.com/fungenomics/G34-gliomas#github--version-control)


## Directory organization

* `reference_datasets`: data, code, and processed outputs for scRNAseq
normal brain reference datasets, with one sub-directory per publication.
For external publications, the process
of obtaining or deriving cluster markers used in later analysis is recorded here.
* `bulk_transcriptome_epigenome`: data, code, figures and output for 
bulk RNA-seq and ChIP-seq analysis of patient tumors and tumor-derived cell lines.
* `singlecell_normal`: data, code, figures and output for analysis of the scRNAseq
normal brain data
* `renv`: directory maintained by the R package `renv`, containing the isolated
project specific library
* `include`: shared templates, Rmd/HTML headers/footers, and R functions used
throughout the analysis

## Map from figures to code

This section contains a pointer from each figure in the paper to the section (§) where it's generated in the code.
For each figure panel, I provide a partial path to the RMD/MD/HTML files within this repository/directory,
and then the section in the rendered HTML which specifically produces that panel.
As described below in the section on preproducibility, the source data for the figure is
typically saved alongside the figure itself.

### Figure 2

- **Figure 2A**: `bulk_transcriptome_epigenome/02-GSEA...`, § 6.1.1 Forebrain reference
- **Figure 2C**: `bulk_transcriptome_epigenome/02-GSEA...`, § 6.1.2 Striatal SVZ
- **Figure 2D**: `bulk_transcriptome_epigenome/01-bulk_RNAseq_pipeline...` § 4.3.2 Lineage specific TFs
- **Figure 2E**: `bulk_transcriptome_epigenome/03-ChIPseq...`, § 4.2.1 DGE
- **Figure 2F**: `singlecell_normal/analysis/01-interneuron_pseudotime...`,
    - top panel § 4.1 Cell type density along normal interneuron differentiation trajectory
    - bottom panel § 4.2 Plot genes of interest along pseudotime

### Figure S2

- **Figure S2A**: `bulk_transcriptome_epigenome/02-GSEA...`, § 6.2 GSEA enrichment plots
- **Figure S2B**: `bulk_transcriptome_epigenome/02-GSEA...`, § 6.4 Confirmation of signal by direct expression of gene programs
- **Figure S2D**: `bulk_transcriptome_epigenome/02-GSEA...`, § 6.1.3 Adult V-SVZ
- **Figure S2E**: `bulk_transcriptome_epigenome/01-bulk_RNAseq_pipeline...` § 4.3.2 Lineage specific TFs

### Figure 3

- **Figure 3A**: `bulk_transcriptome_epigenome/analysis/04-isogenic_cell_lines...`,
    - top panel § 4.4.3 Targeted DGE, for stem condition - GBM002
    - bottom panel § 4.4.2 Targeted DGE, for serum condition - GBM002
- **Figure 3C**: `bulk_transcriptome_epigenome/analysis/04-isogenic_cell_lines...` § 4.5.4 Visualize results
- **Figure 3D**: `singlecell_normal/analysis/02-gene_bubbleplots...`,
    - left § 4.1 Mouse developing forebrain
    - right § 4.3 Striatal SVZ

### Figure S5

- **Figure S5A**: `bulk_transcriptome_epigenome/01-bulk_RNAseq_pipeline...` § 4.4.1 G34 mutants
- **Figure S5B**: `singlecell_normal/analysis/02-gene_bubbleplots...`, § 4.2 Human



### Figure 5

- **Figure 5D**: `singlecell_normal/analysis/03-astrocyte_interneuron_coexpression...`, § 4.4 Developing mouse forebrain and 4.5 H3G34R/V tumors


### Figure S6

- **Figure S6C**: `singlecell_normal/analysis/03-astrocyte_interneuron_coexpression...`, § 4.3 Human fetal telencephalon


## Notes for reproducibility


### `rr` template & helpers

This repository uses the [`rr`](https://github.com/sjessa/rr) template, which contains
a set of R markdown templates to help me ensure reproducibility. Secondly, this also
provides a set of helper functions (located in `rr_helpers.R` and prefixed by `rr_` in the
function name) to help encourage documentation.


### R and R package versions

The R library for this project is managed with the package `renv`,
which:

1. maintains an isolated project-specific library in the `renv` folder,
2. stores packages according to version
3. records the R, Bioconductor, and package versions in the file `renv.lock`, which
can be used to reproduce the R package environment elsewhere

The R version used is 3.5.1.


### R Markdown

Each markdown/HTML file has a "Reproducibility report" at the bottom, indicating
when the document was last rendered, the most recent git commit when it was rendered,
the seed, and the R session info.


### Figures and source data

For most figures, the source data underlying the plot is saved along side the figure
in the respective `figures` directory. If so, a message is displayed
in the markdown/HTML files underneath the chunk whih produces the plot,
giving the path for the figures/source data within this project directory.

e.g. `[figure/source data @ G34-gliomas/bulk_transcriptome_epigenome/figures/01/gsx2_pdgfra_correlation…]`


### Analysis outputs

For most text file & R object outputs, there is a text file saved next to the object
with the extension `.desc`, with a very brief one-line description of what's contained in the file.

e.g. for the output file `bulk_transcriptome_epigenome/output/02/fgsea_df.tsv`,
there is an associated description file `bulk_transcriptome_epigenome/output/02/fgsea_df.desc`



## GitHub / version control

This directory is tracked with git and has an associated GitHub repository in the Kleinman
lab account at https://github.com/fungenomics/G34-gliomas.

The following are tracked / available on GitHub:

* `.Rmd` files, containing the code, and `.md` files, containing code and outputs
* Figures in `png` format
* Certain output files (`tsv`/`Rda`/`Rds`), if they're small
* The brief `desc` files for otputs
* The lockfile produced by the `renv` package

The following are not tracked / available on GitHub:

* Figures in `pdf` format, and figure source data
* Raw data and large analysis output / processed data files
* The actual packages in the R library 

## Contact

Selin Jessa
selin.jessa at mail.mcgill.ca
