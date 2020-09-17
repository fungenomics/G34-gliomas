
# G34-gliomas

Selin Jessa
selin.jessa@mail.mcgill.ca

This repository contains the code & data for the bulk analysis included
the G34R/V HGG manuscript (Chen, Deshmukh, Jessa, Hadjadj, et al),
for the analysis that was performed by our lab.

Contents:

* [Directory organization](https://github.com/fungenomics/G34-gliomas#directory-organization)
* [Notes for reproducibility](https://github.com/fungenomics/G34-gliomas#notes-for-reproducibility)
* [GitHub / version control](https://github.com/fungenomics/G34-gliomas#github--version-control)
* [Paper analysis not included in this repo](https://github.com/fungenomics/G34-gliomas#paper-analysis-not-included-in-this-repository)


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




## Notes for reproducibility


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
* `HTML` files, because these can be large, but they can be regenerated from the intermediate `.md` files saved
* Figures in `pdf` format, and figure source data
* Raw data and large analysis output / processed data files
* The actual packages in the R library 

An HTML file can be regenerated from a markdown `.md` file, like so on
the command line with the `rmarkdown` package:

`R --no-save -e "rmarkdown::render(input = 'bulk_transcriptome_epigenome/analysis/01-bulk_RNAseq_pipeline.md', output_format = 'html_document')"`


## Paper analysis not included in this repository

* In-house bulk RNAseq pipeline. The standard bulk RNA-seq analysis performed by the Kleinman Lab in-house pipeline
is run in the standard location on Beluga, with associated level 3 analysis at
`/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples`.
That directory contains a number of iterations as samples were added, etc;
the R Markdown files contain the paths for the exact outputs used in the final manuscript.
* scRNAseq analysis of the patient tumor samples and external datasets. For the processing and QC of individual
samples, the backup is at `/lustre03/project/6004736/sjessa/from_hydra/HGG-G34/samples`. The analysis performed by Veronique Lisi is
stored at `/lustre03/project/6004736/vlisi/fromHydra/SCRATCH/vlisi/LEGACY/G34`.
* Any analysis performed by the Jabado Lab (Carol Chen, Djihad Hadjadj, Shriya Deshmukh)


