
# G34-gliomas

Selin Jessa
selin.jessa@mail.mcgill.ca

This repository contains the code & data for the bulk analysis included
the G34R/V HGG manuscript (Chen, Deshmukh, Jessa, Hadjadj, et al),
for the analysis that was performed by our lab.



## Directory organization

* `reference_datasets`: data, code, and processed outputs for scRNAseq
normal brain reference datasets, with one sub-directory per publication. The process
of obtaining or deriving cluster markers used in later analysis is recorded here.
* `bulk_transcriptome_epigenome`: data, code, figures and output for 
bulk RNA-seq and ChIP-seq analysis of patient tumors and tumor-derived cell lines.
* `singlecell_normal`: data, code, figures and output for 
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
in the respective `figures` directory. If so, a message like this is displayed
in the markdown/HTML files underneath the chunk whih produces the plot,
giving the path for the figures/source data within this project directory:

`[figure/source data @ G34-gliomas/bulk_transcriptome_epigenome/figures/01/gsx2_pdgfra_correlation…]`


### Analysis outputs






## GitHub repository

This directory is tracked with git and has an associated GitHub repository in the Kleinman
lab account at https://github.com/fungenomics/G34-gliomas.

The following are tracked / available on GitHub:

* `.Rmd` files, containing the code, and `.md` files, containing code and outputs
* Figures in `png` format
* Certain output files (`tsv`/`Rda`/`Rds`), if they're small
* The lockfile produced by the `renv` package

The following are not tracked / available on GitHub:
* `HTML` files, because these can be large, but they can be regenerated from the intermediate `.md` files saved
* Figures in `pdf` format
* Raw data and large processed data files
* The actual packages in the R library 


