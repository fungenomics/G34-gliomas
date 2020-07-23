---
title: "03 - Co-expression of astrocyte and interneuron programs in normal brain"
date: "20 July, 2020"
output:
  html_document:
    keep_md: true
    theme: flatly
    css: ../../include/style.css
    toc: yes
    toc_depth: 4
    number_sections: true
    df_print: paged
    includes:
      before_body: ../../include/header.html
      after_body:  ../../include/footer.html
---

<!-- FRONT MATTER, insert configuration info -->


<!-- Load custom CSS/JS for code folding -->
<link rel="stylesheet" type="text/css" href="../../include/hideOutput.css">
<script src="../../include/hideOutput.js"></script>

***

# Configuration

Configuration of project directory & analysis outputs:

<details><summary>Show full config</summary>

```r
library(here)

# Set up outputs
message("Document index: ", doc_id)
```

```
## Document index: 03
```

```r
# Specify where to save outputs
out        <- here(subdir, "output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here(subdir, "figures", doc_id, "/"); dir.create(figout, recursive = TRUE)
cache      <- paste0("~/tmp/", basename(here()), "/", subdir, "/", doc_id, "/")
```

</details>

The root directory of this project is:

```
## /mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas-repo/singlecell_normal/output/03
```

```
## G34-gliomas-repo/singlecell_normal/figures/03//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Given that we a dual neuronal-astrocytic gene expression program in the tumors, here we interrogate the expression of interneuron and astrocytic
gene signatures in the normal brain, to understand if this is an
aberrant state, or also observed under normal conditions.


# Libraries


```r
# Load liraries here
library(here)

library(tidyverse)
library(Seurat)

source("../../../../misc/sjlib.R")
```


# Load data

Human fetal brain (5-37 PCW, Nowakowski et al, Science, 2017):


```r
load(here("singlecell_normal/data/Nowakowski2017_seurat.Rda"))
```

Prepare the palette for this dataset, so that astrocytes, internurons, and
radial glial cells are distinct:


```r
palette_nowakowski <- read_csv(here("bulk/data/Nowakowski2017_Table_4_cluster_labels_with_colours.csv")) %>%
  select(Cluster = `Cluster Name`, colour = Colour) %>%
  tibble::deframe()
```

```
## Parsed with column specification:
## cols(
##   `Cluster Name` = col_character(),
##   `Cluster Number (Fig 1B)` = col_double(),
##   `Cluster Interpretation` = col_character(),
##   Colour = col_character()
## )
```

```r
palette_nowakowski["Astrocyte"] <- "#FF3F38"
palette_nowakowski[grepl("RG", names(palette_nowakowski))] <- "#ffcc00"
```

Mouse developing forebrain (E12-P6)


```r
load(here("singlecell_normal/data/Jessa2019_forebrain_seurat.Rda"))
```

Prepare palette:


```r
palette_atlas <- read_tsv(here("bulk/data/Jessa2019_Table_2a.tsv")) %>% 
  select(Cluster, Colour) %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  tibble::deframe()
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   Sample = col_character(),
##   Age = col_character(),
##   Species = col_character(),
##   Structure = col_character(),
##   Cell_type = col_character(),
##   Cluster = col_character(),
##   Cell_class = col_character(),
##   Alias = col_character(),
##   Colour = col_character(),
##   Signature = col_character(),
##   Signature_note = col_character(),
##   Refinement_pons_progenitors = col_logical(),
##   Refinement_forebrain_embryonal_progenitors = col_logical(),
##   Refinement_forebrain_P0_progenitors = col_logical(),
##   Refinement_hg19w_astrocytes = col_logical(),
##   Pons_progenitors_trajectory = col_logical()
## )
```

```
## See spec(...) for full column specifications.
```

```r
palette_atlas[grepl("ASTR",      names(palette_atlas))] <- "#FF3F38"
palette_atlas[grepl("VRGC",      names(palette_atlas))] <- "#ffcc00"
palette_atlas[grepl("MGIN|CINH", names(palette_atlas))] <- "#084081"
```



# Analysis

## Get gene signatures

Here, we'll get the leading edge genes from the GSEA analysis, for one human
interneuron signature and one astrocyte signature, both enriched in G34R/V tumors.


```r
dge_le <- read_tsv(here("bulk/output/02/fgsea_results_padj<0.01_IDH_top_le.tsv"))
```

```
## Parsed with column specification:
## cols(
##   Comparison = col_character(),
##   Signature = col_character(),
##   Sample = col_character(),
##   Age = col_character(),
##   Cell_type = col_character(),
##   pval = col_double(),
##   padj = col_double(),
##   ES = col_double(),
##   NES = col_double(),
##   nMoreExtreme = col_double(),
##   size = col_double(),
##   leadingEdge = col_character(),
##   Dataset = col_character(),
##   leadingEdge_top = col_character()
## )
```

We load the genes which are detected in the tumors, in order to be able to restrict
our signatures to genes detected in the tumor scRNAseq data:


```r
load(here("singlecell_normal/data/G34_malignant_scRNAseq_genes.Rda"))
```

Intersect the signatures with genes which are detected in the single-cell tumor data:


```r
signatures <- readRDS(here("bulk/output/02/signatures_sym.Rds"))

# Newborn interneuron signature
interneuron_genes <- base::intersect(signatures$`HF nIN1`,      tumor_genes)

# Astrocyte signature
astro_genes       <- base::intersect(signatures$`HF Astrocyte`, tumor_genes)
```



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-07-20 17:13:08
```

The git repository and last commit:

```
## Local:    master /mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [2bc3284] 2020-07-09: Update README
```

The random seed was set with `set.seed(100)`

The R session info:
<details>

```
## R version 3.5.0 (2018-04-23)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /var/chroots/hydraex-centos-7/usr/lib64/R/lib/libRblas.so
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] stats     graphics  grDevices datasets  utils     methods   base     
## 
## other attached packages:
##  [1] bindrcpp_0.2.2  Seurat_2.3.4    Matrix_1.2-14   cowplot_0.9.4  
##  [5] forcats_0.3.0   stringr_1.3.1   dplyr_0.7.7     purrr_0.3.4    
##  [9] readr_1.3.1     tidyr_0.8.2     tibble_3.0.1    ggplot2_3.1.0  
## [13] tidyverse_1.2.1 here_0.1       
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.2.0        snow_0.4-3          backports_1.1.3    
##   [4] Hmisc_4.2-0         plyr_1.8.4          igraph_1.2.2       
##   [7] lazyeval_0.2.1      splines_3.5.0       digest_0.6.16      
##  [10] foreach_1.4.4       htmltools_0.3.6     lars_1.2           
##  [13] gdata_2.18.0        magrittr_1.5        checkmate_1.9.1    
##  [16] cluster_2.0.7-1     mixtools_1.1.0      ROCR_1.0-7         
##  [19] modelr_0.1.3        R.utils_2.7.0       colorspace_1.4-0   
##  [22] rvest_0.3.2         haven_2.0.0         xfun_0.12          
##  [25] crayon_1.3.4        jsonlite_1.6        bindr_0.1.1        
##  [28] survival_2.41-3     zoo_1.8-4           iterators_1.0.10   
##  [31] ape_5.2             glue_1.4.1          gtable_0.2.0       
##  [34] kernlab_0.9-27      prabclus_2.2-7      DEoptimR_1.0-8     
##  [37] scales_1.0.0        mvtnorm_1.0-10      bibtex_0.4.2       
##  [40] Rcpp_1.0.0          metap_1.1           dtw_1.20-1         
##  [43] htmlTable_1.13.1    reticulate_1.10     foreign_0.8-70     
##  [46] bit_1.1-14          proxy_0.4-22        mclust_5.4.2       
##  [49] SDMTools_1.1-221    Formula_1.2-3       tsne_0.1-3         
##  [52] stats4_3.5.0        htmlwidgets_1.3     httr_1.4.0         
##  [55] gplots_3.0.1.1      RColorBrewer_1.1-2  fpc_2.1-11.1       
##  [58] acepack_1.4.1       modeltools_0.2-22   ellipsis_0.2.0.1   
##  [61] ica_1.0-2           pkgconfig_2.0.2     R.methodsS3_1.7.1  
##  [64] flexmix_2.3-14      nnet_7.3-12         reshape2_1.4.3     
##  [67] tidyselect_1.1.0    rlang_0.4.6         munsell_0.5.0      
##  [70] cellranger_1.1.0    tools_3.5.0         cli_1.0.1          
##  [73] generics_0.0.2      broom_0.5.1         ggridges_0.5.1     
##  [76] evaluate_0.12       yaml_2.2.0          npsurv_0.4-0       
##  [79] knitr_1.21          bit64_0.9-7         fitdistrplus_1.0-14
##  [82] robustbase_0.93-2   caTools_1.17.1.1    RANN_2.6           
##  [85] pbapply_1.4-0       nlme_3.1-137        R.oo_1.22.0        
##  [88] xml2_1.2.0          hdf5r_1.0.0         compiler_3.5.0     
##  [91] rstudioapi_0.9.0    png_0.1-7           lsei_1.2-0         
##  [94] stringi_1.2.4       lattice_0.20-35     trimcluster_0.1-2.1
##  [97] vctrs_0.3.1         pillar_1.4.4        lifecycle_0.2.0    
## [100] BiocManager_1.30.10 Rdpack_0.10-1       lmtest_0.9-36      
## [103] data.table_1.12.0   bitops_1.0-6        irlba_2.3.3        
## [106] gbRd_0.4-11         R6_2.3.0            latticeExtra_0.6-28
## [109] renv_0.10.0         KernSmooth_2.23-15  gridExtra_2.3      
## [112] codetools_0.2-15    MASS_7.3-49         gtools_3.8.1       
## [115] assertthat_0.2.0    rprojroot_1.3-2     withr_2.1.2        
## [118] diptest_0.75-7      parallel_3.5.0      doSNOW_1.0.16      
## [121] hms_0.4.2           grid_3.5.0          rpart_4.1-13       
## [124] class_7.3-14        rmarkdown_1.11      segmented_0.5-3.0  
## [127] Rtsne_0.15          git2r_0.27.1        lubridate_1.7.4    
## [130] base64enc_0.1-3
```

</details>


***

<!-- END OF END MATTER -->
