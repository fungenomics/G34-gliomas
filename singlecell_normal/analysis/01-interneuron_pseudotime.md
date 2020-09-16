---
title: "01 - Interneuron pseudotemporal trajectory"
date: "16 September, 2020"
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
## Document index: 01
```

```r
# Specify where to save outputs
out        <- here(subdir, "output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here(subdir, "figures", doc_id, "/"); dir.create(figout, recursive = TRUE)
cache      <- paste0("~/tmp/", basename(here()), "/", subdir, "/", doc_id, "/")

message("Cache: ", cache)
```

```
## Cache: ~/tmp/G34-gliomas/singlecell_normal/01/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/singlecell_normal/output/01
```

```
## G34-gliomas/singlecell_normal/figures/01//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

To pinpoint the where along interneuron differentiation
G34 HGGs may arise, we retrieve genes which are (1) upregulated by RNA-seq and
H3K27-acetylated in G34 gliomas compared to other HGGs, and (2) downregulated
by RNA-seq and H3K27-trimethylated in G34 gliomas compared to other HGGs,
from the multi-omics analysis in the `bulk_transcriptome_epigenome` directory.
We then interrogate the expression of these genes along the normal interneuron
differentiation trajectory that we've constructed previously with monocle.


# Libraries


```r
# Load libraries here
library(tidyr)
library(dplyr)
library(readr)
library(glue)
library(purrr)
library(cytobox)
library(magrittr)
library(pheatmap)
library(monocle)

source(here(subdir, "analysis/functions.R"))
```


# Analysis

## Cell type density along normal interneuron differentiation trajectory

First we load a table of pseudotime values per cell:


```r
# Clean up the cell names, which got mangled a bit due to additional prefixes
# when the Seurat objects were merged
tidy_cellnames <- function(pseudotime) {
    
    pseudotime %>%
        mutate(cell = case_when(
            grepl("e12|e15", cell) ~ paste0("__", cell),
            grepl("p0", cell) ~ paste0("__", cell),
            TRUE ~ cell
        )) %>%
        mutate(pseudotime = Pseudotime)
    
}

inhibitory_pseudotime <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_inhibitory_neurons_pseudotime.txt")) %>%
    tidy_cellnames()
```

```
## Parsed with column specification:
## cols(
##   cell = col_character(),
##   Pseudotime = col_double()
## )
```

Next load the cell type identities for each cell:


```r
inhibitory_pseudotime$pseudotime %<>% scales::rescale(to = c(0, 100))

# Load the metadata for the joint forebrain and set the identities
forebrain_meta <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_joint_forebrain_metadata.tsv"))
```

```
## Parsed with column specification:
## cols(
##   Cell = col_character(),
##   nGene = col_double(),
##   nUMI = col_double(),
##   orig.ident = col_character(),
##   percent_mito = col_double(),
##   res.0.8 = col_double(),
##   ID_20190715_with_blacklist_and_refined = col_character(),
##   ID_20190715_joint_clustering = col_character(),
##   ID_20190730_with_blacklist_and_refined = col_character(),
##   Cluster = col_double()
## )
```

```r
inhibitory_pseudotime <- left_join(inhibitory_pseudotime, forebrain_meta, by = c("cell" = "Cell"))
```

Finally, generate a density plot of each cell type along pseudotime, to represent
the relative proportions of each cell type as differentiation proceeds:


```r
# Grab cells in the clusters included in this Monocle analysis,
# and plot their density in pseuodtime
inhibitory_pseudotime %>%
    filter(grepl("VRGC|INIP|MGINH|CINH", ID_20190730_with_blacklist_and_refined)) %>%
    separate(ID_20190730_with_blacklist_and_refined,
             sep = "_",
             into = c("Sample", "Cell_type")) %>%
    rr_ggplot(aes(x = Pseudotime), plot_num = 1) +
    geom_density(aes(fill = Cell_type), alpha = 0.8, size = 0.4) +
    scale_fill_manual(values = c(
        "VRGC"   = "#ffb219",
        "INIP-P" = "#8798C8",
        "MGINH"  = "#4a1cc9",
        "CINHN"  = "#135ca0"
    )) +
    theme_min()
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/01//interneuron_differentiation_density_plot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/01//interneuron_differentiation_density_plot...*]~</span>

## Plot genes of interest along pseudotime


```r
load(here("reference_datasets/2019_Jessa/Jessa2019_inhibitory_neurons.monocle.Rda"))

pal_pseudotime <- colorRampPalette(c("navy", "white", "red"))(n = 200)
```

One way to evaluate whether the genes which drive enrichment of fetal neuronal signatures
vs. depletion of pediatric/adult signatures in G34 gliomas by GSEA are related to maturity is by checking
their expression during normal differentiation. Get the genes from the GSEA analysis:


```r
signatures <- c("MGE newborn neurons 1", "Parvalbumin interneurons")

fgsea_df_idh_filt <- read_tsv(here("bulk_transcriptome_epigenome/output/02/fgsea_results_padj<0.01_IDH_top_le.tsv"))
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

```r
leading_edge <- fgsea_df_idh_filt %>% 
    filter(Cell_type %in% signatures) %>% 
    select(Cell_type, leadingEdge_top) %>% 
    tibble::deframe() %>% 
    map(~ .x %>% stringr::str_split(", ") %>% unlist() %>% hg2mm() %>% filter_NA())

leading_edge
```

```
## $`MGE newborn neurons 1`
##  [1] "Arx"           "Dlx6"          "Smarcd1"       "Plxna2"       
##  [5] "Dcx"           "St8sia5"       "E030030I06Rik" "Dlx5"         
##  [9] "Fam65b"        "Dlx2"          "Tmem2"         "Gad2"         
## [13] "Dlx1"          "Pdgfd"         "Metap1d"       "Arl4d"        
## [17] "Sox1"          "Ccser1"       
## 
## $`Parvalbumin interneurons`
##  [1] "Epha6"    "Fgf12"    "Nxph1"    "Gria4"    "Sparcl1"  "Cntnap5c"
##  [7] "Cntnap5b" "Cntnap5a" "Adcy8"    "Kcnab1"   "Frmd5"    "Plcxd3"  
## [13] "Caln1"    "Kcnd2"    "Dner"     "Vwc2"     "Fam19a2"  "Shisa9"  
## [19] "Ptchd4"
```

However, we only want to plot genes with reasonable mean exp in the clusters of
this lineage, so that the smoothing makes sense:


```r
all_genes <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_forebrain_gene_names.tsv")) %>%
    pull(genes)
```

```
## Parsed with column specification:
## cols(
##   genes = col_character()
## )
```

```r
filter_by_mean_expr <- function(genes) {
    
    genes <- genes[genes %in% all_genes] %>% unique()
    
    # Check expression in all the clusters which comprise the cortical
    # inhibitory neuron lineages
    feather::read_feather(here("reference_datasets/2019_Jessa/Jessa2019_mean_expression.feather"),
                          c("Cluster", genes)) %>%
        gather(Gene, Mean_expression, 2:ncol(.)) %>%
        filter(grepl("^F", Cluster) & grepl("e12|e15|p0", Cluster)) %>%
        filter(grepl("INIP|INH|VRGC", Cluster)) %>%
        filter(Mean_expression > 0.1) %>%
        pull(Gene) %>%
        unique() %>%
        sort()
    
}

(leading_edge_filt <- map(leading_edge, ~ filter_by_mean_expr(.x)))
```

```
## $`MGE newborn neurons 1`
##  [1] "Arl4d"   "Arx"     "Dcx"     "Dlx1"    "Dlx2"    "Dlx5"    "Dlx6"   
##  [8] "Fam65b"  "Gad2"    "Metap1d" "Plxna2"  "Smarcd1" "Sox1"    "Tmem2"  
## 
## $`Parvalbumin interneurons`
## [1] "Caln1"   "Dner"    "Fam19a2" "Frmd5"   "Gria4"   "Kcnd2"   "Nxph1"  
## [8] "Plcxd3"  "Sparcl1"
```

```r
names(leading_edge_filt) <- gsub(" ", "_", names(leading_edge_filt))
```

Extract genes to examine in pseudotime from the transcriptome/epigenome analyses -- we'll
gate genes based on high zK27Ac and low zK27me3, or vice versa:


```r
# Read in output of ChIP-seq analysis
multiomic_genes <- read_tsv(here("bulk_transcriptome_epigenome/output/03/DGE_and_histone_marks.tsv"))
```

```
## Parsed with column specification:
## cols(
##   ENSID = col_character(),
##   symbol = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double(),
##   G34R.zK27ac = col_double(),
##   G34R.zK27m3 = col_double(),
##   G34RV.K27ac.median = col_double(),
##   G34RV.K27me3.median = col_double(),
##   nonG34RV.K27ac.median = col_double(),
##   nonG34RV.K27me3.median = col_double()
## )
```

```r
multiomic_genes_top <- list("UP" = multiomic_genes %>%
                             filter(G34R.zK27ac > 0.6) %>%
                             pull(symbol) %>% hg2mm() %>% filter_NA(),
                         "DOWN" = multiomic_genes %>%
                             filter(G34R.zK27m3 > 0.6) %>%
                             pull(symbol) %>% hg2mm() %>% filter_NA())

(multiomic_genes_top_filt <- map(multiomic_genes_top, ~ filter_by_mean_expr(.x)))
```

```
## $UP
##  [1] "Cdca7l"  "Chic2"   "Dlx1"    "Dlx2"    "Dlx6"    "Efnb1"   "Elovl2" 
##  [8] "Ephb2"   "Foxg1"   "Foxm1"   "Gad2"    "Gsx2"    "Ildr2"   "Insm1"  
## [15] "Jag1"    "Npy"     "Nr2e1"   "Otx1"    "Pax6"    "Plk2"    "Prdm16" 
## [22] "Sel1l3"  "Skida1"  "Slc10a4" "Sox1"    "Sp8"     "Tacc3"   "Tead2"  
## 
## $DOWN
##  [1] "Arhgdig" "Cadps"   "Chd5"    "Dlgap3"  "Igsf21"  "Jph3"    "Jph4"   
##  [8] "Npas1"   "Ntrk2"   "Nxph1"   "Sertad4" "Sez6l2"  "Shank1"  "St8sia1"
```


Intersect RNA-seq leading genes and multi-omic genes:


```r
pseudotime_union <- list(
    # Newborn interneuron & upregulated genes
    "UP"   = unique(c(leading_edge_filt[[1]], multiomic_genes_top_filt[[1]])),
                         
    # Postnatal interneuron & downregulated genes
    "DOWN" = unique(c(leading_edge_filt[[2]], multiomic_genes_top_filt[[2]])))

imap(pseudotime_union, ~ ordered_pseudotime_heatmap(monocle[.x,],
                                                    show_rownames = TRUE,
                                                    order_by_pseudotime = TRUE,
                                                    hmcols = pal_pseudotime,
                                                    prefix = glue("{figout}/pseudotime_heatmap{.y}")))
```

```
## $UP
## 
## $DOWN
```

```r
knitr::include_graphics(glue("{figout}/pseudotime_heatmapUP.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/01///pseudotime_heatmapUP.png" width="2100" />

```r
knitr::include_graphics(glue("{figout}/pseudotime_heatmapDOWN.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/01///pseudotime_heatmapDOWN.png" width="2100" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/01//pseudotime_heatmap...*]~</span>


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-09-16 10:51:08
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [feed0d9] 2020-09-16: Update lockfile (after installing monocle)
```

The random seed was set with `set.seed(100)`

The R session info:
<details>

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/2018.3.222/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so
## 
## locale:
##  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
##  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
##  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] splines   stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] monocle_2.10.1      DDRTree_0.1.5       irlba_2.3.3        
##  [4] VGAM_1.1-3          ggplot2_3.1.0       Biobase_2.42.0     
##  [7] BiocGenerics_0.28.0 Matrix_1.2-14       pheatmap_1.0.12    
## [10] magrittr_1.5        cytobox_0.6.1       purrr_0.3.4        
## [13] glue_1.4.2          readr_1.3.1         dplyr_0.8.0        
## [16] tidyr_0.8.2         here_0.1           
## 
## loaded via a namespace (and not attached):
##   [1] snow_0.4-3           backports_1.1.9      Hmisc_4.2-0         
##   [4] sn_1.6-2             plyr_1.8.6           igraph_1.2.5        
##   [7] lazyeval_0.2.2       densityClust_0.3     fastICA_1.2-2       
##  [10] TH.data_1.0-10       digest_0.6.25        foreach_1.5.0       
##  [13] htmltools_0.5.0      viridis_0.5.1        lars_1.2            
##  [16] gdata_2.18.0         checkmate_2.0.0      cluster_2.0.7-1     
##  [19] mixtools_1.2.0       ROCR_1.0-7           limma_3.38.3        
##  [22] matrixStats_0.56.0   docopt_0.7.1         R.utils_2.10.1      
##  [25] sandwich_2.5-1       colorspace_1.4-1     ggrepel_0.8.0       
##  [28] xfun_0.17            sparsesvd_0.2        jsonlite_1.7.1      
##  [31] crayon_1.3.4         survival_2.41-3      zoo_1.8-8           
##  [34] iterators_1.0.12     ape_5.4-1            gtable_0.3.0        
##  [37] kernlab_0.9-29       prabclus_2.3-2       DEoptimR_1.0-8      
##  [40] scales_1.1.1         mvtnorm_1.1-1        bibtex_0.4.2.2      
##  [43] Rcpp_1.0.5           metap_1.4            dtw_1.21-3          
##  [46] plotrix_3.7-8        viridisLite_0.3.0    htmlTable_2.0.1     
##  [49] tmvnsim_1.0-2        reticulate_1.16      foreign_0.8-70      
##  [52] bit_4.0.4            proxy_0.4-24         mclust_5.4.6        
##  [55] SDMTools_1.1-221.2   Formula_1.2-3        tsne_0.1-3          
##  [58] htmlwidgets_1.5.1    httr_1.4.2           FNN_1.1.3           
##  [61] gplots_3.0.1.1       RColorBrewer_1.1-2   fpc_2.2-7           
##  [64] acepack_1.4.1        TFisher_0.2.0        modeltools_0.2-23   
##  [67] ellipsis_0.3.1       Seurat_2.3.4         ica_1.0-2           
##  [70] pkgconfig_2.0.3      R.methodsS3_1.8.1    flexmix_2.3-15      
##  [73] nnet_7.3-12          tidyselect_1.1.0     rlang_0.4.7         
##  [76] reshape2_1.4.4       munsell_0.5.0        tools_3.5.1         
##  [79] mathjaxr_1.0-1       ggridges_0.5.2       evaluate_0.14       
##  [82] stringr_1.4.0        yaml_2.2.1           knitr_1.29          
##  [85] bit64_4.0.5          fitdistrplus_1.1-1   robustbase_0.93-6   
##  [88] caTools_1.17.1.1     randomForest_4.6-14  RANN_2.6.1          
##  [91] pbapply_1.4-3        nlme_3.1-137         slam_0.1-47         
##  [94] R.oo_1.24.0          hdf5r_1.3.3          compiler_3.5.1      
##  [97] rstudioapi_0.11      png_0.1-7            tibble_3.0.3        
## [100] stringi_1.5.3        lattice_0.20-35      HSMMSingleCell_1.2.0
## [103] multtest_2.38.0      vctrs_0.3.4          mutoss_0.1-12       
## [106] pillar_1.4.6         lifecycle_0.2.0      combinat_0.0-8      
## [109] Rdpack_1.0.0         lmtest_0.9-38        data.table_1.13.0   
## [112] cowplot_0.9.4        bitops_1.0-6         gbRd_0.4-11         
## [115] R6_2.4.1             latticeExtra_0.6-28  KernSmooth_2.23-15  
## [118] gridExtra_2.3        codetools_0.2-15     MASS_7.3-50         
## [121] gtools_3.8.2         assertthat_0.2.1     rprojroot_1.3-2     
## [124] withr_2.2.0          qlcMatrix_0.9.7      mnormt_2.0.2        
## [127] multcomp_1.4-13      diptest_0.75-7       doSNOW_1.0.18       
## [130] hms_0.5.3            grid_3.5.1           rpart_4.1-13        
## [133] class_7.3-14         rmarkdown_1.11       segmented_1.2-0     
## [136] Rtsne_0.15           git2r_0.27.1         numDeriv_2016.8-1.1 
## [139] base64enc_0.1-3
```

</details>


***

<!-- END OF END MATTER -->
