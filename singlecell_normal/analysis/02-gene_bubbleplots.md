---
title: "02 - Expression of individual genes during development"
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
## Document index: 02
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
## G34-gliomas-repo/singlecell_normal/output/02
```

```
## G34-gliomas-repo/singlecell_normal/figures/02//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Demonstrate the expression patterns of select genes during normal development / 
under normal conditions, as a reference for expression patterns in the G34R/V tumors.

# Libraries


```r
# Load liraries here
library(here)

library(tidyverse)
library(Seurat)
library(pvclust)
library(dendextend)
library(glue)

source("../../../../misc/sjlib.R")
source("functions.R")
```

# Analysis

## Mouse developing forebrain

Get expression & signatures:


```r
atlas_signatures <- readRDS(here("bulk/output/02/atlas_signatures.Rds"))

# Get genes for neuro-ectoderm mouse forebrain clusters
atlas_signatures <- atlas_signatures$mm_sym[grepl("^F", names(atlas_signatures$mm_sym)) & !grepl("MGL|MAC|ENDO|PERI|MNG", names(atlas_signatures$mm_sym))]

blacklist <- read_tsv(here("bulk/data/Jessa2019_Table_2a.tsv")) %>%
  select(Cluster, Signature) %>%
  filter(is.na(Signature)) %>%
  pull(Cluster) %>%
  {gsub("_", " ", .)}
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
# Removed one mixed progenitor cluster which is not precisely labeled
blacklist <- c(blacklist, "F-p6 MIXP-P")
```

Construct dendrogram by computing the correlations:


```r
load(here("singlecell_normal/data/Jessa2019_cluster_mean_expr_matrix.Rda"))

# Only keep forebrain neuroectoderm clusters
meanexp <- meanexp[, grepl("^F", colnames(meanexp)) & !grepl("MGL|MAC|ENDO|PERI|MNG", colnames(meanexp))]

# Get the unique sig genes from mouse only
atlas_unique_genes_mm <- atlas_signatures %>% unlist %>% unique

# Subset to these genes, and filter out blacklist clusters
meanexp_mm_uniq <- meanexp[atlas_unique_genes_mm,
                           !(colnames(meanexp) %in% blacklist)]

# Take correlations
cor_mm <- cytokit::correlateExpression(meanexp_mm_uniq, meanexp_mm_uniq,
                                       genes = atlas_unique_genes_mm,
                                       from_sp = "mm",
                                       to_sp = "mm")
```

```
## Computing correlations...
```

```r
# Get matrix without NAs
meanexp_mm_uniq_no_NA <- meanexp_mm_uniq
meanexp_mm_uniq_no_NA[is.na(meanexp_mm_uniq_no_NA)] <- 0
```

Next we bootstrap hierarchical clustering over the correlations,
using the `pvclust` package.


```r
result <- pvclust(meanexp_mm_uniq_no_NA, method.dist = spearman,
                  method.hclust = "complete", nboot = 100)
```

```
## Bootstrap (r = 0.5)... Done.
## Bootstrap (r = 0.6)... Done.
## Bootstrap (r = 0.7)... Done.
## Bootstrap (r = 0.8)... Done.
## Bootstrap (r = 0.9)... Done.
## Bootstrap (r = 1.0)... Done.
## Bootstrap (r = 1.1)... Done.
## Bootstrap (r = 1.2)... Done.
## Bootstrap (r = 1.3)... Done.
## Bootstrap (r = 1.4)... Done.
```

```r
dend_mm <- as.dendrogram(result)

dendrogram_order_atlas <- colnames(cor_mm)[order.dendrogram(dend_mm)]

rr_saveRDS(dendrogram_order_atlas,
           file = glue("{out}/dendrogram_order_atlas.Rda"),
           desc = "Order of clusters in dendrogram constructed based on gene expression correlation.")
```

```
## ...writing description of dendrogram_order_atlas.Rda to G34-gliomas-repo/singlecell_normal/output/02/dendrogram_order_atlas.Rda
```



```r
par(mar = c(20,2,2,2))
plot(dend_mm)
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/singlecell_normal/figures/02//plot_pvclust_complete_atlas-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/singlecell_normal/figures/02//plot_pvclust_complete_atlas...*]~</span>

Generate the plot by loading the mean expression per cluster for select
genes, and the proportion of cells in each cluster expressing the gene.


```r
genes <- c("Gapdh", "Pdgfra", "Gsx2", "Dlx1", "Dlx2")

goi <- feather::read_feather(here("singlecell_normal/data/Jessa2019_mean_expression.feather"),
                             c("Cluster", genes)) %>%
  tibble::column_to_rownames(var = "Cluster") %>%
  apply(2, scales::rescale) %>%
  as.data.frame %>%
  tibble::add_column(Cluster = rownames(.), .before = 1) %>%
  gather(Gene, Expression, 2:ncol(.)) %>%
  filter(Cluster %in% colnames(meanexp_mm_uniq_no_NA))

pct1 <- feather::read_feather(here("singlecell_normal/data/Jessa2019_pct1.feather"),
                              c("Cluster", genes)) %>%
  gather(Gene, Pct1, 2:ncol(.)) %>%
  filter(!is.na(Pct1)) %>%
  filter(Cluster %in% colnames(meanexp_mm_uniq_no_NA))

bubbleplot_data <- left_join(goi, pct1, by = c("Cluster", "Gene")) %>%
  filter(Pct1 > 0) %>%
  mutate(Cluster = factor(Cluster, levels = dendrogram_order_atlas)) %>%
  mutate(Gene = factor(Gene, levels = rev(genes))) %>%
  filter(!is.na(Cluster))

bubbleplot_data %>%
  rr_ggplot(aes(x = Cluster, y = Gene), plot_num = 1) +
  geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
  scale_radius() +
  scale_color_gradientn(colours = tail(rdbu, 70)) +
  theme_min() +
  rotateX() +
  theme(panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 13))
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/singlecell_normal/figures/02//atlas_bubbleplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/singlecell_normal/figures/02//atlas_bubbleplot...*]~</span>


<!-- ## Human -->

<!-- ### Get expression / signatures -->

<!-- ```{r get_nowakowski} -->

<!-- load("../../../../ext_datasets/2017-nowakowski_human_fetal_cortex/20191016-processing_SJ/output/01/nowakowski_seurat.Rda") -->
<!-- load("output/02/atlas_signatures_formatted.Rda") -->

<!-- nowakowski_meanexp <- cytobox::meanClusterExpression(nowa, genes = c("GAPDH", unique(unlist(nowakowski_signatures$hg_sym)))) -->

<!-- nowakowski_meanexp <- nowakowski_meanexp[, !(colnames(nowakowski_meanexp) %in% c("Glyc", "U1", "U2", "U3", "U4", "None"))] -->

<!-- cor_nowa <- cytokit::correlateExpression(nowakowski_meanexp, nowakowski_meanexp, -->
<!--                                          genes = unique(unlist(nowakowski_signatures$hg_sym)), -->
<!--                                          from_sp = "nowa", -->
<!--                                          to_sp = "nowa") -->

<!-- # Get matrix without NAs -->
<!-- meanexp_nowa_uniq_no_NA <- nowakowski_meanexp -->
<!-- meanexp_nowa_uniq_no_NA[is.na(meanexp_nowa_uniq_no_NA)] <- 0 -->

<!-- # Save inputs -->
<!-- save(meanexp_nowa_uniq_no_NA, cor_nowa, file = glue("{out}/pvclust_inputs_nowakowski.Rda")) -->

<!-- ``` -->

<!-- ```{r pvclust_complete2, dependson = 'get_nowakowski'} -->

<!-- result_nowa <- pvclust(meanexp_nowa_uniq_no_NA, method.dist = spearman, -->
<!--                        method.hclust = "complete", nboot = 100) -->

<!-- dend_nowa <- as.dendrogram(result_nowa) -->
<!-- labels_colors(dend_nowa) <- nowa@misc$colours[colnames(meanexp_nowa_uniq_no_NA)][order.dendrogram(dend_nowa)] -->

<!-- dendrogram_order_nowa <- colnames(cor_nowa)[order.dendrogram(dend_nowa)] -->
<!-- save(dendrogram_order_nowa, file = glue("{out}/complete_dendrogram_order_nowa.Rda")) -->

<!-- ``` -->


<!-- ```{r plot_pvclust_complete2, fig.height = 7, fig.width = 10, dependson = 'get_nowakowski'} -->

<!-- par(mar = c(20,2,2,2)) -->
<!-- plot(dend_nowa) -->

<!-- ``` -->

<!-- ### Bubbleplot -->

<!-- ```{r nowa_bubbleplot_prep, dependson = 'get_nowakowski'} -->

<!-- prop <- function(x) { -->

<!--   sum(x)/length(x) -->

<!-- } -->

<!-- calc_pct1 <- function(seurat, -->
<!--                       cluster_col, -->
<!--                       goi) { -->

<!--   require(data.table) -->

<!--   message("## Processing ", seurat@project.name) -->

<!--   seurat <- SetAllIdent(seurat, cluster_col) -->

<!--   x <- as.matrix(seurat@data)[goi, ] -->
<!--   # Binarize -->
<!--   x[x > 0] <- 1 -->
<!--   x <- as.data.table(t(x)) -->
<!--   # Add cluster info for cells -->
<!--   x[, Cluster := as.character(seurat@ident)] -->
<!--   # Get prop of cells expressing a gene, within each cluster -->
<!--   x[, lapply(.SD, prop), by = Cluster] -->

<!-- } -->

<!-- nowa_pct1 <- calc_pct1(nowa, cluster_col = "WGCNAcluster", -->
<!--                        goi = c("GAPDH", "MOXD1", "PDGFRA", -->
<!--                                "GSX2", "DLX1", "DLX2", -->
<!--                                "EOMES", "NEUROD2")) %>%  -->
<!--   data.frame() %>%  -->
<!--   gather(Gene, Pct1, 2:ncol(.)) -->

<!-- nowa_meanexp_goi <- nowakowski_meanexp[c("GAPDH", "MOXD1", "PDGFRA", "DLX1", "DLX2", "EOMES", "NEUROD2"), ] %>%  -->
<!--   t() %>%  -->
<!--   apply(2, scales::rescale) %>%  -->
<!--   data.frame() %>%  -->
<!--   tibble::rownames_to_column(var = "Cluster") %>%  -->
<!--   gather(Gene, Expression, 2:ncol(.)) -->

<!-- ``` -->

<!-- ```{r nowakowski_bubbleplot, dependson = c('pvclust_complete2', 'nowa_bubbleplot_prep'), fig.width = 10, fig.height = 3.5} -->

<!-- genes <- c("GAPDH", "MOXD1", "PDGFRA", "DLX1", "DLX2", "EOMES", "NEUROD2") -->

<!-- bubbleplot_data <- left_join(nowa_meanexp_goi, nowa_pct1, by = c("Cluster", "Gene")) %>% -->
<!--   filter(Pct1 > 0) %>% -->
<!--   mutate(Cluster = factor(Cluster, levels = dendrogram_order_nowa)) %>% -->
<!--   mutate(Gene = factor(Gene, levels = rev(genes))) %>% -->
<!--   filter(!is.na(Cluster)) %T>% -->
<!--   write_tsv(glue("{out}/bubbleplot_data_nowakowski.tsv")) -->

<!-- bubbleplot_data %>% -->
<!--   ggplot(aes(x = Cluster, y = Gene)) + -->
<!--   geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) + -->
<!--   scale_radius() + -->
<!--   scale_color_gradientn(colours = tail(rdbu, 70)) + -->
<!--   theme_min() + -->
<!--   rotateX() + -->
<!--   theme(panel.grid.major.x = element_line(colour = "grey90"), -->
<!--         panel.border = element_blank(), -->
<!--         axis.ticks.x = element_blank(), -->
<!--         axis.ticks.y = element_blank(), -->
<!--         axis.text.y = element_text(size = 13)) -->

<!-- ``` -->




## Striatal SVZ

Repeat the same analysis for the striatal SVZ dataset from Anderson et al:


```r
seurat_anderson     <- readRDS(here("singlecell_normal/data/Anderson2020_seurat.Rds"))
anderson_signatures <- readRDS(here("bulk/data/Anderson2020_signatures.Rds"))

meanexp_anderson    <- cytobox::meanClusterExpression(seurat_anderson,
                                                      genes = c(genes, unique(unlist(anderson_signatures$mm_sym))))
```

```
## Computing cluster means for Anderson_Foxp1_striatum...
```

```r
cor_anderson <- cytokit::correlateExpression(meanexp_anderson, meanexp_anderson,
                                             genes = unique(unlist(anderson_signatures$mm_sym)),
                                             from_sp = "mm",
                                             to_sp = "mm")
```

```
## Computing correlations...
```

```r
# Get matrix without NAs
meanexp_anderson_uniq_no_NA <- meanexp_anderson
meanexp_anderson_uniq_no_NA[is.na(meanexp_anderson_uniq_no_NA)] <- 0
```

Create the dendrogram:


```r
result_anderson <- pvclust(meanexp_anderson_uniq_no_NA, method.dist = spearman,
                           method.hclust = "complete", nboot = 100)
```

```
## Bootstrap (r = 0.5)... Done.
## Bootstrap (r = 0.6)... Done.
## Bootstrap (r = 0.7)... Done.
## Bootstrap (r = 0.8)... Done.
## Bootstrap (r = 0.9)... Done.
## Bootstrap (r = 1.0)... Done.
## Bootstrap (r = 1.1)... Done.
## Bootstrap (r = 1.2)... Done.
## Bootstrap (r = 1.3)... Done.
## Bootstrap (r = 1.4)... Done.
```

```r
dend_anderson <- as.dendrogram(result_anderson)

dendrogram_order_anderson <- colnames(cor_anderson)[order.dendrogram(dend_anderson)]

rr_saveRDS(dendrogram_order_anderson,
           file = glue("{out}/dendrogram_order_anderson.Rda"),
           desc = "Order of clusters in dendrogram constructed based on gene expression correlation.")
```

```
## ...writing description of dendrogram_order_anderson.Rda to G34-gliomas-repo/singlecell_normal/output/02/dendrogram_order_anderson.Rda
```

Show the dendrogram:


```r
par(mar = c(20,2,2,2))
plot(dend_anderson)
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/singlecell_normal/figures/02//plot_pvclust_complete_anderson-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/singlecell_normal/figures/02//plot_pvclust_complete_anderson...*]~</span>

Construct the data for the bubbleplot, including the proportion of cells in each cluster
where each gene is detected, and the mean expression of our genes of interest:


```r
pct1_anderson <- calc_pct1(seurat_anderson, cluster_col = "Cluster",
                           genes = genes) %>%
  data.frame() %>%
  gather(Gene, Pct1, 2:ncol(.))
```

```
## Loading required package: data.table
```

```
## 
## Attaching package: 'data.table'
```

```
## The following object is masked from 'package:dendextend':
## 
##     set
```

```
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
```

```
## The following object is masked from 'package:purrr':
## 
##     transpose
```

```
## ## Processing Anderson_Foxp1_striatum
```

```r
meanexp_anderson_goi <- meanexp_anderson[genes, ] %>%
  t() %>%
  apply(2, scales::rescale) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "Cluster") %>%
  gather(Gene, Expression, 2:ncol(.))
```



```r
bubbleplot_data <- left_join(meanexp_anderson_goi, pct1_anderson, by = c("Cluster", "Gene")) %>%
  filter(Pct1 > 0) %>%
  mutate(Cluster = factor(Cluster, levels = dendrogram_order_anderson)) %>%
  mutate(Gene = factor(Gene, levels = rev(genes))) %>%
  filter(!is.na(Cluster))

bubbleplot_data %>%
  rr_ggplot(aes(x = Cluster, y = Gene), plot_num = 1) +
  geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
  scale_radius() +
  scale_color_gradientn(colours = tail(rdbu, 70)) +
  theme_min() +
  rotateX() +
  theme(panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 13))
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/singlecell_normal/figures/02//anderson_bubbleplot-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-07-20 15:09:00
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
##  [1] data.table_1.12.0 bindrcpp_0.2.2    glue_1.4.1       
##  [4] dendextend_1.9.0  pvclust_2.0-0     Seurat_2.3.4     
##  [7] Matrix_1.2-14     cowplot_0.9.4     forcats_0.3.0    
## [10] stringr_1.3.1     dplyr_0.7.7       purrr_0.3.4      
## [13] readr_1.3.1       tidyr_0.8.2       tibble_3.0.1     
## [16] ggplot2_3.1.0     tidyverse_1.2.1   here_0.1         
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.2.0        snow_0.4-3          backports_1.1.3    
##   [4] Hmisc_4.2-0         plyr_1.8.4          igraph_1.2.2       
##   [7] lazyeval_0.2.1      splines_3.5.0       digest_0.6.16      
##  [10] foreach_1.4.4       htmltools_0.3.6     viridis_0.5.1      
##  [13] lars_1.2            gdata_2.18.0        magrittr_1.5       
##  [16] checkmate_1.9.1     cluster_2.0.7-1     mixtools_1.1.0     
##  [19] ROCR_1.0-7          modelr_0.1.3        R.utils_2.7.0      
##  [22] colorspace_1.4-0    rvest_0.3.2         haven_2.0.0        
##  [25] xfun_0.12           crayon_1.3.4        jsonlite_1.6       
##  [28] bindr_0.1.1         survival_2.41-3     zoo_1.8-4          
##  [31] iterators_1.0.10    ape_5.2             gtable_0.2.0       
##  [34] kernlab_0.9-27      prabclus_2.2-7      DEoptimR_1.0-8     
##  [37] scales_1.0.0        mvtnorm_1.0-10      bibtex_0.4.2       
##  [40] Rcpp_1.0.0          metap_1.1           dtw_1.20-1         
##  [43] viridisLite_0.3.0   htmlTable_1.13.1    reticulate_1.10    
##  [46] foreign_0.8-70      bit_1.1-14          proxy_0.4-22       
##  [49] mclust_5.4.2        SDMTools_1.1-221    Formula_1.2-3      
##  [52] tsne_0.1-3          stats4_3.5.0        htmlwidgets_1.3    
##  [55] httr_1.4.0          gplots_3.0.1.1      RColorBrewer_1.1-2 
##  [58] fpc_2.1-11.1        acepack_1.4.1       modeltools_0.2-22  
##  [61] ellipsis_0.2.0.1    ica_1.0-2           pkgconfig_2.0.2    
##  [64] R.methodsS3_1.7.1   flexmix_2.3-14      nnet_7.3-12        
##  [67] labeling_0.3        reshape2_1.4.3      tidyselect_1.1.0   
##  [70] rlang_0.4.6         munsell_0.5.0       cellranger_1.1.0   
##  [73] tools_3.5.0         cli_1.0.1           generics_0.0.2     
##  [76] broom_0.5.1         ggridges_0.5.1      evaluate_0.12      
##  [79] yaml_2.2.0          npsurv_0.4-0        knitr_1.21         
##  [82] bit64_0.9-7         fitdistrplus_1.0-14 robustbase_0.93-2  
##  [85] caTools_1.17.1.1    RANN_2.6            pbapply_1.4-0      
##  [88] nlme_3.1-137        whisker_0.3-2       R.oo_1.22.0        
##  [91] xml2_1.2.0          hdf5r_1.0.0         compiler_3.5.0     
##  [94] rstudioapi_0.9.0    png_0.1-7           lsei_1.2-0         
##  [97] stringi_1.2.4       lattice_0.20-35     trimcluster_0.1-2.1
## [100] vctrs_0.3.1         pillar_1.4.4        lifecycle_0.2.0    
## [103] BiocManager_1.30.10 Rdpack_0.10-1       lmtest_0.9-36      
## [106] bitops_1.0-6        irlba_2.3.3         gbRd_0.4-11        
## [109] R6_2.3.0            latticeExtra_0.6-28 renv_0.10.0        
## [112] KernSmooth_2.23-15  gridExtra_2.3       codetools_0.2-15   
## [115] MASS_7.3-49         gtools_3.8.1        assertthat_0.2.0   
## [118] rprojroot_1.3-2     withr_2.1.2         diptest_0.75-7     
## [121] parallel_3.5.0      doSNOW_1.0.16       hms_0.4.2          
## [124] grid_3.5.0          rpart_4.1-13        class_7.3-14       
## [127] rmarkdown_1.11      segmented_0.5-3.0   Rtsne_0.15         
## [130] git2r_0.27.1        lubridate_1.7.4     base64enc_0.1-3
```

</details>


***

<!-- END OF END MATTER -->
