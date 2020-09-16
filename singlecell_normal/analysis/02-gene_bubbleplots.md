---
title: "02 - Expression of individual genes during development"
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
## Document index: 02
```

```r
# Specify where to save outputs
out        <- here(subdir, "output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here(subdir, "figures", doc_id, "/"); dir.create(figout, recursive = TRUE)
cache      <- paste0("~/tmp/", basename(here()), "/", subdir, "/", doc_id, "/")

message("Cache: ", cache)
```

```
## Cache: ~/tmp/G34-gliomas/singlecell_normal/02/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/singlecell_normal/output/02
```

```
## G34-gliomas/singlecell_normal/figures/02//
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
# Load libraries here
library(tidyr)
library(dplyr)
library(readr)
library(Seurat)
library(pvclust)
library(dendextend)
library(glue)
library(feather)

source(here(subdir, "analysis/functions.R"))
```

# Analysis

## Mouse developing forebrain

Get expression & signatures:


```r
atlas_signatures <- readRDS(here("bulk_transcriptome_epigenome/output/02/atlas_signatures.Rds"))

# Get genes for neuro-ectoderm mouse forebrain clusters
atlas_signatures <- atlas_signatures$mm_sym[grepl("^F", names(atlas_signatures$mm_sym)) & !grepl("MGL|MAC|ENDO|PERI|MNG", names(atlas_signatures$mm_sym))]

blacklist <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_Table_2a.tsv")) %>%
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
load(here("reference_datasets/2019_Jessa/Jessa2019_cluster_mean_expr_matrix.Rda"))

# Only keep forebrain neuroectoderm clusters
meanexp <- meanexp[, grepl("^F", colnames(meanexp)) & !grepl("MGL|MAC|ENDO|PERI|MNG", colnames(meanexp))]

# Get the unique sig genes from mouse only
atlas_unique_genes_mm <- atlas_signatures %>% unlist %>% unique

# Subset to these genes, and filter out blacklist clusters
meanexp_mm_uniq <- meanexp[atlas_unique_genes_mm,
                           !(colnames(meanexp) %in% blacklist)]

# Take correlations
cor_mm <- cytobox::correlateExpression(meanexp_mm_uniq, meanexp_mm_uniq,
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
# Run pvclust
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
# Convert to a dendrogram
dend_mm <- as.dendrogram(result)

dendrogram_order_atlas <- colnames(cor_mm)[order.dendrogram(dend_mm)]

# Save the leaf order to use for plotting
rr_saveRDS(dendrogram_order_atlas,
           file = glue("{out}/dendrogram_order_atlas.Rda"),
           desc = "Order of clusters in dendrogram constructed based on gene expression correlation.")
```

```
## ...writing description of dendrogram_order_atlas.Rda to G34-gliomas/singlecell_normal/output/02/dendrogram_order_atlas.Rda
```



```r
par(mar = c(20,2,2,2))
plot(dend_mm)
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//plot_pvclust_complete_atlas-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//plot_pvclust_complete_atlas...*]~</span>

Generate the plot by loading the mean expression per cluster for select
genes, and the proportion of cells in each cluster expressing the gene.


```r
# The reason we include Gapdh here is because it's expressed in all clusters,
# and hence will force plotting all clusters, even if they don't express the other
# genes of interest which are actually cell-type specific
genes <- c("Gapdh", "Pdgfra", "Gsx2", "Dlx1", "Dlx2")

# Mean expression per cluster
goi <- feather::read_feather(here("reference_datasets/2019_Jessa/Jessa2019_mean_expression.feather"),
                             c("Cluster", genes)) %>%
  tibble::column_to_rownames(var = "Cluster") %>%
  apply(2, scales::rescale) %>%
  as.data.frame %>%
  tibble::add_column(Cluster = rownames(.), .before = 1) %>%
  gather(Gene, Expression, 2:ncol(.)) %>%
  filter(Cluster %in% colnames(meanexp_mm_uniq_no_NA))

# pct1 = proportion of cells in each cluster in which each gene is detected
pct1 <- feather::read_feather(here("reference_datasets/2019_Jessa/Jessa2019_pct1.feather"),
                              c("Cluster", genes)) %>%
  gather(Gene, Pct1, 2:ncol(.)) %>%
  filter(!is.na(Pct1)) %>%
  filter(Cluster %in% colnames(meanexp_mm_uniq_no_NA))

# Combine to generate a table with both values
bubbleplot_data <- left_join(goi, pct1, by = c("Cluster", "Gene")) %>%
  filter(Pct1 > 0) %>%
  mutate(Cluster = factor(Cluster, levels = dendrogram_order_atlas)) %>%
  mutate(Gene = factor(Gene, levels = rev(genes))) %>%
  filter(!is.na(Cluster))

# Generate bubbleplot, using dendrogram order as computed
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

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//atlas_bubbleplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//atlas_bubbleplot...*]~</span>


## Human

Construct the mean expression matrix and compute correlations for the human
fetal brain dataset:


```r
load(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski_seurat.Rda"))
load(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski_signatures.Rda"))

nowakowski_meanexp <- cytobox::meanClusterExpression(nowa, genes = c("GAPDH", unique(unlist(nowakowski_signatures$hg_sym))))
```

```
## Computing cluster means for Nowakowski human fetal telencephalon...
```

```r
# Remove unknown clusters
nowakowski_meanexp <- nowakowski_meanexp[, !(colnames(nowakowski_meanexp) %in% c("Glyc", "U1", "U2", "U3", "U4", "None"))]

cor_nowa <- cytobox::correlateExpression(nowakowski_meanexp, nowakowski_meanexp,
                                         genes = unique(unlist(nowakowski_signatures$hg_sym)),
                                         from_sp = "nowa",
                                         to_sp = "nowa")
```

```
## Computing correlations...
```

```r
# Get matrix without NAs
meanexp_nowa_uniq_no_NA <- nowakowski_meanexp
meanexp_nowa_uniq_no_NA[is.na(meanexp_nowa_uniq_no_NA)] <- 0
```


```r
result_nowa <- pvclust(meanexp_nowa_uniq_no_NA, method.dist = spearman,
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
dend_nowa <- as.dendrogram(result_nowa)
dendrogram_order_nowa <- colnames(cor_nowa)[order.dendrogram(dend_nowa)]

rr_saveRDS(dendrogram_order_nowa,
           file = glue("{out}/dendrogram_order_nowakowski.Rda"),
           desc = "Order of clusters in dendrogram constructed based on gene expression correlation.")
```

```
## ...writing description of dendrogram_order_nowakowski.Rda to G34-gliomas/singlecell_normal/output/02/dendrogram_order_nowakowski.Rda
```



```r
par(mar = c(20,2,2,2))
plot(dend_nowa)
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//plot_pvclust_complete_nowakowski-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//plot_pvclust_complete_nowakowski...*]~</span>

Calculate the inputs for the bubbleplot:


```r
nowa_pct1 <- calc_pct1(nowa, cluster_col = "WGCNAcluster",
                       genes = c("GAPDH", "PDGFRA", "GSX2", "DLX1", "DLX2")) %>%
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
## ## Processing Nowakowski human fetal telencephalon
```

```r
# Scaling expression to [0,1]
nowa_meanexp_goi <- nowakowski_meanexp[c("GAPDH", "PDGFRA", "DLX1", "DLX2"), ] %>%
  t() %>%
  apply(2, scales::rescale) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "Cluster") %>%
  gather(Gene, Expression, 2:ncol(.))
```

Create bubbleplot:


```r
genes <- c("GAPDH", "MOXD1", "PDGFRA", "DLX1", "DLX2")

bubbleplot_data <- left_join(nowa_meanexp_goi, nowa_pct1, by = c("Cluster", "Gene")) %>%
  filter(Pct1 > 0) %>%
  mutate(Cluster = factor(Cluster, levels = dendrogram_order_nowa)) %>%
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

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//nowakowski_bubbleplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//nowakowski_bubbleplot...*]~</span>


## Striatal SVZ

Repeat the same analysis for the striatal SVZ dataset from Anderson et al:


```r
genes <- c("Gapdh", "Pdgfra", "Gsx2", "Dlx1", "Dlx2")

seurat_anderson     <- readRDS(here("reference_datasets/2020_Anderson/processed_data/01-seurat.Rds"))
anderson_signatures <- readRDS(here("reference_datasets/2020_Anderson/processed_data/02-cluster_gene_signatures.Rds"))

meanexp_anderson    <- cytobox::meanClusterExpression(seurat_anderson,
                                                      genes = c(genes, unique(unlist(anderson_signatures$mm_sym))))
```

```
## Computing cluster means for Anderson_Foxp1_striatum...
```

```r
cor_anderson <- cytobox::correlateExpression(meanexp_anderson, meanexp_anderson,
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
## ...writing description of dendrogram_order_anderson.Rda to G34-gliomas/singlecell_normal/output/02/dendrogram_order_anderson.Rda
```

Show the dendrogram:


```r
par(mar = c(20,2,2,2))
plot(dend_anderson)
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//plot_pvclust_complete_anderson-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//plot_pvclust_complete_anderson...*]~</span>

Construct the data for the bubbleplot, including the proportion of cells in each cluster
where each gene is detected, and the mean expression of our genes of interest:


```r
pct1_anderson <- calc_pct1(seurat_anderson, cluster_col = "Cluster",
                           genes = genes) %>%
  data.frame() %>%
  gather(Gene, Pct1, 2:ncol(.))
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

Finally generate the bubbleplot in the same order as the dendrogram constructed above:


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

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//anderson_bubbleplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//anderson_bubbleplot...*]~</span>

Separately, generate a bubbleplot focused on clusters of interest:

Finally generate the bubbleplot in the same order as the dendrogram constructed above:


```r
select_cluster_order <- c("15-OPCs", "31-OPCs", "4-OPCs", "16-OPCs",
                          "33-Neurogenic progenitor",
                          "7-Neurogenic progenitor",
                          "5-Neurogenic progenitor", 
                          "8-Neurogenic progenitor",
                          "19-Neurogenic progenitor")

bubbleplot_data <- left_join(meanexp_anderson_goi, pct1_anderson, by = c("Cluster", "Gene")) %>%
  filter(Pct1 > 0) %>%
  mutate(Cluster = factor(Cluster, levels = select_cluster_order)) %>%
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

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/02//anderson_bubbleplot_select-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/02//anderson_bubbleplot_select...*]~</span>


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-09-16 07:39:53
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [73f10e6] 2020-09-16: Regenerate HTMLs for bulk analysis
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
## [1] stats     graphics  grDevices datasets  utils     methods   base     
## 
## other attached packages:
##  [1] data.table_1.13.0 feather_0.3.5     glue_1.4.2        dendextend_1.14.0
##  [5] pvclust_2.2-0     Seurat_2.3.4      Matrix_1.2-14     cowplot_0.9.4    
##  [9] ggplot2_3.1.0     readr_1.3.1       dplyr_0.8.0       tidyr_0.8.2      
## [13] here_0.1         
## 
## loaded via a namespace (and not attached):
##   [1] snow_0.4-3          backports_1.1.9     Hmisc_4.2-0        
##   [4] sn_1.6-2            plyr_1.8.6          igraph_1.2.5       
##   [7] lazyeval_0.2.2      splines_3.5.1       TH.data_1.0-10     
##  [10] digest_0.6.25       foreach_1.5.0       htmltools_0.5.0    
##  [13] viridis_0.5.1       lars_1.2            gdata_2.18.0       
##  [16] magrittr_1.5        checkmate_2.0.0     cluster_2.0.7-1    
##  [19] mixtools_1.2.0      ROCR_1.0-7          R.utils_2.10.1     
##  [22] sandwich_2.5-1      colorspace_1.4-1    xfun_0.17          
##  [25] jsonlite_1.7.1      crayon_1.3.4        survival_2.41-3    
##  [28] zoo_1.8-8           iterators_1.0.12    ape_5.4-1          
##  [31] gtable_0.3.0        kernlab_0.9-29      prabclus_2.3-2     
##  [34] BiocGenerics_0.28.0 DEoptimR_1.0-8      scales_1.1.1       
##  [37] mvtnorm_1.1-1       bibtex_0.4.2.2      Rcpp_1.0.5         
##  [40] metap_1.4           dtw_1.21-3          plotrix_3.7-8      
##  [43] viridisLite_0.3.0   htmlTable_2.0.1     tmvnsim_1.0-2      
##  [46] reticulate_1.16     foreign_0.8-70      bit_4.0.4          
##  [49] proxy_0.4-24        mclust_5.4.6        SDMTools_1.1-221.2 
##  [52] Formula_1.2-3       tsne_0.1-3          stats4_3.5.1       
##  [55] htmlwidgets_1.5.1   httr_1.4.2          gplots_3.0.1.1     
##  [58] RColorBrewer_1.1-2  fpc_2.2-7           acepack_1.4.1      
##  [61] TFisher_0.2.0       modeltools_0.2-23   ellipsis_0.3.1     
##  [64] ica_1.0-2           farver_2.0.3        pkgconfig_2.0.3    
##  [67] R.methodsS3_1.8.1   flexmix_2.3-15      nnet_7.3-12        
##  [70] labeling_0.3        tidyselect_1.1.0    rlang_0.4.7        
##  [73] reshape2_1.4.4      munsell_0.5.0       tools_3.5.1        
##  [76] mathjaxr_1.0-1      ggridges_0.5.2      evaluate_0.14      
##  [79] stringr_1.4.0       yaml_2.2.1          knitr_1.29         
##  [82] bit64_4.0.5         fitdistrplus_1.1-1  robustbase_0.93-6  
##  [85] caTools_1.17.1.1    purrr_0.3.4         RANN_2.6.1         
##  [88] pbapply_1.4-3       nlme_3.1-137        R.oo_1.24.0        
##  [91] hdf5r_1.3.3         compiler_3.5.1      rstudioapi_0.11    
##  [94] png_0.1-7           tibble_3.0.3        stringi_1.5.3      
##  [97] lattice_0.20-35     multtest_2.38.0     vctrs_0.3.4        
## [100] mutoss_0.1-12       pillar_1.4.6        lifecycle_0.2.0    
## [103] BiocManager_1.30.10 Rdpack_1.0.0        lmtest_0.9-38      
## [106] bitops_1.0-6        irlba_2.3.3         gbRd_0.4-11        
## [109] R6_2.4.1            latticeExtra_0.6-28 renv_0.10.0        
## [112] KernSmooth_2.23-15  gridExtra_2.3       codetools_0.2-15   
## [115] MASS_7.3-50         gtools_3.8.2        assertthat_0.2.1   
## [118] rprojroot_1.3-2     withr_2.2.0         mnormt_2.0.2       
## [121] multcomp_1.4-13     diptest_0.75-7      parallel_3.5.1     
## [124] doSNOW_1.0.18       hms_0.5.3           grid_3.5.1         
## [127] rpart_4.1-13        class_7.3-14        rmarkdown_1.11     
## [130] segmented_1.2-0     Rtsne_0.15          git2r_0.27.1       
## [133] numDeriv_2016.8-1.1 Biobase_2.42.0      base64enc_0.1-3
```

</details>


***

<!-- END OF END MATTER -->
