---
title: "01 - Bulk RNAseq"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "09 July, 2020"
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
```

</details>

The root directory of this project is:

```
## /mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas-repo/bulk/output/01
```

```
## G34-gliomas-repo/bulk/figures/01//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Here we perform basic RNAseq analysis for G34R/V mutant high-grade gliomas.

This document:

* loads metadata for bulk RNAseq of G34R/V mutant tumors as well as tumors of
other comparator groups
* loads metadata for bulk RNAseq of tumor-derived cell lines
* prepares the inputs for our in-house bulk RNAseq pipeline (code not included in this
repository)
* generates basic plots of gene expression levels

More complex downstream analysis is performed in subsequent documents.


# Libraries


```r
# General
library(tidyverse)
library(magrittr)
library(glue)

# For plotting
library(ggplot2)
library(scales)
library(ggrepel)
library(cytobox)
ggplot2::theme_set(cytobox::theme_min(base_size = 13))

# Custom
source(here::here("bulk/analysis/functions.R"))
```

# Analysis

## Prepare inputs for RNAseq pipeline - cell lines

Load sample metadata:


```r
bulk_samples <- readxl::read_xlsx(
  here::here("bulk", "data", "2020-02-25_bulkRNAseq_metadata.xlsx")) %>% 
  mutate(G34_batch = ifelse(is.na(G34_batch), "Old", G34_batch))

# Count number of samples per group
table(bulk_samples$Group)
```

```
## 
##             HGG-G34R/V     HGG-H3.3-K27M-pons HGG-H3.3-K27M-thalamic 
##                     18                      8                      5 
##                HGG-IDH                 HGG-WT        PNET-HGNET-BCOR 
##                     10                      9                      5
```

Next, we prepare the inputs for the RNAseq. For each comparison (i.e. differential
gene expression), we need two files:

* `info.samples.tsv` - describing the membership of each sample to each group
* `info.groups.tsv`  - describing the two groups in the comparison


We perform differential gene expression analysis between G34 tumors and each other
“comparator” group:


```r
# Get the tumor groups (other than G34R/V)
(groups <- unique(bulk_samples$Group) %>% tail(5))
```

```
## [1] "HGG-H3.3-K27M-pons"     "HGG-H3.3-K27M-thalamic"
## [3] "HGG-IDH"                "HGG-WT"                
## [5] "PNET-HGNET-BCOR"
```

```r
# Generate inputs for each group
for (group in groups) {
  
  # Create a name describing the comparison
  comp <- paste0(gsub("/", ".", "HGG-G34R/V"), "_vs_", gsub("/", ".", group))
  dir.create(file.path(out, comp), showWarnings = FALSE)
  
  # Create the info.samples.tsv, which needs three columns: ID (path
  # to the alignment and counts on Hydra), Nickname (short sample name), and Group
  bulk_samples %>%
    select(ID, Nickname = Sample, Group, Batch = G34_batch) %>%
    filter(Group %in% c("HGG-G34R/V", group)) %>%
    write_tsv(file.path(out, comp, "info.samples.tsv"))
  
  # Create the info.groups.tsv
  data.frame(Group = c("HGG-G34R/V", group),
             Label = c("HGG-G34R/V", group),
             Order = c(2, 1),
             stringsAsFactors = FALSE) %>%
    write_tsv(file.path(out, comp, "info.groups.tsv"))
  
}
```


We will also need one version with all samples together, since normalization is
dataset-dependent, and we will compare expression levels across groups.


```r
# As above, but without batch correction since we won't be doing differential
# expression analysis here
dir.create(file.path(out, "All_tumor_samples"), showWarnings = FALSE)

bulk_samples %>%
  select(ID, Nickname = IC_Sample, Group) %>%
  write_tsv(file.path(out, "All_tumor_samples", "info.samples.tsv"))

data.frame(Group = unique(bulk_samples$Group),
           Label = unique(bulk_samples$Group),
           Order = seq_along(unique(bulk_samples$Group)),
           stringsAsFactors = FALSE) %>%
  write_tsv(file.path(out, "All_tumor_samples", "info.groups.tsv"))
```


Finally, we peform one analysis among G34R/V tumors, between PDGFRA mutants
vs. WT:


```r
# Exactly as above
dir.create(file.path(out, "G34RV_by_PDGFRA"), showWarnings = FALSE)

bulk_samples %>%
  filter(Group == "HGG-G34R/V") %>% 
  select(ID, Nickname = Sample, Group = Genotype_PDGFRA, G34_batch, Genotype_G34) %>%
  write_tsv(file.path(out, "G34RV_by_PDGFRA", "info.samples.tsv"))

data.frame(Group = c("Mutant", "WT"),
           Label = c("Mutant", "WT"),
           Order = c(2, 1),
           stringsAsFactors = FALSE) %>%
  write_tsv(file.path(out, "G34RV_by_PDGFRA", "info.groups.tsv"))
```

The in-house RNAseq pipeline is then run on these inputs, generating
normalized expression counts as well as the output of differential gene expression
analysis using DESeq2.


## Prepare inputs for RNAseq pipeline - cell lines



## Plot gene expression

Prepare colour palettes for visualization:


```r
palette_groups <- c(
  "G34R/V"            = "cyan4",
  "G34R/V-PDGFRA mut" = "#26A27D",
  "G34R/V-PDGFRA WT"  = "#145448",
  "K27M-pons"         = "gold1",
  "K27M-thal"         = "goldenrod3",
  "IDH"               = "#883067",
  "WT"                = "navy",
  "HGNET-BCOR"        = "lavender")

palette_pdgfra <- c("Mutant" = "red",
                    "WT"     = "gray70")

palette_g34    <- c("G34R"   = "#3d8a71",
                    "G34V"   = "#6bbfad") 
```



Load normalized counts for genes of interest from the pipeline run:


```r
pipeline_path <- "../../../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/"

# Define genes which we'll plot in this document
genes_of_interest <- c("GSX2", "DLX1", "DLX2", "PDGFRA", "MOXD1", "EOMES", "NEUROD2")

# Load counts and convert to tidy format
counts_subset <- extract_pipeline_counts(path = file.path(pipeline_path,
                                                          "All_tumor_samples/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                         goi = genes_of_interest) %>% 
  left_join(bulk_samples, by = c("sample" = "Sample"))
```

```
## Joining, by = "gene_symbol"
```

```r
# Tidy up and simplify the names of the tumor groups
counts_subset <- counts_subset %>%
  mutate(Group2 = recode(Group,
                         "HGG-G34R/V" = "G34R/V",
                         "HGG-H3.3-K27M-pons" = "K27M-pons",
                         "HGG-H3.3-K27M-thalamic" = "K27M-thal",
                         "HGG-IDH" = "IDH",
                         "HGG-WT" = "WT",
                         "PNET-HGNET-BCOR" = "HGNET-BCOR")) %>%
  mutate(Group3 = case_when(
    Genotype_PDGFRA == "Mutant" ~ paste0(Group2, "-", "PDGFRA mut"),
    Genotype_PDGFRA == "WT" ~ paste0(Group2, "-", "PDGFRA WT"),
    TRUE ~ Group2
  )) %>% 
  mutate(Group3 = factor(Group3, levels = c(c("G34R/V-PDGFRA mut",
                                              "G34R/V-PDGFRA WT",
                                              "IDH",
                                              "K27M-pons",
                                              "K27M-thal",
                                              "WT",
                                              "HGNET-BCOR"))))
```

### PDGFRA

Plot PDGFRA expression levels across tumor groups:


```r
# Plot PDGFRA expression and save source data alongside
counts_subset %>% 
  filter(gene_symbol == "PDGFRA") %>% 
  select(sample, group = Group3, Pdgfra_normalized_expression = gene_expression) %>% 
  rr_ggplot(plot_num = 1,
            aes(x = group, y = Pdgfra_normalized_expression)) +
  geom_boxplot(aes(fill = group), width = 0.5) +
  scale_fill_manual(values = palette_groups) +
  geom_point() +
  scale_y_continuous(labels = comma) +
  rotateX() +
  ylab("Normalized expression") + xlab("Group") + ggtitle("PDGFRA") +
  noLegend()
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/01//bulk_RNAseq_PDGFRA_expression-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/01//bulk_RNAseq_PDGFRA_expression...*]~</span>

### Lineage specific TFs

Plot expression levels of certain interneuron and excitatory neuron
lineage specific TFs/genes:



```r
# Subset data to save as source data
boxplot_genes <- c("MOXD1", "GSX2", "DLX1", "DLX2", "EOMES", "NEUROD2")

counts_subset_boxplot <- counts_subset %>% 
  filter(gene_symbol %in% boxplot_genes) %>% 
  select(sample, gene_symbol, group = Group2, Normalized_expression = gene_expression)

# Helper function
boxplot_tumor_groups <- function(gene) {
  
  p1 <- counts_subset_boxplot %>%
    # Specify the order of tumour groups along the x-axis
    mutate(group = factor(group, levels = c(c("G34R/V", "IDH", "K27M-pons", "K27M-thal", "WT", "HGNET-BCOR")))) %>% 
    filter(gene_symbol == gene) %>%
    ggplot(aes(x = group, y = Normalized_expression)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA) +
    geom_point(size = 1.5) +
    scale_fill_manual(values = palette_groups) +
    rotateX() +
    ylab("Normalized expression") + xlab("Group") + ggtitle(gene) +
    noLegend()
  
  p1
  
}

# Loop over genes, generate boxplot, and assemble into one figure
boxplots <- map(boxplot_genes, boxplot_tumor_groups)
rr_plot_grid(df = counts_subset_boxplot, plot_num = 1, plotlist = boxplots, ncol = 3)
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/01//bulk_RNAseq_TF_expression-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/01//bulk_RNAseq_TF_expression...*]~</span>




<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-07-09 08:18:52
```

The git repository and last commit:

```
## Local:    master /mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [c2ba647] 2020-06-30: Update rr infrastructure to handle analyses in subdirectories
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
## BLAS/LAPACK: /var/chroots/hydrars-centos-7/usr/lib64/R/lib/libRblas.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C         LC_TIME=C           
##  [4] LC_COLLATE=C         LC_MONETARY=C        LC_MESSAGES=C       
##  [7] LC_PAPER=C           LC_NAME=C            LC_ADDRESS=C        
## [10] LC_TELEPHONE=C       LC_MEASUREMENT=C     LC_IDENTIFICATION=C 
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] cowplot_0.9.4   bindrcpp_0.2.2  cytobox_0.6.1   ggrepel_0.8.0  
##  [5] scales_1.0.0    glue_1.4.1      magrittr_1.5    forcats_0.3.0  
##  [9] stringr_1.3.1   dplyr_0.7.7     purrr_0.3.4     readr_1.3.1    
## [13] tidyr_0.8.2     tibble_3.0.1    ggplot2_3.1.0   tidyverse_1.2.1
## [17] here_0.1       
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.2.0        snow_0.4-3          backports_1.1.3    
##   [4] Hmisc_4.2-0         plyr_1.8.4          igraph_1.2.2       
##   [7] lazyeval_0.2.1      splines_3.5.0       digest_0.6.16      
##  [10] foreach_1.4.4       htmltools_0.3.6     viridis_0.5.1      
##  [13] lars_1.2            gdata_2.18.0        checkmate_1.9.1    
##  [16] cluster_2.0.7-1     mixtools_1.1.0      ROCR_1.0-7         
##  [19] modelr_0.1.3        R.utils_2.7.0       colorspace_1.4-0   
##  [22] rvest_0.3.2         haven_2.0.0         xfun_0.12          
##  [25] crayon_1.3.4        jsonlite_1.6        bindr_0.1.1        
##  [28] survival_2.41-3     zoo_1.8-4           iterators_1.0.10   
##  [31] ape_5.2             gtable_0.2.0        kernlab_0.9-27     
##  [34] prabclus_2.2-7      DEoptimR_1.0-8      mvtnorm_1.0-10     
##  [37] bibtex_0.4.2        Rcpp_1.0.0          metap_1.1          
##  [40] dtw_1.20-1          viridisLite_0.3.0   htmlTable_1.13.1   
##  [43] reticulate_1.10     foreign_0.8-70      bit_1.1-14         
##  [46] proxy_0.4-22        mclust_5.4.2        SDMTools_1.1-221   
##  [49] Formula_1.2-3       tsne_0.1-3          stats4_3.5.0       
##  [52] htmlwidgets_1.3     httr_1.4.0          gplots_3.0.1.1     
##  [55] RColorBrewer_1.1-2  fpc_2.1-11.1        acepack_1.4.1      
##  [58] modeltools_0.2-22   ellipsis_0.2.0.1    Seurat_2.3.4       
##  [61] ica_1.0-2           pkgconfig_2.0.2     R.methodsS3_1.7.1  
##  [64] flexmix_2.3-14      nnet_7.3-12         labeling_0.3       
##  [67] reshape2_1.4.3      tidyselect_1.1.0    rlang_0.4.6        
##  [70] munsell_0.5.0       cellranger_1.1.0    tools_3.5.0        
##  [73] cli_1.0.1           generics_0.0.2      broom_0.5.1        
##  [76] ggridges_0.5.1      evaluate_0.12       yaml_2.2.0         
##  [79] npsurv_0.4-0        knitr_1.21          bit64_0.9-7        
##  [82] fitdistrplus_1.0-14 robustbase_0.93-2   caTools_1.17.1.1   
##  [85] randomForest_4.6-14 RANN_2.6            pbapply_1.4-0      
##  [88] nlme_3.1-137        R.oo_1.22.0         xml2_1.2.0         
##  [91] hdf5r_1.0.0         compiler_3.5.0      rstudioapi_0.9.0   
##  [94] png_0.1-7           lsei_1.2-0          stringi_1.2.4      
##  [97] lattice_0.20-35     trimcluster_0.1-2.1 Matrix_1.2-14      
## [100] vctrs_0.3.1         pillar_1.4.4        lifecycle_0.2.0    
## [103] Rdpack_0.10-1       lmtest_0.9-36       data.table_1.12.0  
## [106] bitops_1.0-6        irlba_2.3.3         gbRd_0.4-11        
## [109] R6_2.3.0            latticeExtra_0.6-28 KernSmooth_2.23-15 
## [112] gridExtra_2.3       codetools_0.2-15    MASS_7.3-49        
## [115] gtools_3.8.1        assertthat_0.2.0    rprojroot_1.3-2    
## [118] withr_2.1.2         diptest_0.75-7      parallel_3.5.0     
## [121] doSNOW_1.0.16       hms_0.4.2           grid_3.5.0         
## [124] rpart_4.1-13        class_7.3-14        rmarkdown_1.11     
## [127] segmented_0.5-3.0   Rtsne_0.15          git2r_0.27.1       
## [130] lubridate_1.7.4     base64enc_0.1-3
```

</details>


***

<!-- END OF END MATTER -->
