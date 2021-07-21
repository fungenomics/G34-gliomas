---
title: "03 - Co-expression of astrocyte and interneuron programs in normal brain and tumor cells"
date: "21 July, 2021"
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

message("Cache: ", cache)
```

```
## Cache: ~/tmp/G34-gliomas/singlecell_normal/03/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/singlecell_normal/output/03
```

```
## G34-gliomas/singlecell_normal/figures/03//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Given that we a dual neuronal-astrocytic gene expression program in the tumors,
here we interrogate the expression of interneuron and astrocytic
gene signatures in the normal brain, and then evaluate the same signatures
in the tumors, to understand if this is an
aberrant state, or also observed under normal conditions.

For each dataset:

1. Take the mean expression each of an interneuron gene signature and an astrocyte signature
2. Subset to cells which are astrocytes, interneurons, or radial glia
3. Generate a scatterplot of cells based on their expression of the two signatures

# Libraries


```r
# Load liraries here
library(here)
library(tidyr)
library(dplyr)
library(glue)
library(ggplot2)
library(purrr)
library(readr)
library(Seurat)
```


# Load data

Human fetal brain (5-37 PCW, Nowakowski et al, Science, 2017):


```r
load(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski_seurat.Rda"))
```

Prepare the palette for this dataset, so that astrocytes, internurons, and
radial glial cells are distinct:


```r
palette_nowakowski <- read_csv(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski2017_tableS4_cluster_labels_with_colours.csv")) %>%
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
load(here("reference_datasets/2019_Jessa/Jessa2019_forebrain_seurat.Rda"))
```

Prepare palette:


```r
palette_atlas <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_Table_2a.tsv")) %>%
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
dge_le <- read_tsv(here("bulk_transcriptome_epigenome/output/02/fgsea_results_padj<0.01_IDH_top_le.tsv"))
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
load(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski_signatures.Rda"))

# Newborn interneuron signature
interneuron_genes    <- nowakowski_signatures$hg_sym$nIN1 
interneuron_genes_mm <- nowakowski_signatures$mm_sym$nIN1

# Astrocyte signature
astro_genes          <- nowakowski_signatures$hg_sym$Astrocyte
astro_genes_mm       <- nowakowski_signatures$mm_sym$Astrocyte
```


## Human fetal telencephalon (Nowakowski et al, 2017)


```r
data_nowakowski <- data.frame(
  interneuron_signature = cytobox::meanGeneExpression(nowa, interneuron_genes)$Mean_marker_expression,
  astrocyte_signature   = cytobox::meanGeneExpression(nowa, astro_genes)$Mean_marker_expression,
  cluster               = nowa@meta.data$WGCNAcluster,
  celltype              = nowa@meta.data$Cell_type
)

# Filter to cell types of interest for this analysis: astrocytes, inhibitory neurons
# and ganglionic eminence radial glia
data_nowakowski <- data_nowakowski %>% 
  filter(cluster %in% c("Astrocyte",
                        "nIN1", "nIN2", "nIN3", "nIN4",
                        "IN-CTX-CGE1", "IN-CTX-MGE1",
                        "MGE-IPC1", "MGE-IPC2", "MGE-IPC3",
                        "MGE-RG1", "MGE-RG2", "MGE-div"))

# Clean up
rm(nowa)
```

Generate the scatterplot:


```r
data_nowakowski %>% 
  rr_ggplot(aes(x = interneuron_signature, y = astrocyte_signature, colour = cluster), plot_num = 1) +
  geom_point(shape = 16, size = 1, alpha = 0.2) +
  scale_color_manual(values = palette_nowakowski, name = "Cell type") +
  geom_density_2d() +
  ggtitle("Human fetal telencephalon")
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/03//nowakowski_scatterplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/03//nowakowski_scatterplot...*]~</span>



## Developing mouse forebrain (Jessa et al, 2019)

In the mouse atlas:


```r
data_atlas <- data.frame(
  interneuron_signature = cytobox::meanGeneExpression(joint_cortex_small, interneuron_genes_mm)$Mean_marker_expression,
  astrocyte_signature   = cytobox::meanGeneExpression(joint_cortex_small, astro_genes_mm)$Mean_marker_expression,
  cluster               = gsub("_", " ", joint_cortex_small@meta.data$ID_20190730_with_blacklist_and_refined)
)
```

```
## NOTE: [Best3, 5830418K08Rik] undetected in Forebrain E12.5
```

```
## NOTE: [NA, Sirpb1c, Gm9733, Sirpb1a, Gm5150, B3galtl, Pmp2] undetected in Forebrain E12.5
```

```r
# Filter to relevant cell tpyes
data_atlas <- data_atlas %>% 
  filter(grepl("ASTR|VRGC|MGIN|CINH|INIP", cluster) & !grepl("BLACKLIST", cluster))

# Clean up
rm(joint_cortex_small)
```

Plot:


```r
data_atlas %>% 
  rr_ggplot(aes(x = interneuron_signature, y = astrocyte_signature, colour = cluster), plot_num = 1) +
  geom_point(shape = 16, size = 1, alpha = 0.2) +
  scale_color_manual(values = palette_atlas, name = "Cell type") +
  geom_density_2d() +
  ggtitle("Normal mouse brain")
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/03//mouse_forebrain_scatter-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/03//mouse_forebrain_scatter...*]~</span>

## H3G34R/V tumors

Load Seurat object with all samples combined, and re-constitute a V2 object from the V3 object:


```r
load("/lustre03/project/6004736/vlisi/fromHydra/SCRATCH/vlisi/LEGACY/G34/inHouseDataCombined/data/integrated_g34.seurat_v3.save_v2.Rda")

# get cell type projections used in Figure 5
load("/lustre03/project/6004736/vlisi/fromHydra/SCRATCH/vlisi/LEGACY/G34/allSamplesCombined/data/ScSnMetaDataWithNormals.RDat")
cell_meta <- ScSnDataWithNormals

# convert the object to V3
seurat_joint <- CreateSeuratObject(raw.data = integrated_g34@assays$RNA@counts,
                                   meta.data = integrated_g34@meta.data)
seurat_joint <- NormalizeData(seurat_joint)
seurat_joint <- SetDimReduction(seurat_joint, reduction.type = "umap", slot = "cell.embeddings",
                                integrated_g34@reductions[["umap"]]@cell.embeddings)
seurat_joint <- SetDimReduction(seurat_joint, reduction.type = "umap", slot = "key", "UMAP_")

# add in the cell type projections used in Figure 5
seurat_joint <- AddMetaData(seurat_joint, ScSnDataWithNormals[, c("ConsensusLineageRefinedConsensus"), drop = FALSE])
seurat_joint@meta.data$Cell_type_projection <- seurat_joint@meta.data$ConsensusHybrid2L3NP.smoothk10
seurat_joint@meta.data$Cell_type_broad_Fig5 <- seurat_joint@meta.data$ConsensusLineageRefinedConsensus
seurat_joint@meta.data$Status_normal_malignant <- seurat_joint@meta.data$CNVnormals

seurat_joint@meta.data <- seurat_joint@meta.data[, c(colnames(seurat_joint@meta.data)[1:22],
                                                     "Status_normal_malignant",
                                                     "Cell_type_projection",
                                                     "Cell_type_broad_Fig5")]

# clean up
rm(integrated_g34)
rm(ScSnDataWithNormals)

save(seurat_joint, file = glue("{out}/integrated_seurat.Rda"))
```


```r
load(glue("{out}/integrated_seurat.Rda"))
```

Subset the same signatures used in the normal brain to the genes detected in the tumors:


```r
interneuron_genes_keep <- base::intersect(interneuron_genes, rownames(seurat_joint@data))
astro_genes_keep <- base::intersect(astro_genes, rownames(seurat_joint@data))
```

Calculate mean expression:


```r
data_tumors <- data.frame(
  interneuron_signature = cytobox::meanGeneExpression(seurat_joint, interneuron_genes_keep)$Mean_marker_expression,
  astrocyte_signature   = cytobox::meanGeneExpression(seurat_joint, astro_genes_keep)$Mean_marker_expression,
  cluster_broad         = seurat_joint@meta.data$Cell_type_broad_Fig5
)
```


Plot:


```r
data_tumors %>% 
  filter(cluster_broad %in% c("neurons", "astro")) %>% 
  rr_ggplot(aes(x = interneuron_signature, y = astrocyte_signature, colour = cluster_broad), plot_num = 1) +
  geom_point(shape = 16, size = 1, alpha = 0.2) +
  scale_color_manual(values = c("astro" = "#FF3F38", "neurons" = "#084081"), name = "Cell type") +
  geom_density_2d() +
  ggtitle("G34R/V HGGs")
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/singlecell_normal/figures/03//tumor_scatter-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/singlecell_normal/figures/03//tumor_scatter...*]~</span>



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2021-07-21 10:30:46
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [677ad1b] 2020-12-11: Update & clean up - Update paths - Add pointers to figures in the README - Output final set of signatures used
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] Seurat_2.3.4  Matrix_1.2-14 cowplot_0.9.4 readr_1.3.1   purrr_0.3.4  
##  [6] ggplot2_3.1.0 glue_1.4.2    dplyr_0.8.0   tidyr_0.8.2   here_0.1     
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
##  [25] crayon_1.3.4        jsonlite_1.7.1      survival_2.41-3    
##  [28] zoo_1.8-8           iterators_1.0.12    ape_5.4-1          
##  [31] gtable_0.3.0        kernlab_0.9-29      prabclus_2.3-2     
##  [34] BiocGenerics_0.28.0 DEoptimR_1.0-8      scales_1.1.1       
##  [37] mvtnorm_1.1-1       bibtex_0.4.2.2      Rcpp_1.0.5         
##  [40] metap_1.4           dtw_1.21-3          plotrix_3.7-8      
##  [43] viridisLite_0.3.0   htmlTable_2.0.1     tmvnsim_1.0-2      
##  [46] reticulate_1.16     foreign_0.8-70      bit_4.0.4          
##  [49] proxy_0.4-24        mclust_5.4.6        SDMTools_1.1-221.2 
##  [52] Formula_1.2-3       stats4_3.5.1        tsne_0.1-3         
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
##  [85] caTools_1.17.1.1    randomForest_4.6-14 RANN_2.6.1         
##  [88] pbapply_1.4-3       nlme_3.1-137        cytobox_0.6.1      
##  [91] R.oo_1.24.0         hdf5r_1.3.3         compiler_3.5.1     
##  [94] rstudioapi_0.11     png_0.1-7           tibble_3.0.3       
##  [97] stringi_1.5.3       lattice_0.20-35     multtest_2.38.0    
## [100] vctrs_0.3.4         mutoss_0.1-12       pillar_1.4.6       
## [103] lifecycle_0.2.0     Rdpack_1.0.0        lmtest_0.9-38      
## [106] data.table_1.13.0   bitops_1.0-6        irlba_2.3.3        
## [109] gbRd_0.4-11         R6_2.4.1            latticeExtra_0.6-28
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
