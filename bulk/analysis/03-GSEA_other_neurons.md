---
title: "03 - GSEA with cell type gene sigantures"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "13 July, 2020"
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
## G34-gliomas-repo/bulk/output/03
```

```
## G34-gliomas-repo/bulk/figures/03//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

In this document, we gather a reference of developmental gene signatures from several brain cell types and studies, and evaluate their enrichment in the differentially-expressed genes between G34 tumors and other groups. The enrichment analysis is performed with GSEA, using the package fgsea.

We obtained reference datasets for the forebrain across species and the lifespan:

* Jessa et al, 2019
* Nowakowski et al, 2017
* Mayer et al, 2018
* Velmeshev et al, 2019

We also obtained two reference datasets capturing the adult SVZ:

* Anderson et al, 2020
* Mizrak et al, 2019

# Libraries


```r
# General
library(tidyverse)
library(magrittr)
library(glue)

# For analysis
library(fgsea)

# For plotting
library(ggplot2)
library(scales)
library(cytobox)
library(pheatmap)
library(ggrepel)
library(cowplot)
ggplot2::theme_set(cytobox::theme_min(base_size = 13))

# Custom
source(here::here("bulk/analysis/functions.R"))
```

# Load data

This sections loads the reference datasets, including cell-type specific gene
signatures, corresponding palettes, and cluster annotations.

## Fetal mouse brain

Load signatures from our Nature Genetics paper:


```r
# Annotation
jessa_table_2a <- read_tsv(here("bulk/data/Jessa2019_Table_2a.tsv")) %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  mutate(Age = paste0("Mouse ", Age))

# Get data from paper, filter to clusters with a valid signature, from mouse samples
jessa_table_2a_mouse_noblacklist <- jessa_table_2a %>% 
  filter(Species == "Mouse" & Structure == "Forebrain" & !is.na(Signature))

# Parse the signature column to get gene signatures as vectors
atlas_signatures <- list("mm_sym" = map(jessa_table_2a_mouse_noblacklist$Signature,
                                        ~ stringr::str_split(.x, ", ") %>%
                                          # Get the first element, otherwise we have
                                          # a list of one-element lists
                                          .[[1]]))
names(atlas_signatures$mm_sym) <- jessa_table_2a_mouse_noblacklist$Cluster

# Convert to human symbols / Ensembl IDs
atlas_signatures$hg_sym <- map(atlas_signatures$mm_sym,
                               # Filter out NAs and empty strings
                               ~ .x %>% cytobox::mm2hg() %>% .[!is.na(.)] %>% .[. != ""])
atlas_signatures$hg_ens <- map(atlas_signatures$hg_sym,
                               ~ .x %>% cytobox::symbols2ensembl() %>% .[!is.na(.)] %>% .[. != ""])

rr_saveRDS(object = atlas_signatures,
           file = glue("{out}/atlas_signatures.Rds"),
           desc = "List of gene signatures from Jessa et al 2019, with mouse and human gene symbols and ENSEMBL IDs")
```

## Fetal human brain

We've previously extracted gene signatures from this study (provided in Table S5) and prepared them in a
similar format:


```r
# Annotation
nowakowski_cluster_anno <- read_csv(here("bulk/data/Nowakowski2017_Table_4_cluster_labels_with_colours.csv")) %>%
  select(Cluster = `Cluster Name`,
         Cell_type = `Cluster Interpretation`) %>% 
  mutate(Cluster = paste0("HF ",  Cluster)) %>% 
  mutate(Age = "Human fetal",
         Sample = "Human fetal cortex/MGE")
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
# Signatures
load(here("bulk/data/Nowakowski2017_signatures.Rda"))

# Denote as human fetal
names(nowakowski_signatures$hg_ens) <- paste0("HF ", names(nowakowski_signatures$hg_ens))
names(nowakowski_signatures$hg_sym) <- paste0("HF ", names(nowakowski_signatures$hg_sym))
```


## Pediatric/adult human brain


```r
# Annotation
velmeshev_anno <- read_tsv(here("bulk/data/Velmeshev2019_cluster_names.tsv")) %>% 
  mutate(Age = "Human ped/adult",
         Cluster = paste0("HP ", Cluster))
```

```
## Parsed with column specification:
## cols(
##   Cluster = col_character(),
##   Cell_type = col_character()
## )
```

```r
# Signatures
velmeshev_signatures <- readRDS(here("bulk/data/Velmeshev2019_signatures.Rda"))

# Denote as human postnatal/pediatric
names(velmeshev_signatures$hg_ens) <- paste0("HP ", names(velmeshev_signatures$hg_ens))
names(velmeshev_signatures$hg_sym) <- paste0("HP ", names(velmeshev_signatures$hg_sym))
```

## Mouse ganglionic eminences


```r
# Annotation
mayer_anno <- bind_rows(
  read_tsv(here("bulk/data/Mayer2018_10x_dataset_cluster_sample_labels.tsv"))     %>% mutate(Cluster2 = paste0("GEs_Ctx ", Cluster)),
  read_tsv(here("bulk/data/Mayer2018_dropseq_dataset_cluster_sample_labels.tsv")) %>% mutate(Cluster2 = paste0("GEs ", Cluster))
) %>% 
  mutate(Cell_type = paste0(Sample, " ", Cluster),
         Age = ifelse(is.na(Age), "Mouse E13-14", paste0("Mouse ", Age))) %>%
  select(Cluster = Cluster2, Cell_type, Sample, Age)
```

```
## Parsed with column specification:
## cols(
##   Cluster = col_double(),
##   Sample = col_character(),
##   Age = col_character()
## )
```

```
## Parsed with column specification:
## cols(
##   Cluster = col_double(),
##   Sample = col_character()
## )
```

```r
# Signatures
mayer10x_signatures     <- readRDS(here("bulk/data/Mayer2018_10x_signatures.Rds"))
mayerdropseq_signatures <- readRDS(here("bulk/data/Mayer2018_dropseq_signatures.Rds"))

mayer10x_signatures     <- map(mayer10x_signatures, ~ .x %>% set_names(paste0("GEs_Ctx ", names(.x))))
mayerdropseq_signatures <- map(mayerdropseq_signatures, ~ .x %>% set_names(paste0("GEs ", names(.x))))
```


## Adult V-SVZ


```r
# Signatures
mizrak_rep1_signatures <- readRDS(here("bulk/data/Mizrak2019_Rep1_signatures.Rds"))
mizrak_rep2_signatures <- readRDS(here("bulk/data/Mizrak2019_Rep2_signatures.Rds"))

# Add a prefix to the cluster names
mizrak_rep1_signatures <- map(mizrak_rep1_signatures, ~ .x %>% set_names(paste0("VSVZ-1 ", names(.x))))
mizrak_rep2_signatures <- map(mizrak_rep2_signatures, ~ .x %>% set_names(paste0("VSVZ-2 ", names(.x))))

mizrak_anno <- data.frame(Cluster   = c(names(mizrak_rep1_signatures$hg_sym),
                                        names(mizrak_rep2_signatures$hg_sym)),
                          Cell_type = c(names(mizrak_rep1_signatures$hg_sym),
                                        names(mizrak_rep2_signatures$hg_sym)),
                          Age = "Mouse adult",
                          Sample = "V-SVZ")
```


## Striatal SVZ


```r
# Signatures
anderson_signatures <- readRDS(here("bulk/data/Anderson2020_signatures.Rds"))
anderson_signatures <- map(anderson_signatures, ~ .x %>% set_names(paste0("Str/SVZ ", names(.x))))

# Annotation
anderson_anno <- data.frame(Cluster = names(anderson_signatures$hg_sym),
                            Cell_type = names(anderson_signatures$hg_sym),
                            Age = "Mouse P9",
                            Sample = "Striatum")
```

## Gather all signatures


```r
# Human gene symbols
signatures_sym <- c(atlas_signatures$hg_sym,
                    nowakowski_signatures$hg_sym,
                    velmeshev_signatures$hg_sym,
                    mayer10x_signatures$hg_sym,
                    mayerdropseq_signatures$hg_sym,
                    anderson_signatures$hg_sym,
                    mizrak_rep1_signatures$hg_sym,
                    mizrak_rep2_signatures$hg_sym)

# Human Ensembl IDs
signatures_ens <- c(atlas_signatures$hg_ens,
                    nowakowski_signatures$hg_ens,
                    velmeshev_signatures$hg_ens,
                    mayer10x_signatures$hg_ens,
                    mayerdropseq_signatures$hg_ens,
                    anderson_signatures$hg_ens,
                    mizrak_rep1_signatures$hg_ens,
                    mizrak_rep2_signatures$hg_ens)

rr_saveRDS(file = glue("{out}/signatures_ens.Rds"),
           desc = "Cell type specific gene signatures from several forebrain references (list of vectors, using Ensemble gene IDs)",
           signatures_ens)
```

```
## ...writing description of signatures_ens.Rds to G34-gliomas-repo/bulk/output/03/signatures_ens.desc
```


## Assemble cluster annotations

Put together the cell type annotations for each of the datasets, since they are often
labelled by abbreviations. Create uniform dataframes, with at least four columns:
Sample, Cluster, Cell_type, and Age. The Cluster column matches the names of
the signatures in the `signatures_*` objects created above.


```r
cell_type_anno <- bind_rows(jessa_table_2a %>% select(Sample, Age, Cell_type, Cluster),
                            nowakowski_cluster_anno,
                            velmeshev_anno,
                            mizrak_anno,
                            mayer_anno,
                            anderson_anno)
```


## Palettes

Create colour palettes for each dataset, for visualization:


```r
palette_atlas <- jessa_table_2a %>% select(Cluster, Colour) %>% deframe()

palette_nowakowski <- read_csv(here("bulk/data/Nowakowski2017_Table_4_cluster_labels_with_colours.csv")) %>%
  mutate(Cluster = paste0("HF ", `Cluster Name`)) %>%
  select(Cluster, colour = Colour) %>%
  deframe()
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
palette_type <- c("RGC" = "#ffcc00",
                  "RGC (prolif.)" = "#e8a805",
                  "Oligodendrocyte precursors" = "#e0de53",
                  "Oligodendrocytes" = "#b4e04e",
                  "Astrocytes" = "#00a385",
                  "Neuronal progenitors" = "#ffbda3",
                  "Prenatal inhib. neurons" = "#8798C8",
                  "Pediatric/adult inhib. neurons" = "#135ca0",
                  "Other neurons" = "darkred",
                  "Other" = "gray90",
                  "Immune" = "#aca2b2",
                  "Glial progenitors" = "#d5d98b",
                  "Choroid/ependymal" = "#8ee5cf",
                  "Non-neuroectoderm" = "gray70")

palette_age <- c("Human fetal" = "#458cde",
                 "Mouse E12.5" = "#73a8e6",
                 "Mouse E13-14" = "#4b7db8",
                 "Mouse E13.5" = "#4b7db8",
                 "Mouse E15.5" = "#a4c6ed",
                 "Mouse E18.5" = "#6889b0",
                 "Mouse P0" = "#f2b39d",
                 "Mouse Embryonic/Postnatal" = "#eba288",
                 "Mouse P3" = "tomato",
                 "Mouse P6" = "red",
                 "Mouse P9" = "#b01010",
                 "Mouse P10" = "#b01010",
                 "Human ped/adult" = "darkred",
                 "Mouse adult" = "#4a0303")
```

# Analysis

With the cell type signatures derived from the various scRNAseq studies, we can now
use these as input to GSEA, and evaluate their enrichment in differentially expressed
genes between G34R/V tumors vs. other pediatric brain tumor types.

## Prep stats


```r
pipeline_path <- "../../../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-07_other_neuronal_tumors/"

# Keep all the DGE analyses between tumor groups, accounting for batch
comps <- list.files(pipeline_path, full.names = TRUE)
names(comps) <- basename(comps)

# Get the file with the DESeq2 DGE output, and create a vector with the Ensembl
# gene IDs and stat for that gene, in each comparison
stats <- map(comps, ~ read_tsv(Sys.glob(file.path(.x, "diff/Ensembl.ensGene.exon/*.tsv"))) %>%
               filter(!is.na(stat)) %>%
               separate(ID, into = c("Ens", "Sym"), sep = ":") %>%
               select(Ens, stat) %>%
               deframe())

rr_saveRDS(file = glue("{out}/stats.Rds"),
           desc = "List of DESeq2 stats from differential expression analyses between tumor groups, to use for GSEA",
           stats)
```

## Run GSEA

Running GSEA analysis for each differential expression comparison, using the GSVA package:


```r
bulk_fgsea <- map(stats, ~ fgsea(pathways = signatures_ens,
                                 stats   = .x,
                                 # Default params
                                 minSize = 15,
                                 maxSize = 500,
                                 nperm   = 10000))

rr_saveRDS(file = glue("{out}/bulk_fgsea.Rds"),
           desc = "List of results from fgsea for enrichment of cell type signatures in tumor differential expression analyses",
           object = bulk_fgsea)
```

```
## ...writing description of bulk_fgsea.Rds to G34-gliomas-repo/bulk/output/03/bulk_fgsea.desc
```

## Filter gene signatures

For the forebrain reference, remove clusters which are proliferating based on either labeled as such, or G2M score,
as cell-cycle genes may confound GSEA enrichment scores particularly in these tumor comparisons.
I am not removing any clusters from the striatum/SVZ reference.


```r
# 1. Atlas
# Load signatures with cell cycle scores, to quantitatively set a filter
infosig_with_cc <- read_tsv(here("bulk/data/Jessa2019_cluster_cell_cycle_scores.tsv")) %>%
  arrange(desc(Median_G2M))
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
##   Colour = col_character()
## )
```

```
## See spec(...) for full column specifications.
```

```r
infosig_with_cc %>%
  filter(Species == "Mouse") %>%
  mutate(Cluster = factor(Cluster, levels = .$Cluster)) %>%
  ggplot(aes(x = Cluster, y = Median_G2M)) +
  geom_point(aes(colour = Cluster)) +
  scale_colour_manual(values = palette_atlas) +
  # Set an arbitrary threshold, where the G2M score starts tapering off
  geom_hline(yintercept = 0.35) +
  theme_min() +
  rotateX() +
  noLegend()
```

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/03//filter_signatures-1.png)<!-- -->

```r
atlas_sigs_drop <- infosig_with_cc %>%
  # Remove anything with G2M scores above the threshold, OR which is labelled
  # as proliferating, as indicated by a trailing "-P" in the cluster abbreviation
  filter(Median_G2M > 0.35 | grepl("-P$", Cluster) |
           # Remove anything labeled Other/Other neurons
           (Cell_class %in% c("Other", "Other neurons")) |
           # Filter to only forebrain signatures
           (Structure != "Forebrain")) %>%
  pull(Cluster)

# 2. Nowakowski
nowakowski_sigs_drop <- nowakowski_cluster_anno %>%
  # Filter out proliferating clusters, 
  filter(grepl("-P$|div|G2M/S|MGE-IPC1|MGE-IPC2|prolif", Cluster) |
           # Filter out Unknown
           grepl("Unknown", Cell_type)) %>%
  pull(Cluster)

bulk_fgsea_filt <- map(bulk_fgsea, ~ .x %>%
                         filter(!(pathway %in% c(atlas_sigs_drop, nowakowski_sigs_drop))))
```

## Tidy output


```r
# A function to convert the leading edge from ENSG id to symbols
flatten_leading_edge <- function(x) {
  
  # Convert leading edge genes to gene symbols
  for (i in 1:nrow(x)) {
    
    x[i, ]$leadingEdge <- x[i, ]$leadingEdge %>% ensembl2symbols_safe() %>% glue_collapse(sep = ", ")
    
  }
  
  x$leadingEdge <- as.character(x$leadingEdge)
  x
  
}

fgsea_df <- imap_dfr(bulk_fgsea_filt, ~ .x %>%
                       as.data.frame() %>%
                       flatten_leading_edge() %>%
                       tibble::add_column(Comparison = .y, .before = 1) %>%
                       arrange(desc(NES))) %>%
  rename(Signature = pathway) %>%
  left_join(cell_type_anno, by = c("Signature" = "Cluster")) %>%
  select(Comparison, Signature, Sample, Age, Cell_type, everything())

# Save output
rr_write_tsv(df = fgsea_df,
             path = glue("{out}/fgsea_df.tsv"),
             desc = "Tidied dataframe with GSEA results for bulk tumor group comparisons")
```

```
## ...writing description of fgsea_df.tsv to G34-gliomas-repo/bulk/output/03/fgsea_df.desc
```


## Heatmap across comparisons

Helpers for creating the heatmaps:


```r
load(here("bulk/output/02/signif_in_idh.Rda"))

# Helper function to actually generate the heatmap
hm_fun <- partial(pheatmap,
                  border_color = NA,
                  color = colorRampPalette(c("blue","white","red"))(100),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  cellwidth = 12,
                  cellheight = 12,
                  number_color = "white",
                  # Prepare the cluster annotation tracks
                  annotation_colors = list("Age"  = palette_age,
                                           "Type" = palette_type,
                                           "Species" = c("Human" = "black", "Mouse" = "gray70")))

# Include only the high-grade gliomas at first
row_order <- c("ETMR_vs_HGG-IDH",
               "NB-FOXR2_vs_HGG-IDH",
               "HGG-G34R.V_vs_ETMR",
               "HGG-G34R.V_vs_NB-FOXR2")
```

### Forebrain reference


```r
# Prep input
heatmap_inputs_fb <- prep_gsea_heatmap(fgsea_df,
                                       signatures = signif_in_idh,
                                       filters = quos(!(Sample %in% c("Striatum", "V-SVZ"))),
                                       row_order = row_order)

# Generate heatmap
hm_fun(mat = heatmap_inputs_fb$heatmap_data_wide,
       annotation_col = heatmap_inputs_fb$col_anno,
       display_number = heatmap_inputs_fb$signif_data_wide,
       main = "Forebrain references",
       filename = glue("{figout}/gsea_heatmap_hgg.pdf"))

hm_fun(mat = heatmap_inputs_fb$heatmap_data_wide,
       annotation_col = heatmap_inputs_fb$col_anno,
       display_number = heatmap_inputs_fb$signif_data_wide,
       main = "Forebrain references",
       filename = glue("{figout}/gsea_heatmap_hgg.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg.png"))
```

<img src="/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/03///gsea_heatmap_hgg.png" width="4599" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/03//gsea_heatmap_hgg...*]~</span>



### Striatal SVZ

Generate one heatmap for the striatal SVZ signatures only:


```r
heatmap_inputs_str <- prep_gsea_heatmap(fgsea_df,
                                        signatures = signif_in_idh,
                                        filters = quos(Sample == "Striatum"),
                                        row_order = row_order)

# Generate heatmap
hm_fun(mat = heatmap_inputs_str$heatmap_data_wide,
       annotation_col = heatmap_inputs_str$col_anno,
       display_number = heatmap_inputs_str$signif_data_wide,
       main = "P9 Striatum",
       filename = glue("{figout}/gsea_heatmap_hgg_striatum.pdf"))

hm_fun(mat = heatmap_inputs_str$heatmap_data_wide,
       annotation_col = heatmap_inputs_str$col_anno,
       display_number = heatmap_inputs_str$signif_data_wide,
       main = "P9 Striatum",
       filename = glue("{figout}/gsea_heatmap_hgg_striatum.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg_striatum.png"))
```

<img src="/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/03///gsea_heatmap_hgg_striatum.png" width="2717" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/03//gsea_heatmap_hgg_striatum...*]~</span>

### Adult V-SVZ

Finally, generate one for the adult V-SVZ signatures only:


```r
heatmap_inputs_svz <- prep_gsea_heatmap(fgsea_df,
                                        signatures = signif_in_idh,
                                        filters = quos(Sample == "V-SVZ"),
                                        row_order = row_order)
# Generate heatmap
hm_fun(mat = heatmap_inputs_svz$heatmap_data_wide,
       annotation_col = heatmap_inputs_svz$col_anno,
       display_number = heatmap_inputs_svz$signif_data_wide,
       main = "Adult V-SVZ",
       filename = glue("{figout}/gsea_heatmap_hgg_svz.pdf"))

hm_fun(mat = heatmap_inputs_svz$heatmap_data_wide,
       annotation_col = heatmap_inputs_svz$col_anno,
       display_number = heatmap_inputs_svz$signif_data_wide,
       main = "Adult V-SVZ",
       filename = glue("{figout}/gsea_heatmap_hgg_svz.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg_svz.png"))
```

<img src="/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/03///gsea_heatmap_hgg_svz.png" width="2067" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/03//gsea_heatmap_hgg_svz...*]~</span>

<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-07-13 17:27:33
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
##  [1] bindrcpp_0.2.2  cowplot_0.9.4   ggrepel_0.8.0   pheatmap_1.0.12
##  [5] cytobox_0.6.1   scales_1.0.0    fgsea_1.8.0     Rcpp_1.0.0     
##  [9] glue_1.4.1      magrittr_1.5    forcats_0.3.0   stringr_1.3.1  
## [13] dplyr_0.7.7     purrr_0.3.4     readr_1.3.1     tidyr_0.8.2    
## [17] tibble_3.0.1    ggplot2_3.1.0   tidyverse_1.2.1 here_0.1       
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.2.0        snow_0.4-3          backports_1.1.3    
##   [4] Hmisc_4.2-0         fastmatch_1.1-0     plyr_1.8.4         
##   [7] igraph_1.2.2        lazyeval_0.2.1      splines_3.5.0      
##  [10] BiocParallel_1.16.5 digest_0.6.16       foreach_1.4.4      
##  [13] htmltools_0.3.6     viridis_0.5.1       lars_1.2           
##  [16] gdata_2.18.0        checkmate_1.9.1     cluster_2.0.7-1    
##  [19] mixtools_1.1.0      ROCR_1.0-7          modelr_0.1.3       
##  [22] R.utils_2.7.0       colorspace_1.4-0    rvest_0.3.2        
##  [25] haven_2.0.0         xfun_0.12           crayon_1.3.4       
##  [28] jsonlite_1.6        bindr_0.1.1         survival_2.41-3    
##  [31] zoo_1.8-4           iterators_1.0.10    ape_5.2            
##  [34] gtable_0.2.0        kernlab_0.9-27      prabclus_2.2-7     
##  [37] DEoptimR_1.0-8      mvtnorm_1.0-10      bibtex_0.4.2       
##  [40] metap_1.1           dtw_1.20-1          viridisLite_0.3.0  
##  [43] htmlTable_1.13.1    reticulate_1.10     foreign_0.8-70     
##  [46] bit_1.1-14          proxy_0.4-22        mclust_5.4.2       
##  [49] SDMTools_1.1-221    Formula_1.2-3       tsne_0.1-3         
##  [52] stats4_3.5.0        htmlwidgets_1.3     httr_1.4.0         
##  [55] gplots_3.0.1.1      RColorBrewer_1.1-2  fpc_2.1-11.1       
##  [58] acepack_1.4.1       modeltools_0.2-22   ellipsis_0.2.0.1   
##  [61] Seurat_2.3.4        ica_1.0-2           pkgconfig_2.0.2    
##  [64] R.methodsS3_1.7.1   flexmix_2.3-14      nnet_7.3-12        
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
