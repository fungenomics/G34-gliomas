---
title: "02 - GSEA with cell type gene sigantures"
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
## G34-gliomas-repo/bulk/output/02
```

```
## G34-gliomas-repo/bulk/figures/02//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

In this document, we gather a reference of developmental gene signatures from several brain cell types and studies, and evaluate their enrichment in the differentially-expressed genes between G34 tumors and other groups. The enrichment analysis is performed with GSEA, using the package fgsea.

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

## Cluster annotations

First we load cell type annotations for each of the datasets, since they are often
labelled by abbreviations. Create uniform dataframes, with at least three columns:
Cluster, Cell_type, and Age.


```r
jessa_table_2a <- read_tsv(here("bulk/data/Jessa2019_Table_2a.tsv")) %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  mutate(Age = paste0("Mouse ", Age))
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
nowakowski_cluster_anno <- read_csv(here("bulk/data/Nowakowski2017_Table_4_cluster_labels_with_colours.csv")) %>%
  select(Cluster = `Cluster Name`,
         Cell_type = `Cluster Interpretation`) %>% 
  mutate(Cluster = paste0("HF ",  Cluster)) %>% 
  mutate(Age = "Human fetal")
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
# Give these clusters distinct numbers since they have the same "cell type"
nowakowski_cluster_anno[nowakowski_cluster_anno$Cell_type == "MGE newborn neurons", ]$Cell_type <- paste0(nowakowski_cluster_anno[nowakowski_cluster_anno$Cell_type == "MGE newborn neurons", ]$Cell_type, 1:5)

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
cell_type_anno <- bind_rows(jessa_table_2a %>% select(Sample, Age, Cell_type, Cluster),
                            nowakowski_cluster_anno,
                            velmeshev_anno)
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
                  "Oligodendrocyte precursors" = "#e0de53",
                  "Oligodendrocytes" = "#b4e04e",
                  "Astrocytes" = "#00a385",
                  "Neuronal progenitors" = "#ffbda3",
                  "Prenatal inhib. neurons" = "#8798C8",
                  "Pediatric/adult inhib. neurons" = "#135ca0",
                  "Other neurons" = "darkred")

palette_age <- c("Human fetal" = "#458cde",
                 "Mouse E12.5" = "#73a8e6",
                 "Mouse E15.5" = "#a4c6ed",
                 "Mouse P0" = "#f2b39d",
                 "Mouse P3" = "tomato",
                 "Mouse P6" = "red",
                 "Human ped/adult" = "darkred")
```

## Fetal mouse brain

Load signatures from our Nature Genetics paper:


```r
# Get data from paper, filter to clusters with a valid signature, from mouse samples
jessa_table_2a_mouse_noblacklist <- jessa_table_2a %>% 
  filter(Species == "Mouse" & !is.na(Signature))

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

```
## ...writing description of atlas_signatures.Rds to G34-gliomas-repo/bulk/output/02/atlas_signatures.desc
```

## Fetal human brain

We've previously extracted gene signatures from this study (provided in Table S5) and prepared them in a
similar format:


```r
load(here("bulk/data/Nowakowski2017_signatures.Rda"))

# Denote as human fetal
names(nowakowski_signatures$hg_ens) <- paste0("HF ", names(nowakowski_signatures$hg_ens))
names(nowakowski_signatures$hg_sym) <- paste0("HF ", names(nowakowski_signatures$hg_sym))
```


## Pediatric/adult human brain


```r
velmeshev_signatures <- readRDS(here("bulk/data/Velmeshev2019_signatures.Rda"))

# Denote as human postnatal/pediatric
names(velmeshev_signatures$hg_ens) <- paste0("HP ", names(velmeshev_signatures$hg_ens))
names(velmeshev_signatures$hg_sym) <- paste0("HP ", names(velmeshev_signatures$hg_sym))
```


# Analysis

With the cell type signatures derived from the various scRNAseq studies, we can now
use these as input to GSEA, and evaluate their enrichment in differentially expressed
genes between G34R/V tumors vs. other pediatric brain tumor types.

## Prep stats


```r
pipeline_path <- "../../../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/"

# Keep all the DGE analyses between tumor groups, accounting for batch
comps <- list.files(pipeline_path, full.names = TRUE)[
  grepl("^HGG", list.files(pipeline_path)) &
    grepl("vs", list.files(pipeline_path)) &
    grepl("batch_covariate", list.files(pipeline_path))]

names(comps) <- gsub("_batch_covariate", "", basename(comps))

# Get the file with the DESeq2 DGE output, and create a vector with the Ensembl
# gene IDs and stat for that gene, in each comparison
stats <- map(comps, ~ read_tsv(Sys.glob(file.path(.x, "diff/Ensembl.ensGene.exon/*.tsv"))) %>%
               filter(!is.na(stat)) %>%
               separate(ID, into = c("Ens", "Sym"), sep = ":") %>%
               select(Ens, stat) %>%
               deframe())
```

```
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
```

```r
rr_saveRDS(file = glue("{out}/stats.Rds"),
           desc = "List of DESeq2 stats from differential expression analyses between tumor groups, to use for GSEA",
           stats)
```

```
## ...writing description of stats.Rds to G34-gliomas-repo/bulk/output/02/stats.desc
```

## Run


```r
bulk_fgsea <- map(stats, ~ fgsea(pathways = c(atlas_signatures$hg_ens,
                                              nowakowski_signatures$hg_ens,
                                              velmeshev_signatures$hg_ens),
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
## ...writing description of bulk_fgsea.Rds to G34-gliomas-repo/bulk/output/02/bulk_fgsea.desc
```

## Filter gene signatures

Remove clusters which are proliferating based on either labeled as such, or G2M score,
as cell-cycle genes may confound GSEA enrichment scores particularly in these tumor comparisons


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

![](/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/02//filter_signatures-1.png)<!-- -->

```r
atlas_sigs_keep <- infosig_with_cc %>%
  # Remove anything with G2M scores above the threshold, OR which is labelled
  # as proliferating, as indicated by a trailing "-P" in the cluster abbreviation
  filter(Median_G2M < 0.35 & !grepl("-P$", Cluster)) %>%
  # Remove anything labeled Other/Other neurons
  filter(!(Cell_class %in% c("Other", "Other neurons"))) %>%
  # Filter to only forebrain signatures
  filter(Structure == "Forebrain") %>%
  pull(Cluster)

# 2. Nowakowski
nowakowski_sigs_keep <- nowakowski_cluster_anno %>%
  # Filter out proliferating clusters
  filter(!grepl("-P$|div|G2M/S|MGE-IPC1|MGE-IPC2|prolif", Cluster)) %>%
  # Filter out Unknown
  filter(!grepl("Unknown", Cell_type)) %>%
  pull(Cluster)

# 3. Velmeshev - keep all signatures
velmeshev_sigs_keep <- names(velmeshev_signatures$hg_sym)

bulk_fgsea_filt <- map(bulk_fgsea, ~ .x %>%
                         filter(pathway %in% c(atlas_sigs_keep, nowakowski_sigs_keep, velmeshev_sigs_keep)))
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
             path = glue("{out}/fgsea_df.Rda"),
             desc = "Tidied dataframe with GSEA results for bulk tumor group comparisons")
```

```
## ...writing description of fgsea_df.Rda to G34-gliomas-repo/bulk/output/02/fgsea_df.Rda
```

## Top leading edge

For the IDH comparison, let's filter down the leading edge genes to those which
are the top-ranked ones i.e. rank < 2000, to ensure we have good lists to work with.


```r
# Generate stats based on the gene symbols
stats_idh_sym <- read_tsv(Sys.glob(file.path(comps["HGG-G34R.V_vs_HGG-IDH"],
                                             "diff/Ensembl.ensGene.exon/*.tsv"))) %>%
  filter(!is.na(stat)) %>%
  separate(ID, into = c("Ens", "Sym"), sep = ":") %>%
  select(Sym, stat) %>%
  deframe()
```

```
## Parsed with column specification:
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
```

```r
# Apply filters and order by score
fgsea_df_idh_filt <- fgsea_df %>%
  filter(Comparison == "HGG-G34R.V_vs_HGG-IDH") %>%
  filter(padj < 0.01)

# Gather the pathways
pathways <- c(atlas_signatures$hg_sym,
              nowakowski_signatures$hg_sym,
              velmeshev_signatures$hg_sym)

fgsea_df_idh_filt$leadingEdge_top <- NA

# Get the top-ranked leading edge genes
for (i in 1:nrow(fgsea_df_idh_filt)) {

  genes <- str_split(fgsea_df_idh_filt[i, ]$leadingEdge, ", ") %>% getElement(1)
  sig <- fgsea_df_idh_filt[i, ]$Signature
  nes <- fgsea_df_idh_filt[i, ]$NES

  fgsea_df_idh_filt[i, ]$leadingEdge_top <- get_gene_ranks(pathway   = pathways[[sig]],
                                                           stats     = stats_idh_sym,
                                                           direction = ifelse(nes > 0, "enriched", "depleted")) %>%
    filter(gene %in% genes & rank < 2000) %>%
    pull(gene) %>%
    as.character() %>%
    glue_collapse(", ")

}
```

<!-- # Visualization -->

<!-- ```{r palettes_tumors} -->

<!-- palette_groups <- c("HGG-G34R/V"             = "#41a334", -->
<!--                     "HGG-H3.3-K27M-pons"     = "#a30300", -->
<!--                     "HGG-H3.3-K27M-thalamic" = "#f72d2a", -->
<!--                     "HGG-IDH"                = "#da00f1", -->
<!--                     "HGG-WT"                 = "#ada7a6", -->
<!--                     "PNET-HGNET-BCOR"        = "#1b24d1") -->

<!-- ``` -->

<!-- ## GSEA enrichment plots -->

<!-- The most enriched or depleted signature for each cell type: -->

<!-- ```{r enrich_plots_simple, fig.width = 10, fig.height = 6, dependson = 'prep_stats', echo_label = TRUE} -->

<!-- # A function to add the NES and adj p-value to the enrichment plot -->
<!-- plot_enrich_with_stats <- function(sig, title, source, colour) { -->

<!--   signatures = switch(source, -->
<!--                       "nowakowski" = nowakowski_signatures$hg_ens, -->
<!--                       "atlas" = atlas_signatures$hg_ens, -->
<!--                       "velmeshev" = velmeshev_signatures$hg_ens) -->

<!--   plotEnrichment2(signatures[[sig]], -->
<!--                   stats$`HGG-G34R.V_vs_HGG-IDH`, -->
<!--                   colour = colour) + -->
<!--     labs(title = paste0(title, "\nNES: ", -->
<!--                         fgsea_df_idh_filt %>% filter(Signature == sig) %>% pull(NES) %>% round(2), -->
<!--                         ", p_adj: ", -->
<!--                         fgsea_df_idh_filt %>% filter(Signature == sig) %>% pull(padj) %>% scales::scientific())) -->
<!-- } -->

<!-- # Plot for select signatures -->
<!-- p1 <- plot_enrich_with_stats("HF RG-early", "Fetal radial glia", "nowakowski", "red") -->

<!-- p2 <- plot_enrich_with_stats("HF nIN1", "Newborn inhibitory neurons", "nowakowski", "red") -->

<!-- p3 <- plot_enrich_with_stats("HP IN-PV", "Mature inhibitory neurons", "velmeshev", "blue") -->

<!-- p4 <- plot_enrich_with_stats("F-p6 OPC", "Oligodendrocyte precursors", "atlas", "blue") -->

<!-- p5 <- plot_enrich_with_stats("HF EN-PFC1", "Fetal excitatory neurons", "nowakowski", "blue") -->

<!-- # Plot the enrichment of astrocytes against HGNET-BCOR -->
<!-- p6 <- plotEnrichment2(atlas_signatures$hg_ens[["F-p0 ASTR2"]], -->
<!--                       stats$`HGG-G34R.V_vs_PNET-HGNET-BCOR`, -->
<!--                       colour = "red") + -->
<!--   labs(title = paste0("Astrocytes", "\nNES: ", -->
<!--                       fgsea_df %>% filter(Comparison == "HGG-G34R.V_vs_PNET-HGNET-BCOR" & Signature == sig) %>% -->
<!--                         pull(NES) %>% round(2), -->
<!--                       ", p_adj: ", -->
<!--                       fgsea_df %>% filter(Comparison == "HGG-G34R.V_vs_PNET-HGNET-BCOR" & Signature == sig) %>% -->
<!--                         pull(padj) %>% scales::scientific())) -->

<!-- plot_grid(p1, p2, p6, p3, p4, p5, ncol = 3) -->

<!-- ``` -->

## Heatmap across comparisons

Define a function to summarize the cell types into broader "classes",
distinguishing between different types of inhibitory neurons.


```r
summarizeCellTypes3 <- function(df, cluster_col) {

  cc_quo <- enquo(cluster_col)

  df %>%
    mutate(
      # Define some broader cell type classes
      Type = case_when(
        grepl("RG|[Rr]adial", !!cc_quo) & !grepl("NRGN", !!cc_quo) & !grepl("-P.{0,1}$|prolif", !!cc_quo) ~ "RGC",
        grepl("EXIP|INIP|NEURP|IP|[Ii]ntermediate", !!cc_quo) & !grepl("VIP", !!cc_quo) ~ "Neuronal progenitors",
        grepl("MGE|CGE|SST|PV|[Ii]nhib", !!cc_quo) ~ "Prenatal inhib. neurons",
        grepl("PV|SST|VIP|SV2C|Somato|Parv", !!cc_quo) ~ "Pediatric/adult inhib. neurons",
        grepl("CEX|PEX|[Ee]xcit|NRGN|[Nn]eu", !!cc_quo) ~ "Other neurons",
        grepl("OPC|Oligodendrocyte progenitor cell|Oligodendrocyte precursor", !!cc_quo) ~ "Oligodendrocyte precursors",
        grepl("NFOL|MOL|[Oo]ligo", !!cc_quo) ~ "Oligodendrocytes",
        grepl("ASTR|Astr", !!cc_quo) ~ "Astrocytes",
        TRUE ~ "Other")) %>%
    mutate(Type = factor(Type, levels = c("RGC",
                                          "RGC (prolif.)",
                                          "Neuronal progenitors",
                                          "Prenatal inhib. neurons",
                                          "Pediatric/adult inhib. neurons",
                                          "Other neurons",
                                          "Glial progenitors",
                                          "Oligodendrocyte precursors",
                                          "Oligodendrocytes",
                                          "Astrocytes",
                                          "Choroid/ependymal",
                                          "Immune",
                                          "Non-neuroectoderm",
                                          "Other")))

}
```

Generate the input data for the heatmap:


```r
# Select only signatures which are signif enriched/depleted in the IDH comparison
signif_in_idh <- fgsea_df_idh_filt$Signature

heatmap_data_long <- fgsea_df %>%
  filter(Signature %in% signif_in_idh) %>%
  # Subset to relevant cell types
  summarizeCellTypes3(Cell_type) %>%
  filter(Type %in% c("RGC",
                     "Prenatal inhib. neurons",
                     "Pediatric/adult inhib. neurons",
                     "Neuronal progenitors",
                     "Other neurons",
                     "Oligodendrocytes",
                     "Oligodendrocyte precursors",
                     "Astrocytes")) %>%
  mutate(Cell_type = paste0(Age, " ", Cell_type))

col_order <- heatmap_data_long %>%
  distinct(Cell_type, Type) %>%
  arrange(Type) %>%
  pull(Cell_type)

col_anno <- heatmap_data_long %>%
  distinct(Cell_type, Age, Type) %>%
  mutate(Species = ifelse(grepl("Mouse", Cell_type), "Mouse", "Human")) %>%
  data.frame() %>%
  tibble::column_to_rownames(var = "Cell_type")

# Include only the high-grade gliomas at first
row_order <- c("HGG-G34R.V_vs_HGG-IDH",
               "HGG-G34R.V_vs_HGG-H3.3-K27M-pons",
               "HGG-G34R.V_vs_HGG-H3.3-K27M-thalamic",
               "HGG-G34R.V_vs_HGG-WT")

heatmap_data_wide <- heatmap_data_long %>%
  select(Comparison, Cell_type, NES) %>%
  spread(Cell_type, NES) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "Comparison") %>%
  .[row_order, col_order]
```


```r
hm_fun <- partial(pheatmap,
                  border_color = NA,
                  color = colorRampPalette(c("blue","white","red"))(100),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  cellwidth = 12,
                  cellheight = 12,
                  annotation_colors = list("Age"  = palette_age,
                                           "Type" = palette_type,
                                           "Species" = c("Human" = "black", "Mouse" = "gray70")))

# Generate heatmap
hm_fun(mat = heatmap_data_wide,
       annotation_col = col_anno,
       filename = glue("{figout}/gsea_heatmap_hgg.pdf"))

hm_fun(mat = heatmap_data_wide,
       annotation_col = col_anno,
       filename = glue("{figout}/gsea_heatmap_hgg.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg.png"))
```

<img src="/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/02///gsea_heatmap_hgg.png" width="4000" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/02//gsea_heatmap_hgg...*]~</span>

Generate a second for HGNET-BCOR:


```r
signif_in_bcor <- fgsea_df %>%
  filter(Comparison == "HGG-G34R.V_vs_PNET-HGNET-BCOR" & padj < 0.01) %>%
  pull(Signature)

heatmap_data_long <- fgsea_df %>%
  filter(Signature %in% signif_in_bcor) %>%
  # Subset to relevant cell types
  summarizeCellTypes3(Cell_type) %>%
  filter(Type %in% c("RGC",
                     "Prenatal inhib. neurons",
                     "Pediatric/adult inhib. neurons",
                     "Neuronal progenitors",
                     "Other neurons",
                     "Oligodendrocytes",
                     "Oligodendrocyte precursors",
                     "Astrocytes")) %>%
  mutate(Cell_type = paste0(Age, " ", Cell_type)) %>%
  filter(Comparison == "HGG-G34R.V_vs_PNET-HGNET-BCOR")

col_order <- heatmap_data_long %>%
  distinct(Cell_type, Type) %>%
  arrange(Type) %>%
  pull(Cell_type)

col_anno <- heatmap_data_long %>%
  distinct(Cell_type, Age, Type) %>%
  data.frame() %>%
  tibble::column_to_rownames(var = "Cell_type")

heatmap_data_wide <- heatmap_data_long %>%
  distinct(Cell_type, Type, .keep_all = TRUE)  %>%
  select(Cell_type, NES) %>%
  spread(Cell_type, NES) %>%
  as.data.frame() %>%
  .[ , col_order]

hm_fun(mat = heatmap_data_wide,
       annotation_col = col_anno,
       filename = glue("{figout}/gsea_heatmap_bcor.pdf"))

hm_fun(mat = heatmap_data_wide,
       annotation_col = col_anno,
       filename = glue("{figout}/gsea_heatmap_bcor.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_bcor.png"))
```

<img src="/mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/HGG-G34/G34-gliomas-repo/bulk/figures/02///gsea_heatmap_bcor.png" width="2613" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas-repo/bulk/figures/02//gsea_heatmap_bcor...*]~</span>



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-07-09 08:18:16
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
##  [67] labeling_0.3        reshape2_1.4.3      tidyselect_1.1.0   
##  [70] rlang_0.4.6         munsell_0.5.0       cellranger_1.1.0   
##  [73] tools_3.5.0         cli_1.0.1           generics_0.0.2     
##  [76] broom_0.5.1         ggridges_0.5.1      evaluate_0.12      
##  [79] yaml_2.2.0          npsurv_0.4-0        knitr_1.21         
##  [82] bit64_0.9-7         fitdistrplus_1.0-14 robustbase_0.93-2  
##  [85] caTools_1.17.1.1    randomForest_4.6-14 RANN_2.6           
##  [88] pbapply_1.4-0       nlme_3.1-137        R.oo_1.22.0        
##  [91] xml2_1.2.0          hdf5r_1.0.0         compiler_3.5.0     
##  [94] rstudioapi_0.9.0    png_0.1-7           lsei_1.2-0         
##  [97] stringi_1.2.4       lattice_0.20-35     trimcluster_0.1-2.1
## [100] Matrix_1.2-14       vctrs_0.3.1         pillar_1.4.4       
## [103] lifecycle_0.2.0     Rdpack_0.10-1       lmtest_0.9-36      
## [106] data.table_1.12.0   bitops_1.0-6        irlba_2.3.3        
## [109] gbRd_0.4-11         R6_2.3.0            latticeExtra_0.6-28
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
