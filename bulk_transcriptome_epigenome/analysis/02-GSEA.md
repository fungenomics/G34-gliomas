---
title: "02 - GSEA with cell type gene signatures"
date: "11 December, 2020"
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
## Cache: ~/tmp/G34-gliomas/bulk_transcriptome_epigenome/02/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/bulk_transcriptome_epigenome/output/02
```

```
## G34-gliomas/bulk_transcriptome_epigenome/figures/02//
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

* [Jessa et al, 2019](https://doi.org/10.1038/s41588-019-0531-7)
* [Nowakowski et al, 2017](https://doi.org/10.1126/science.aap8809)
* [Velmeshev et al, 2019](https://doi.org/10.1126/science.aav8130)

We also obtained two reference datasets capturing the adult SVZ:

* [Anderson et al, 2020](https://doi.org/10.1016/j.celrep.2020.02.030)
* [Mizrak et al, 2019](https://doi.org/10.1016/j.celrep.2018.12.044)


# Libraries


```r
# General
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(glue)
library(tibble)

# For analysis
library(fgsea)
library(cytobox)

# For plotting
library(ggplot2)
library(scales)
library(pheatmap)
library(ggrepel)
library(cowplot)
ggplot2::theme_set(theme_min(base_size = 13))

# Custom
source(here::here(subdir, "analysis/functions.R"))
```

# Load data

## Tumor metadata


```r
bulk_samples <- readxl::read_xlsx(
  here::here(subdir, "data/2020-02-25_bulkRNAseq_metadata.xlsx"))
```

The following sections loads the reference datasets, including cell-type specific gene
signatures, corresponding palettes, and cluster annotations.

## Fetal mouse brain

Load signatures from our Nature Genetics paper:


```r
# Annotation
jessa_table_2a <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_Table_2a.tsv")) %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  mutate(Age = paste0("Mouse ", Age),
         Dataset = "Jessa et al 2019")

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

# All signatures are 100-genes
all(map_dbl(atlas_signatures$mm_sym, length) == 100)
```

```
## [1] TRUE
```

```r
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
nowakowski_cluster_anno <- read_csv(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski2017_tableS4_cluster_labels_with_colours.csv")) %>%
  select(Cluster = `Cluster Name`,
         Cell_type = `Cluster Interpretation`) %>% 
  mutate(Cluster = paste0("HF ",  Cluster)) %>% 
  mutate(Age = "Human fetal",
         Sample = "Human fetal cortex/MGE",
         Dataset = "Nowakowski et al 2017")
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
load(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski_signatures.Rda"))

# Denote as human fetal
names(nowakowski_signatures$hg_ens) <- paste0("HF ", names(nowakowski_signatures$hg_ens))
names(nowakowski_signatures$hg_sym) <- paste0("HF ", names(nowakowski_signatures$hg_sym))
```


## Pediatric/adult human brain


```r
# Annotation
velmeshev_anno <- read_tsv(here("reference_datasets/2019_Velmeshev/processed_data/01-cluster_names_with_interpretation.tsv")) %>% 
  mutate(Age = "Human ped/adult",
         Cluster = paste0("HP ", Cluster),
         Dataset = "Velmeshev et al 2019",
         Sample  = "Human ped/adult cortex")
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
velmeshev_signatures <- readRDS(here("reference_datasets/2019_Velmeshev/processed_data/01-cluster_gene_signatures.Rds"))

# Denote as human postnatal/pediatric
names(velmeshev_signatures$hg_ens) <- paste0("HP ", names(velmeshev_signatures$hg_ens))
names(velmeshev_signatures$hg_sym) <- paste0("HP ", names(velmeshev_signatures$hg_sym))
```


## Adult V-SVZ


```r
# Signatures
mizrak_rep1_signatures <- readRDS(here("reference_datasets/2019_Mizrak/processed_data/03-cluster_gene_signatures.Rds"))
mizrak_rep2_signatures <- readRDS(here("reference_datasets/2019_Mizrak/processed_data/04-cluster_gene_signatures.Rds"))

# Add a prefix to the cluster names
mizrak_rep1_signatures <- map(mizrak_rep1_signatures, ~ .x %>% set_names(paste0("VSVZ-1 ", names(.x))))
mizrak_rep2_signatures <- map(mizrak_rep2_signatures, ~ .x %>% set_names(paste0("VSVZ-2 ", names(.x))))

mizrak_anno <- data.frame(Cluster   = c(names(mizrak_rep1_signatures$hg_sym),
                                        names(mizrak_rep2_signatures$hg_sym)),
                          Cell_type = c(names(mizrak_rep1_signatures$hg_sym),
                                        names(mizrak_rep2_signatures$hg_sym)),
                          Age = "Mouse adult",
                          Sample = "V-SVZ",
                          Dataset = "Mizrak et al 2019")
```


## Striatal SVZ


```r
# Signatures
anderson_signatures <- readRDS(here("reference_datasets/2020_Anderson/processed_data/02-cluster_gene_signatures.Rds"))
anderson_signatures <- map(anderson_signatures, ~ .x %>% set_names(paste0("Str/SVZ ", names(.x))))

# Annotation
anderson_anno <- data.frame(Cluster = names(anderson_signatures$hg_sym),
                            Cell_type = names(anderson_signatures$hg_sym),
                            Age = "Mouse P9",
                            Sample = "Striatum",
                            Dataset = "Anderson et al 2020")
```

## Gather all signatures


```r
# Human gene symbols
signatures_sym <- c(atlas_signatures$hg_sym,
                    nowakowski_signatures$hg_sym,
                    velmeshev_signatures$hg_sym,
                    anderson_signatures$hg_sym,
                    mizrak_rep1_signatures$hg_sym,
                    mizrak_rep2_signatures$hg_sym)

length(signatures_sym)
```

```
## [1] 202
```

```r
rr_saveRDS(file = glue("{out}/signatures_sym.Rds"),
           desc = "Cell type specific gene signatures from several forebrain references (list of vectors, using gene symbols)",
           signatures_sym)
```

```
## ...writing description of signatures_sym.Rds to G34-gliomas/bulk_transcriptome_epigenome/output/02/signatures_sym.desc
```

```r
# Human Ensembl IDs
signatures_ens <- c(atlas_signatures$hg_ens,
                    nowakowski_signatures$hg_ens,
                    velmeshev_signatures$hg_ens,
                    anderson_signatures$hg_ens,
                    mizrak_rep1_signatures$hg_ens,
                    mizrak_rep2_signatures$hg_ens)

length(signatures_sym)
```

```
## [1] 202
```

```r
rr_saveRDS(file = glue("{out}/signatures_ens.Rds"),
           desc = "Cell type specific gene signatures from several forebrain references (list of vectors, using Ensemble gene IDs)",
           signatures_ens)
```

```
## ...writing description of signatures_ens.Rds to G34-gliomas/bulk_transcriptome_epigenome/output/02/signatures_ens.desc
```


## Assemble cluster annotations

Put together the cell type annotations for each of the datasets, since they are often
labelled by abbreviations. Create uniform dataframes, with at least four columns:
Sample, Cluster, Cell_type, and Age. The Cluster column matches the names of
the signatures in the `signatures_*` objects created above.


```r
cell_type_anno <- bind_rows(jessa_table_2a %>% select(Sample, Age, Cell_type, Cluster, Dataset),
                            nowakowski_cluster_anno,
                            velmeshev_anno,
                            mizrak_anno,
                            anderson_anno)

rr_write_tsv(df = cell_type_anno, path = glue("{out}/signature_annotations.tsv"),
             desc = "Cluster info for normal signatures including dataset, age, cell class, etc.")
```

```
## ...writing description of signature_annotations.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/02/signature_annotations.desc
```


## Palettes

Create colour palettes for each dataset, for visualization:


```r
palette_atlas <- jessa_table_2a %>% select(Cluster, Colour) %>% deframe()

palette_nowakowski <- read_csv(here("reference_datasets/2017_Nowakowski/processed_data/nowakowski2017_tableS4_cluster_labels_with_colours.csv")) %>%
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
pipeline_path <- "../../../level3/2020-01_G34_submission1_add_samples/"

# Keep all the DGE analyses between tumor groups, accounting for batch
comps <- list.files(pipeline_path, full.names = TRUE)[
  grepl("^HGG", list.files(pipeline_path)) &
    grepl("vs", list.files(pipeline_path)) &
    grepl("batch_covariate", list.files(pipeline_path)) &
    !grepl("IDH_and_WT", list.files(pipeline_path))]

names(comps) <- gsub("_batch_covariate", "", basename(comps))

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
# Run GSEA for each pathway, for each comparison
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
## ...writing description of bulk_fgsea.Rds to G34-gliomas/bulk_transcriptome_epigenome/output/02/bulk_fgsea.desc
```

## Filter gene signatures

For the forebrain reference, remove clusters which are proliferating based on either labeled as such, or Top2a expression,
as cell-cycle genes may confound GSEA enrichment scores particularly in these tumor comparisons. No clusters from the striatum/SVZ reference were removed.


```r
# 1. Atlas
# Load signatures with cell cycle scores, to quantitatively set a filter
infosig_with_cc <- read_tsv(here("reference_datasets/2019_Jessa/Jessa2019_cluster_cell_cycle_scores.tsv")) %>%
  arrange(desc(Median_G2M)) %>% 
  full_join(jessa_table_2a)
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

```
## Joining, by = c("Sample", "Age", "Species", "Structure", "Cluster_number", "Cell_type", "Cluster", "Cell_class", "Alias", "N_cells", "Colour", "nGene_min", "nUMI_min", "percent_mito_min", "nGene_mean", "nUMI_mean", "percent_mito_mean", "nGene_max", "nUMI_max", "percent_mito_max", "nGene_sd", "nUMI_sd", "percent_mito_sd")
```

```r
infosig_with_cc %>%
  filter(Species == "Mouse") %>%
  mutate(Cluster = factor(Cluster, levels = unique(.$Cluster))) %>%
  ggplot(aes(x = Cluster, y = Median_G2M)) +
  geom_point(aes(colour = Cluster)) +
  scale_colour_manual(values = palette_atlas) +
  # Set an arbitrary threshold, where the G2M score starts tapering off
  geom_hline(yintercept = 0.35) +
  theme_min() +
  rotateX() +
  noLegend()
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02//filter_signatures-1.png)<!-- -->

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
           # Filter out Unknown cluster
           grepl("Unknown", Cell_type)) %>%
  pull(Cluster)

bulk_fgsea_filt <- map(bulk_fgsea, ~ .x %>%
                         filter(!(pathway %in% c(atlas_sigs_drop, nowakowski_sigs_drop))))
```

## Save final gene signature set

Having filtered out certain signatures, we'll save the final set used for analysis:


```r
signatures_ens_filt <- signatures_ens[!(names(signatures_ens) %in% c(atlas_sigs_drop, nowakowski_sigs_drop))]

length(signatures_ens_filt)
```

```
## [1] 164
```

```r
rr_saveRDS(file = glue("{out}/signatures_ens_filt.Rds"),
           desc = "Cell type specific gene signatures from several forebrain references (list of vectors, using ENSEMBL IDs), after filtering due to cell cycle/proliferation/irrelevant or technical signal",
           signatures_ens_filt)
```

```
## ...writing description of signatures_ens_filt.Rds to G34-gliomas/bulk_transcriptome_epigenome/output/02/signatures_ens_filt.desc
```

```r
signatures_sym_filt <- signatures_sym[!(names(signatures_sym) %in% c(atlas_sigs_drop, nowakowski_sigs_drop))]

length(signatures_sym_filt)
```

```
## [1] 164
```

```r
rr_saveRDS(file = glue("{out}/signatures_sym_filt.Rds"),
           desc = "Cell type specific gene signatures from several forebrain references (list of vectors, using symEMBL IDs), after filtering due to cell cycle/proliferation/irrelevant or technical signal",
           signatures_sym_filt)
```

```
## ...writing description of signatures_sym_filt.Rds to G34-gliomas/bulk_transcriptome_epigenome/output/02/signatures_sym_filt.desc
```

```r
# make a table
signatures_ens_filt_collapsed <- map(signatures_ens_filt, ~ glue_collapse(.x, ","))
signatures_sym_filt_collapsed <- map(signatures_sym_filt, ~ glue_collapse(.x, ","))
signatures_filt_df <- data.frame(Cluster       = names(signatures_ens_filt_collapsed),
                                 Genes_ENSEMBL = unlist(unname(signatures_ens_filt_collapsed)),
                                 Genes_symbol  = unlist(unname(signatures_sym_filt_collapsed)))

dim(signatures_filt_df)
```

```
## [1] 164   3
```

```r
rr_write_tsv(df = signatures_filt_df,
             path = glue("{out}/signatures_filt.tsv"),
             desc = "Table with all forebrain gene signatures included in analysis after filtering, both gene symbols and ENSEMBl IDs")
```

```
## ...writing description of signatures_filt.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/02/signatures_filt.desc
```

Make GMX format for GSEA input for sharing:


```r
#' Pad a vector x with NAs to length n
pad <- function(x, n = 120) {
    
    if (length(x) < n) c(x, rep(NA, n - length(x)))
    else x
    
}

# pad all signatures
signatures_sym_filt_padded <- map(signatures_sym_filt, pad)
signatures_sym_filt_padded_df <- data.frame(signatures_sym_filt_padded)
colnames(signatures_sym_filt_padded_df) <- names(signatures_sym_filt_padded)

# for GMX format, the first row should be an optional description
dummy_row <- rep("na", length(signatures_sym_filt_padded))
names(dummy_row) = names(signatures_sym_filt_padded)

signatures_sym_filt_padded_df <- bind_rows(dummy_row, signatures_sym_filt_padded_df)

write_tsv(signatures_sym_filt_padded_df, glue("{out}/signatures_filt.gmx"))

dim(signatures_sym_filt_padded_df)
```

```
## [1] 121 164
```

## Tidy output

Here we transform the fgsea output into tidy format:


```r
# A function to convert the leading edge from ENSG id to symbols, and from
# a character vector a flat character string
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
## ...writing description of fgsea_df.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/02/fgsea_df.desc
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

fgsea_df_idh_filt$leadingEdge_top <- ""

# Get the top-ranked leading edge genes
for (i in 1:nrow(fgsea_df_idh_filt)) {
  
  genes <- stringr::str_split(fgsea_df_idh_filt[i, ]$leadingEdge, ", ") %>% getElement(1)
  sig <- fgsea_df_idh_filt[i, ]$Signature
  nes <- fgsea_df_idh_filt[i, ]$NES
  
  fgsea_df_idh_filt[i, ]$leadingEdge_top <- get_gene_ranks(pathway   = signatures_sym[[sig]],
                                                           stats     = stats_idh_sym,
                                                           direction = ifelse(nes > 0, "enriched", "depleted")) %>%
    filter(gene %in% genes & rank < 2000) %>%
    pull(gene) %>%
    as.character() %>%
    glue_collapse(", ")
  
}

rr_write_tsv(fgsea_df_idh_filt,
             glue("{out}/fgsea_results_padj<0.01_IDH_top_le.tsv"),
             "GSEA results for G34-vs-IDH comparison, filtered to significant results, and with leading edge filtered to top-ranked genes.")
```

```
## ...writing description of fgsea_results_padj<0.01_IDH_top_le.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/02/fgsea_results_padj<0.01_IDH_top_le.desc
```

# Visualization


## Heatmap across comparisons

Helpers for creating the heatmaps:


```r
# Select only signatures which are signif enriched/depleted in the IDH comparison
signif_in_idh <- fgsea_df_idh_filt$Signature

save(signif_in_idh, file = glue("{out}/signif_in_idh.Rda"))

# Helper function to actually generate the heatmap, which will be called several
# times to generate different heatmaps with different inputs/outputs
hm_fun <- partial(pheatmap,
                  border_color = NA,
                  color = colorRampPalette(c("blue","white","red"))(100),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  cellwidth = 12,
                  cellheight = 12,
                  # Prepare the cluster annotation tracks
                  annotation_colors = list("Age"  = palette_age,
                                           "Type" = palette_type,
                                           "Species" = c("Human" = "black", "Mouse" = "gray70")))

# Include only the high-grade gliomas at first
row_order <- c("HGG-G34R.V_vs_HGG-IDH",
               "HGG-G34R.V_vs_HGG-H3.3-K27M-pons",
               "HGG-G34R.V_vs_HGG-H3.3-K27M-thalamic",
               "HGG-G34R.V_vs_HGG-WT")
```

### Forebrain reference


```r
# Prep input
heatmap_inputs_fb <- prep_gsea_heatmap(fgsea_df,
                                       signatures = signif_in_idh,
                                       filters = quos(Dataset %in% c("Jessa et al 2019",
                                                                     "Nowakowski et al 2017",
                                                                     "Velmeshev et al 2019")),
                                       row_order = row_order)

# Save source data manually
write_tsv(heatmap_inputs_fb$heatmap_data_wide, glue("{figout}/gsea_heatmap_hgg.source_data.tsv"))

# Generate heatmap in pdf and png
hm_fun(mat = heatmap_inputs_fb$heatmap_data_wide,
       annotation_col = heatmap_inputs_fb$col_anno,
       main = "Forebrain references",
       filename = glue("{figout}/gsea_heatmap_hgg.pdf"))

hm_fun(mat = heatmap_inputs_fb$heatmap_data_wide,
       annotation_col = heatmap_inputs_fb$col_anno,
       main = "Forebrain references",
       filename = glue("{figout}/gsea_heatmap_hgg.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02///gsea_heatmap_hgg.png" width="4000" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/02//gsea_heatmap_hgg...*]~</span>

Generate a second heatmap for HGNET-BCOR, which is a non-glioma entity used to have
a comparator group for G34 gliomas that are not not glial in origin.


```r
signif_in_bcor <- fgsea_df %>%
  filter(Comparison == "HGG-G34R.V_vs_PNET-HGNET-BCOR" & padj < 0.01) %>%
  pull(Signature)

heatmap_inputs_bcor <- prep_gsea_heatmap(fgsea_df,
                                         signatures = signif_in_bcor,
                                         filters = quos(Dataset %in% c("Jessa et al 2019",
                                                                       "Nowakowski et al 2017",
                                                                       "Velmeshev et al 2019"),
                                                        Comparison == "HGG-G34R.V_vs_PNET-HGNET-BCOR"),
                                         row_order = 1)

# Save source data manually
write_tsv(heatmap_inputs_bcor$heatmap_data_wide, glue("{figout}/gsea_heatmap_bcor.source_data.tsv"))

hm_fun(mat = heatmap_inputs_bcor$heatmap_data_wide,
       annotation_col = heatmap_inputs_bcor$col_anno,
       filename = glue("{figout}/gsea_heatmap_bcor.png"))

hm_fun(mat = heatmap_inputs_bcor$heatmap_data_wide,
       annotation_col = heatmap_inputs_bcor$col_anno,
       filename = glue("{figout}/gsea_heatmap_bcor.pdf"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_bcor.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02///gsea_heatmap_bcor.png" width="3456" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/02//gsea_heatmap_bcor...*]~</span>



### Striatal SVZ

Generate one heatmap for the striatal SVZ signatures only:


```r
heatmap_inputs_str <- prep_gsea_heatmap(fgsea_df,
                                        signatures = signif_in_idh,
                                        filters = quos(Dataset == "Anderson et al 2020"),
                                        row_order = row_order)

# Save source data manually
write_tsv(heatmap_inputs_str$heatmap_data_wide, glue("{figout}/gsea_heatmap_hgg_striatum.source_data.tsv"))

# Generate heatmap
hm_fun(mat = heatmap_inputs_str$heatmap_data_wide,
       annotation_col = heatmap_inputs_str$col_anno,
       main = "P9 Striatum",
       filename = glue("{figout}/gsea_heatmap_hgg_striatum.pdf"))

hm_fun(mat = heatmap_inputs_str$heatmap_data_wide,
       annotation_col = heatmap_inputs_str$col_anno,
       main = "P9 Striatum",
       filename = glue("{figout}/gsea_heatmap_hgg_striatum.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg_striatum.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02///gsea_heatmap_hgg_striatum.png" width="3018" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/02//gsea_heatmap_hgg_striatum...*]~</span>

### Adult V-SVZ

Finally, generate one for the adult V-SVZ signatures only:


```r
heatmap_inputs_svz <- prep_gsea_heatmap(fgsea_df,
                                        signatures = signif_in_idh,
                                        filters = quos(Dataset == "Mizrak et al 2019"),
                                        row_order = row_order)

# Save source data manually
write_tsv(heatmap_inputs_svz$heatmap_data_wide, glue("{figout}/gsea_heatmap_hgg_svz.source_data.tsv"))

# Generate heatmap
hm_fun(mat = heatmap_inputs_svz$heatmap_data_wide,
       annotation_col = heatmap_inputs_svz$col_anno,
       main = "Adult V-SVZ",
       filename = glue("{figout}/gsea_heatmap_hgg_svz.pdf"))

hm_fun(mat = heatmap_inputs_svz$heatmap_data_wide,
       annotation_col = heatmap_inputs_svz$col_anno,
       main = "Adult V-SVZ",
       filename = glue("{figout}/gsea_heatmap_hgg_svz.png"))

knitr::include_graphics(glue("{figout}/gsea_heatmap_hgg_svz.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02///gsea_heatmap_hgg_svz.png" width="2368" /><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/02//gsea_heatmap_hgg_svz...*]~</span>

## GSEA enrichment plots

For further analysis of these results, we need specific signatures to investigate.
In order to do this, we'll select for each cell type, the most enriched or
depleted signature.


```r
# A function to add the NES and adj p-value to the enrichment plot, based on
# the functions included in the fgsea package
plot_enrich_with_stats <- function(sig, title, source, colour, plot_num) {
  
  signatures = switch(source,
                      "nowakowski" = nowakowski_signatures$hg_ens,
                      "atlas"      = atlas_signatures$hg_ens,
                      "velmeshev"  = velmeshev_signatures$hg_ens,
                      "anderson"   = anderson_signatures$hg_ens)
  
  plotEnrichment2(signatures[[sig]],
                  stats$`HGG-G34R.V_vs_HGG-IDH`,
                  colour = colour,
                  plot_num = plot_num) +
    labs(title = paste0(title, "\nNES: ",
                        fgsea_df_idh_filt %>% filter(Signature == sig) %>% pull(NES) %>% round(2),
                        ", p_adj: ",
                        fgsea_df_idh_filt %>% filter(Signature == sig) %>% pull(padj) %>% scales::scientific()))
}

# Plot for most enriched or depleted signatures
p1 <- plot_enrich_with_stats("HF RG-early",
                             "Fetal radial glia",
                             "nowakowski",
                             "red",
                             plot_num = 1)

p2 <- plot_enrich_with_stats("HF nIN1",
                             "Newborn inhibitory neurons",
                             "nowakowski",
                             "red",
                             plot_num = 2)

p3 <- plot_enrich_with_stats("Str/SVZ 5-Neurogenic progenitor",
                             "SVZ Gsx2+ progenitors",
                             "anderson",
                             "red",
                             plot_num = 3)

p4 <- plot_enrich_with_stats("HP IN-PV",
                             "Mature inhibitory neurons",
                             "velmeshev",
                             "blue",
                             plot_num = 4)

p5 <- plot_enrich_with_stats("HP OPC",
                             "Oligodendrocyte precursors",
                             "velmeshev",
                             "blue",
                             plot_num = 5)

p6 <- plot_enrich_with_stats("HF EN-PFC1",
                             "Fetal excitatory neurons",
                             "nowakowski",
                             "blue",
                             plot_num = 6)

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02//gsea_enrichment_plots-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/02//gsea_enrichment_plots...*]~</span>


## Prepare supplementary table

Output results of gene set enrichment analysis (GSEA) for G34R/V compared to other tumor entities:


```r
table_gsea <- fgsea_df %>% 
  filter(padj < 0.01) %>% 
  select(Comparison, Cell_type, Age, Dataset, pval, padj, ES, NES, leadingEdge) %>% 
  mutate(Comparison = factor(Comparison, levels = c(
    "HGG-G34R.V_vs_HGG-IDH",
    "HGG-G34R.V_vs_HGG-WT",
    "HGG-G34R.V_vs_HGG-H3.3-K27M-pons",
    "HGG-G34R.V_vs_HGG-H3.3-K27M-thalamic",
    "HGG-G34R.V_vs_PNET-HGNET-BCOR"
  ))) %>% 
  arrange(Comparison, desc(NES))

rr_write_tsv(table_gsea,
             glue("{out}/table_gsea.tsv"),
             desc = "Output of GSEA analysis of tumor groups using several reference datasets, filtered to significant results")
```

```
## ...writing description of table_gsea.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/02/table_gsea.desc
```

## Confirmation of signal by direct expression of gene programs

In this analysis, we have looked at differential signal between tumor groups,
but let's interrogate directly the expression of genes in the interneuron gene
program.

Set up palettes:


```r
palette_groups <- c(
  "G34R/V"            = "cyan4",
  "K27M-pons"         = "gold1",
  "K27M-thal"         = "goldenrod3",
  "IDH"               = "#883067")
```

Select gene signatures to visualize in the heatmap -- we'll use the same ones as for the GSEA enrichment plots:


```r
signatures <- c("HF nIN1",  # Newborn interneuron signature
                "HP OPC") # Oligodendrocyte signature

leading_edge <- fgsea_df_idh_filt %>% 
  filter(Signature %in% signatures) %>% 
  select(Signature, leadingEdge_top) %>% 
  tibble::deframe() %>% 
  map(~ .x %>% stringr::str_split(", ") %>% unlist())

heatmap_genes <- unique(unlist(leading_edge))
```

Load normalized counts for genes of interest from the pipeline run:


```r
# Load counts from the pipeline and subset to genes of interest
counts_subset_heatmap <- extract_pipeline_counts(path = file.path(pipeline_path,
                                                                  "All_tumor_samples/counts/Ensembl.ensGene.exon.vst.tsv.gz"),
                                                 goi = heatmap_genes,
                                                 long = FALSE)
```

```
## Joining, by = "gene_symbol"
```

```r
# gene x sample heatmap
# Set genes with NA expression to 0
counts_subset_heatmap[is.na(counts_subset_heatmap)] <- 0

# Reformat the ENSEMBL_ID:symbol format to get just gene symbol for simplicity
rownames(counts_subset_heatmap) <- rownames(counts_subset_heatmap) %>% strsplit(":") %>% sapply(`[[`, 2)

# Remove samples from the HGG-WT and HGNET-BCOR subgorups
filt_samples <- bulk_samples %>% filter(Group %in% c("HGG-WT", "PNET-HGNET-BCOR")) %>% pull(Sample)
counts_subset_heatmap_filt <- counts_subset_heatmap[, !colnames(counts_subset_heatmap) %in% filt_samples]

# z-score
z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Apply z-scoring across samples, gene x sample
counts_subset_heatmap_filt_z <- apply(counts_subset_heatmap_filt, 1, z_score)

# Construct annotations for gene programs & samples
gene_anno <- data.frame("gene" = rownames(counts_subset_heatmap)) %>% 
  distinct(gene) %>% 
  mutate("Program" = case_when(
    gene %in% signatures_sym[["HP OPC"]] ~ "Oligodendrocyte",
    gene %in% signatures_sym[["HF nIN1"]] ~ "Interneuron",
    TRUE ~ "Other")) %>% 
  tibble::column_to_rownames(var = "gene")

sample_anno <- bulk_samples %>% 
  select(Sample, Group) %>% 
  mutate(Group = recode(Group,
                        "HGG-G34R/V" = "G34R/V",
                        "HGG-H3.3-K27M-pons" = "K27M-pons",
                        "HGG-H3.3-K27M-thalamic" = "K27M-thal",
                        "HGG-IDH" = "IDH",
                        "HGG-WT" = "WT",
                        "PNET-HGNET-BCOR" = "HGNET-BCOR")) %>% 
  filter(!(Group %in% c("WT", "HGNET-BCOR"))) %>% 
  tibble::column_to_rownames(var = "Sample")
```



```r
# Plot as gene x sample
hm_fun <- purrr::partial(pheatmap,
                         mat = t(counts_subset_heatmap_filt_z),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         scale = "none",
                         annotation_row = gene_anno,
                         annotation_col = sample_anno,
                         annotation_colors = list("Group" = palette_groups,
                                                  "Program" = c("Interneuron" = "#8798C8",
                                                                "Oligodendrocyte" = "#b4e04e")),
                         cellwidth = 8,
                         cellheight = 8,
                         border_color = NA,
                         color = rdbu,
                         show_colnames = FALSE,
                         fontsize_row = 5)

hm_fun(filename = glue("{figout}/heatmap_gene_programs_vst.png"))
hm_fun(filename = glue("{figout}/heatmap_gene_programs_vst.pdf"))

knitr::include_graphics(glue("{figout}/heatmap_gene_programs_vst.png"))
```

<img src="/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/02///heatmap_gene_programs_vst.png" width="2237" />



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-12-11 17:43:26
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [71d3b7e] 2020-09-28: Ignore slurm submission script
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
##  [1] cowplot_0.9.4   ggrepel_0.8.0   pheatmap_1.0.12 scales_1.1.1   
##  [5] ggplot2_3.1.0   cytobox_0.6.1   fgsea_1.8.0     Rcpp_1.0.5     
##  [9] tibble_3.0.3    glue_1.4.2      magrittr_1.5    purrr_0.3.4    
## [13] readr_1.3.1     dplyr_0.8.0     tidyr_0.8.2     here_0.1       
## 
## loaded via a namespace (and not attached):
##   [1] snow_0.4-3          backports_1.1.9     Hmisc_4.2-0        
##   [4] fastmatch_1.1-0     sn_1.6-2            plyr_1.8.6         
##   [7] igraph_1.2.5        lazyeval_0.2.2      splines_3.5.1      
##  [10] BiocParallel_1.16.6 TH.data_1.0-10      digest_0.6.25      
##  [13] foreach_1.5.0       htmltools_0.5.0     viridis_0.5.1      
##  [16] lars_1.2            gdata_2.18.0        checkmate_2.0.0    
##  [19] cluster_2.0.7-1     mixtools_1.2.0      ROCR_1.0-7         
##  [22] R.utils_2.10.1      sandwich_2.5-1      colorspace_1.4-1   
##  [25] xfun_0.17           jsonlite_1.7.1      crayon_1.3.4       
##  [28] survival_2.41-3     zoo_1.8-8           iterators_1.0.12   
##  [31] ape_5.4-1           gtable_0.3.0        kernlab_0.9-29     
##  [34] prabclus_2.3-2      BiocGenerics_0.28.0 DEoptimR_1.0-8     
##  [37] mvtnorm_1.1-1       bibtex_0.4.2.2      metap_1.4          
##  [40] dtw_1.21-3          plotrix_3.7-8       viridisLite_0.3.0  
##  [43] htmlTable_2.0.1     tmvnsim_1.0-2       reticulate_1.16    
##  [46] foreign_0.8-70      bit_4.0.4           proxy_0.4-24       
##  [49] mclust_5.4.6        SDMTools_1.1-221.2  Formula_1.2-3      
##  [52] tsne_0.1-3          stats4_3.5.1        htmlwidgets_1.5.1  
##  [55] httr_1.4.2          gplots_3.0.1.1      RColorBrewer_1.1-2 
##  [58] fpc_2.2-7           acepack_1.4.1       TFisher_0.2.0      
##  [61] modeltools_0.2-23   ellipsis_0.3.1      Seurat_2.3.4       
##  [64] ica_1.0-2           pkgconfig_2.0.3     R.methodsS3_1.8.1  
##  [67] flexmix_2.3-15      nnet_7.3-12         reshape2_1.4.4     
##  [70] tidyselect_1.1.0    rlang_0.4.7         munsell_0.5.0      
##  [73] tools_3.5.1         mathjaxr_1.0-1      ggridges_0.5.2     
##  [76] evaluate_0.14       stringr_1.4.0       yaml_2.2.1         
##  [79] knitr_1.29          bit64_4.0.5         fitdistrplus_1.1-1 
##  [82] robustbase_0.93-6   caTools_1.17.1.1    randomForest_4.6-14
##  [85] RANN_2.6.1          pbapply_1.4-3       nlme_3.1-137       
##  [88] R.oo_1.24.0         hdf5r_1.3.3         compiler_3.5.1     
##  [91] rstudioapi_0.11     png_0.1-7           stringi_1.5.3      
##  [94] lattice_0.20-35     Matrix_1.2-14       multtest_2.38.0    
##  [97] vctrs_0.3.4         mutoss_0.1-12       pillar_1.4.6       
## [100] lifecycle_0.2.0     BiocManager_1.30.10 Rdpack_1.0.0       
## [103] lmtest_0.9-38       data.table_1.13.0   bitops_1.0-6       
## [106] irlba_2.3.3         gbRd_0.4-11         R6_2.4.1           
## [109] latticeExtra_0.6-28 renv_0.10.0         KernSmooth_2.23-15 
## [112] gridExtra_2.3       codetools_0.2-15    MASS_7.3-50        
## [115] gtools_3.8.2        assertthat_0.2.1    rprojroot_1.3-2    
## [118] withr_2.2.0         mnormt_2.0.2        multcomp_1.4-13    
## [121] diptest_0.75-7      parallel_3.5.1      doSNOW_1.0.18      
## [124] hms_0.5.3           grid_3.5.1          rpart_4.1-13       
## [127] class_7.3-14        rmarkdown_1.11      segmented_1.2-0    
## [130] Rtsne_0.15          git2r_0.27.1        numDeriv_2016.8-1.1
## [133] Biobase_2.42.0      base64enc_0.1-3
```

</details>


***

<!-- END OF END MATTER -->
