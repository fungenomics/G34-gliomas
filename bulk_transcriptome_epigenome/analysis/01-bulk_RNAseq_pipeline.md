---
title: "01 - Bulk RNAseq"
date: "14 September, 2020"
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
## Cache: ~/tmp/G34-gliomas/bulk_transcriptome_epigenome/01/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/bulk_transcriptome_epigenome/output/01
```

```
## G34-gliomas/bulk_transcriptome_epigenome/figures/01//
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
library(tidyr)
library(readr)
library(dplyr)
library(magrittr)
library(glue)
library(purrr)

# For plotting
library(ggplot2)
library(scales)
library(ggrepel)
ggplot2::theme_set(theme_min(base_size = 13))

# Custom
source(here::here(subdir, "analysis/functions.R"))
```

# Analysis

## Prepare inputs for RNAseq pipeline - cell lines

Load sample metadata:


```r
bulk_samples <- readxl::read_xlsx(
  here::here(subdir, "data/2020-02-25_bulkRNAseq_metadata.xlsx")) %>% 
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
## [1] "HGG-H3.3-K27M-pons"     "HGG-H3.3-K27M-thalamic" "HGG-IDH"               
## [4] "HGG-WT"                 "PNET-HGNET-BCOR"
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

Load sample metadata for parental cell lines:


```r
cl_samples <- readxl::read_xlsx(here::here(subdir, "data/2020-07-23_cellline_parental_metadata.xlsx"))
```


## Plot gene expression: tumours

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
pipeline_path <- "../../../2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/"

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
  rr_ggplot(aes(x = group, y = Pdgfra_normalized_expression), plot_num = 1) +
  geom_boxplot(aes(fill = group), width = 0.5) +
  scale_fill_manual(values = palette_groups) +
  geom_point() +
  scale_y_continuous(labels = comma) +
  rotateX() +
  ylab("Normalized expression") + xlab("Group") + ggtitle("PDGFRA") +
  noLegend()
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/01//bulk_RNAseq_PDGFRA_expression-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/01//bulk_RNAseq_PDGFRA_expression...*]~</span>

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

```
## Loading required package: cowplot
```

```
## 
## Attaching package: 'cowplot'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     ggsave
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/01//bulk_RNAseq_TF_expression-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/01//bulk_RNAseq_TF_expression...*]~</span>

## PDGFRA vs. GSX2 correlation in tumours & cell lines

Here we'll examine the correlation between PDGFRA and GSX2 expression in bulk RNAseq
data for human tumours and patient-derived cell lines.


```r
goi <- c("GSX2", "PDGFRA")

cor_df <- extract_pipeline_counts(file.path(pipeline_path, "G34RV_and_parental_cell_lines/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                  goi)  %>% 
  left_join(bind_rows(bulk_samples, cl_samples), by = c("sample" = "Sample")) %>%
  mutate(Group2 = paste0(Source, " - ", Genotype_PDGFRA)) %>% 
  select(sample, Group2, gene_symbol, gene_expression, Source) %>%
  spread(gene_symbol, gene_expression)
```

```
## Joining, by = "gene_symbol"
```



```r
# Linear model for G34
m <- lm(PDGFRA ~ GSX2, cor_df)

# Check the model
summary(m)
```

```
## 
## Call:
## lm(formula = PDGFRA ~ GSX2, data = cor_df)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -154320  -25556  -15324   32991  236256 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 26239.18   21681.00   1.210 0.239031    
## GSX2           63.00      13.91   4.529 0.000166 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 80730 on 22 degrees of freedom
## Multiple R-squared:  0.4825,	Adjusted R-squared:  0.459 
## F-statistic: 20.51 on 1 and 22 DF,  p-value: 0.0001658
```

```r
cor <- cor_df %>% select(PDGFRA, GSX2) %>% cor() %>% .[2,1] %>% round(2)
rsq <- round(summary(m)$r.squared, 2)
pv  <- round(summary(m)$coefficients[2,4], 4)

# Generate the basic plot
cor_df %>%
  # To avoid errors on log-scale, add a small pseudocount for samples with no Gsx2 expression
  mutate(GSX2 = ifelse(GSX2 == 0, 1, GSX2)) %>% 
  rr_ggplot(aes(x = GSX2, y = PDGFRA), plot_num = 1) +
  geom_point(aes(shape = Group2), size = 2, stroke = 1, colour = "cyan4", fill = "cyan4") +
  # Set the shape to be able to distinguish sample groups
  scale_shape_manual(values = c("Tumor - Mutant"     = "circle open",
                                "Tumor - WT"         = "circle filled",
                                "Cell line - Mutant" = "triangle open",
                                "Cell line - WT"     = "triangle filled")) +
  geom_smooth(data = cor_df, method = "lm", se = FALSE, colour = "cyan4", alpha = 0.5, size = 0.5) +
  # Add the stats
  annotate(geom = "text", label = glue("cor {cor}, R^2 {rsq}, p-value {pv}"), x = 10^2, y = 10^5.9) +
  # Log-transform the scale
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(1, 10^4)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(1, 10^6)) +
  # Zoom into the area of the plot where data is contained
  coord_cartesian(ylim = c(10^3, 10^6))
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/01//gsx2_pdgfra_correlation-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/01//gsx2_pdgfra_correlation...*]~</span>



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-09-14 21:38:19
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [87deb0b] 2020-09-15: Rename bulk folder to be more descriptive for tx/epi analyses
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
##  [1] cowplot_0.9.4 ggrepel_0.8.0 scales_1.1.1  ggplot2_3.1.0 purrr_0.3.4  
##  [6] glue_1.4.2    magrittr_1.5  dplyr_0.8.0   readr_1.3.1   tidyr_0.8.2  
## [11] here_0.1     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.5          git2r_0.27.1        pillar_1.4.6       
##  [4] compiler_3.5.1      BiocManager_1.30.10 plyr_1.8.6         
##  [7] tools_3.5.1         digest_0.6.25       evaluate_0.14      
## [10] lifecycle_0.2.0     tibble_3.0.3        gtable_0.3.0       
## [13] pkgconfig_2.0.3     rlang_0.4.7         yaml_2.2.1         
## [16] xfun_0.17           withr_2.2.0         stringr_1.4.0      
## [19] knitr_1.29          vctrs_0.3.4         hms_0.5.3          
## [22] rprojroot_1.3-2     grid_3.5.1          tidyselect_1.1.0   
## [25] R6_2.4.1            rmarkdown_1.11      reshape2_1.4.4     
## [28] farver_2.0.3        codetools_0.2-15    backports_1.1.9    
## [31] ellipsis_0.3.1      htmltools_0.5.0     assertthat_0.2.1   
## [34] colorspace_1.4-1    renv_0.10.0         labeling_0.3       
## [37] stringi_1.5.3       lazyeval_0.2.2      munsell_0.5.0      
## [40] crayon_1.3.4
```

</details>


***

<!-- END OF END MATTER -->
