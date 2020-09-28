---
title: "01 - Bulk RNAseq"
date: "28 September, 2020"
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

## Prepare inputs for RNAseq pipeline - patient tumors

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

### G34 mutants

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
p1 <- cor_df %>%
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

p1
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/01//gsx2_pdgfra_correlation-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/01//gsx2_pdgfra_correlation...*]~</span>

<!-- Generate a version with labels: -->

<!-- ```{r cor_plot_lab, fig.width = 12, fig.height = 9} -->

<!-- p1 + -->
<!--     geom_label_repel(aes(label = sample)) -->

<!-- ``` -->

### Non-G34

Repeat this computation for non-G34 samples, as a comparator.


```r
cor_df_non_g34 <- counts_subset %>%
  mutate(Group2 = recode(Group2, "K27M-pons" = "K27M", "K27M-thal" = "K27M")) %>% 
  filter(gene_symbol %in% c("GSX2", "PDGFRA")) %>% 
  filter(Group2 != "G34R/V" & Source != "Cell line") %>% 
  select(sample, Group2, Source, gene_symbol, gene_expression) %>% 
  spread(gene_symbol, gene_expression)

# Linear model for non-G34
m2 <- lm(PDGFRA ~ GSX2, cor_df_non_g34)
summary(m2)
```

```
## 
## Call:
## lm(formula = PDGFRA ~ GSX2, data = cor_df_non_g34)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -102237  -16944   -6647    4753  386929 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)   4953.5    15358.1   0.323   0.7490   
## GSX2           699.6      237.5   2.945   0.0057 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 74700 on 35 degrees of freedom
## Multiple R-squared:  0.1986,	Adjusted R-squared:  0.1757 
## F-statistic: 8.676 on 1 and 35 DF,  p-value: 0.0057
```

```r
# R-squared
round(summary(m2)$r.squared, 2)
```

```
## [1] 0.2
```

```r
# p-value
round(summary(m2)$coefficients[2,4], 3)
```

```
## [1] 0.006
```

```r
# Correlation value
sqrt(summary(m2)$r.squared)
```

```
## [1] 0.4456931
```

Let's also perform the calculation per-group:


```r
map_dfr(unique(cor_df_non_g34$Group2), function(grp) {
  
  m <- lm(PDGFRA ~ GSX2, cor_df_non_g34 %>% filter(Group2 == grp))
  data.frame("Group"   = grp,
             "R2"      = round(summary(m)$r.squared, 2),
             "p-value" = round(summary(m)$coefficients[2,4], 3),
             "Correlation" = sqrt(summary(m)$r.squared))
}
)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Group"],"name":[1],"type":["chr"],"align":["left"]},{"label":["R2"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Correlation"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"HGNET-BCOR","2":"0.35","3":"0.297","4":"0.5881594"},{"1":"WT","2":"0.08","3":"0.457","4":"0.2853009"},{"1":"IDH","2":"0.01","3":"0.758","4":"0.1119512"},{"1":"K27M","2":"0.26","3":"0.077","4":"0.5066390"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Extract DGE



```r
read.table(file.path(pipeline_path, "HGG-G34R.V_vs_HGG-IDH_batch_covariate/diff/Ensembl.ensGene.exon/HGG-IDHvsHGG-G34R.V.tsv"), header = T, sep = "\t", check.names = FALSE) %>%
  filter(padj < 0.01 & baseMean > 100 & abs(log2FoldChange) > 2) %>%
  arrange(desc(log2FoldChange)) %T>%
  rr_write_tsv(path = glue("{out}/DGE_HGG-IDHvsHGG-G34R.V_padj<0.01_baseMean>100_absLFC>2.tsv"),
               desc = "Differentially expressed genes for G34R/V vs IDH patient tumors, filtered to have adjusted p-val <0.01, expression baseMean >100, and absolute log2 fold change >2")
```

```
## ...writing description of DGE_HGG-IDHvsHGG-G34R.V_padj<0.01_baseMean>100_absLFC>2.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/01/DGE_HGG-IDHvsHGG-G34R.V_padj<0.01_baseMean>100_absLFC>2.desc
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ID"],"name":[1],"type":["fctr"],"align":["left"]},{"label":["baseMean"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["log2FoldChange"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["lfcSE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["stat"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["pvalue"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["padj"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"ENSG00000230386:RP13-157F18.2","2":"124.3063","3":"9.036575","4":"1.9208618","5":"4.704438","6":"2.545666e-06","7":"2.551452e-04"},{"1":"ENSG00000237735:AP000473.6","2":"166.9447","3":"8.458577","4":"1.0360071","5":"8.164593","6":"3.225195e-16","7":"2.702391e-12"},{"1":"ENSG00000095777:MYO3A","2":"2005.8312","3":"7.670394","4":"0.9619065","5":"7.974157","6":"1.534238e-15","7":"9.182416e-12"},{"1":"ENSG00000079689:SCGN","2":"639.2007","3":"7.338952","4":"1.3227093","5":"5.548424","6":"2.882554e-08","7":"8.848290e-06"},{"1":"ENSG00000185686:PRAME","2":"429.7103","3":"7.149183","4":"1.2123048","5":"5.897183","6":"3.697593e-09","7":"1.702315e-06"},{"1":"ENSG00000124253:PCK1","2":"496.9496","3":"6.990020","4":"0.7635237","5":"9.154948","6":"5.438390e-20","7":"1.139207e-15"},{"1":"ENSG00000197320:AC060834.2","2":"191.5382","3":"6.989506","4":"1.2537353","5":"5.574946","6":"2.476070e-08","7":"7.979613e-06"},{"1":"ENSG00000115844:DLX2","2":"2036.7102","3":"6.796860","4":"0.8194817","5":"8.294096","6":"1.094139e-16","7":"1.145974e-12"},{"1":"ENSG00000144355:DLX1","2":"3996.0621","3":"6.626103","4":"0.7773955","5":"8.523465","6":"1.548499e-17","7":"2.162478e-13"},{"1":"ENSG00000229108:AC005550.4","2":"142.7989","3":"6.545081","4":"1.2744474","5":"5.135623","6":"2.812107e-07","7":"5.151903e-05"},{"1":"ENSG00000180613:GSX2","2":"751.7376","3":"6.351361","4":"0.8644190","5":"7.347549","6":"2.018736e-13","7":"4.027379e-10"},{"1":"ENSG00000164651:SP8","2":"592.8037","3":"6.166849","4":"0.8248334","5":"7.476478","6":"7.634070e-14","7":"1.881349e-10"},{"1":"ENSG00000206588:RNU1-28P","2":"180.9785","3":"6.052195","4":"0.9685376","5":"6.248797","6":"4.136260e-10","7":"2.794978e-07"},{"1":"ENSG00000140798:ABCC12","2":"234.5112","3":"5.986998","4":"1.0058208","5":"5.952350","6":"2.643188e-09","7":"1.350443e-06"},{"1":"ENSG00000106511:MEOX2","2":"1014.2195","3":"5.982819","4":"1.1142680","5":"5.369282","6":"7.905060e-08","7":"1.967982e-05"},{"1":"ENSG00000234710:AC060834.3","2":"102.8190","3":"5.602554","4":"1.2376773","5":"4.526667","6":"5.992111e-06","7":"4.903115e-04"},{"1":"ENSG00000079931:MOXD1","2":"7883.4561","3":"5.577695","4":"0.7279413","5":"7.662287","6":"1.826514e-14","7":"7.652180e-11"},{"1":"ENSG00000235688:AC116614.1","2":"139.5816","3":"5.503076","4":"1.1733396","5":"4.690097","6":"2.730756e-06","7":"2.684075e-04"},{"1":"ENSG00000142700:DMRTA2","2":"632.4912","3":"5.311144","4":"0.7263135","5":"7.312468","6":"2.622799e-13","7":"4.994644e-10"},{"1":"ENSG00000165078:CPA6","2":"225.1601","3":"5.188138","4":"0.9223114","5":"5.625148","6":"1.853488e-08","7":"6.417510e-06"},{"1":"ENSG00000256637:RP11-76I14.1","2":"324.0321","3":"5.174964","4":"0.9154519","5":"5.652906","6":"1.577577e-08","7":"5.648938e-06"},{"1":"ENSG00000253642:RP11-317N12.1","2":"123.2810","3":"5.113232","4":"0.6927080","5":"7.381511","6":"1.565033e-13","7":"3.642615e-10"},{"1":"ENSG00000180251:SLC9A4","2":"141.3534","3":"5.042366","4":"0.7131662","5":"7.070394","6":"1.544940e-12","7":"2.080340e-09"},{"1":"ENSG00000269155:AL009178.1","2":"154.5917","3":"5.016693","4":"0.9988128","5":"5.022656","6":"5.096190e-07","7":"7.792150e-05"},{"1":"ENSG00000172554:SNTG2","2":"367.3406","3":"4.960829","4":"0.8442181","5":"5.876241","6":"4.196873e-09","7":"1.866784e-06"},{"1":"ENSG00000119913:TECTB","2":"108.7082","3":"4.927397","4":"0.9957767","5":"4.948296","6":"7.486613e-07","7":"1.052522e-04"},{"1":"ENSG00000225684:FAM225B","2":"1651.4370","3":"4.916800","4":"0.6273744","5":"7.837106","6":"4.610473e-15","7":"2.414447e-11"},{"1":"ENSG00000231764:DLX6-AS1","2":"13137.4337","3":"4.875184","4":"0.8941102","5":"5.452554","6":"4.965141e-08","7":"1.333427e-05"},{"1":"ENSG00000159212:CLIC6","2":"172.2811","3":"4.840158","4":"0.9807654","5":"4.935083","6":"8.011683e-07","7":"1.107754e-04"},{"1":"ENSG00000138675:FGF5","2":"102.0546","3":"4.807487","4":"1.0296315","5":"4.669133","6":"3.024729e-06","7":"2.847664e-04"},{"1":"ENSG00000176076:KCNE1L","2":"143.5460","3":"4.801079","4":"0.8827226","5":"5.438945","6":"5.359708e-08","7":"1.412233e-05"},{"1":"ENSG00000144407:PTH2R","2":"293.4048","3":"4.748359","4":"0.9903670","5":"4.794545","6":"1.630449e-06","7":"1.816693e-04"},{"1":"ENSG00000234685:NUS1P2","2":"110.5975","3":"4.732100","4":"0.6599377","5":"7.170525","6":"7.471062e-13","7":"1.159260e-09"},{"1":"ENSG00000115507:OTX1","2":"976.4566","3":"4.614460","4":"0.6105138","5":"7.558322","6":"4.083036e-14","7":"1.425490e-10"},{"1":"ENSG00000183671:GPR1","2":"174.9171","3":"4.612630","4":"0.9582584","5":"4.813555","6":"1.482688e-06","7":"1.687968e-04"},{"1":"ENSG00000170381:SEMA3E","2":"1531.0918","3":"4.523987","4":"0.9065074","5":"4.990568","6":"6.020195e-07","7":"8.849686e-05"},{"1":"ENSG00000008196:TFAP2B","2":"157.3933","3":"4.460556","4":"1.1671616","5":"3.821713","6":"1.325279e-04","7":"5.038346e-03"},{"1":"ENSG00000231528:FAM225A","2":"140.9892","3":"4.445208","4":"0.6700680","5":"6.633965","6":"3.267868e-11","7":"3.183892e-08"},{"1":"ENSG00000206585:RNVU1-9","2":"163.8900","3":"4.427262","4":"0.8621134","5":"5.135359","6":"2.816054e-07","7":"5.151903e-05"},{"1":"ENSG00000004848:ARX","2":"782.1428","3":"4.368118","4":"0.8110767","5":"5.385579","6":"7.221172e-08","7":"1.833521e-05"},{"1":"ENSG00000210049:MT-TF","2":"147.7714","3":"4.298040","4":"0.8841616","5":"4.861148","6":"1.167069e-06","7":"1.429659e-04"},{"1":"ENSG00000237070:AC005550.3","2":"197.3844","3":"4.290546","4":"0.9763159","5":"4.394629","6":"1.109620e-05","7":"7.911305e-04"},{"1":"ENSG00000269416:CTB-175P5.4","2":"103.0627","3":"4.267214","4":"0.6103854","5":"6.991016","6":"2.729030e-12","7":"3.464628e-09"},{"1":"ENSG00000166863:TAC3","2":"145.2942","3":"4.206931","4":"0.8648434","5":"4.864385","6":"1.148132e-06","7":"1.423106e-04"},{"1":"ENSG00000221716:SNORA11","2":"444.3235","3":"4.125925","4":"0.9191258","5":"4.488967","6":"7.156953e-06","7":"5.657369e-04"},{"1":"ENSG00000197977:ELOVL2","2":"6015.1216","3":"4.119925","4":"0.7220729","5":"5.705691","6":"1.158717e-08","7":"4.334327e-06"},{"1":"ENSG00000262619:LINC00621","2":"128.2420","3":"4.101435","4":"0.7018823","5":"5.843480","6":"5.112136e-09","7":"2.141729e-06"},{"1":"ENSG00000186007:LEMD1","2":"114.2655","3":"4.092471","4":"0.7495637","5":"5.459804","6":"4.766612e-08","7":"1.296735e-05"},{"1":"ENSG00000113494:PRLR","2":"705.7420","3":"4.082499","4":"0.8375859","5":"4.874126","6":"1.092914e-06","7":"1.377486e-04"},{"1":"ENSG00000122585:NPY","2":"1520.4322","3":"4.028783","4":"0.8318011","5":"4.843445","6":"1.276070e-06","7":"1.505943e-04"},{"1":"ENSG00000144681:STAC","2":"312.1900","3":"3.964790","4":"0.8634246","5":"4.591936","6":"4.391538e-06","7":"3.857096e-04"},{"1":"ENSG00000200816:SNORA38","2":"329.8743","3":"3.951814","4":"0.6395962","5":"6.178608","6":"6.466937e-10","7":"4.043766e-07"},{"1":"ENSG00000166473:PKD1L2","2":"340.4059","3":"3.941925","4":"0.9461844","5":"4.166128","6":"3.098169e-05","7":"1.729142e-03"},{"1":"ENSG00000206634:SNORA22","2":"667.4609","3":"3.929436","4":"0.7188157","5":"5.466541","6":"4.589016e-08","7":"1.273961e-05"},{"1":"ENSG00000223561:AC003090.1","2":"174.8482","3":"3.919037","4":"0.8128776","5":"4.821190","6":"1.427043e-06","7":"1.637972e-04"},{"1":"ENSG00000197415:VEPH1","2":"1870.4426","3":"3.909740","4":"0.7036861","5":"5.556086","6":"2.758914e-08","7":"8.620745e-06"},{"1":"ENSG00000207344:SNORA22","2":"147.5221","3":"3.909257","4":"0.7037508","5":"5.554888","6":"2.777898e-08","7":"8.620745e-06"},{"1":"ENSG00000137441:FGFBP2","2":"451.3922","3":"3.898520","4":"0.7960812","5":"4.897139","6":"9.724222e-07","7":"1.265206e-04"},{"1":"ENSG00000271956:DLX6-AS2","2":"216.0649","3":"3.892941","4":"0.9559779","5":"4.072208","6":"4.656961e-05","7":"2.322659e-03"},{"1":"ENSG00000167588:GPD1","2":"2571.9178","3":"3.891544","4":"0.7721567","5":"5.039836","6":"4.659298e-07","7":"7.478977e-05"},{"1":"ENSG00000108001:EBF3","2":"534.4927","3":"3.847947","4":"0.8183617","5":"4.702013","6":"2.576090e-06","7":"2.575782e-04"},{"1":"ENSG00000006377:DLX6","2":"940.9403","3":"3.821796","4":"0.9730523","5":"3.927637","6":"8.578464e-05","7":"3.667293e-03"},{"1":"ENSG00000221475:SNORA11D","2":"321.8517","3":"3.814320","4":"0.7254309","5":"5.258006","6":"1.456256e-07","7":"3.177596e-05"},{"1":"ENSG00000207513:RNU1-3","2":"52086.6861","3":"3.801542","4":"0.8277811","5":"4.592449","6":"4.380754e-06","7":"3.855708e-04"},{"1":"ENSG00000206652:RNU1-1","2":"52211.3483","3":"3.779214","4":"0.8114647","5":"4.657275","6":"3.204227e-06","7":"2.983135e-04"},{"1":"ENSG00000130751:NPAS1","2":"304.1136","3":"3.762528","4":"0.8207522","5":"4.584244","6":"4.556330e-06","7":"3.976822e-04"},{"1":"ENSG00000145248:SLC10A4","2":"1017.2632","3":"3.720302","4":"0.6323013","5":"5.883750","6":"4.010747e-09","7":"1.806777e-06"},{"1":"ENSG00000201643:SNORA14A","2":"143.7889","3":"3.715696","4":"0.8563513","5":"4.338985","6":"1.431421e-05","7":"9.578310e-04"},{"1":"ENSG00000187475:HIST1H1T","2":"241.0682","3":"3.680121","4":"0.8277360","5":"4.446008","6":"8.748079e-06","7":"6.603617e-04"},{"1":"ENSG00000105825:TFPI2","2":"615.7337","3":"3.663801","4":"0.9808430","5":"3.735359","6":"1.874475e-04","7":"6.525236e-03"},{"1":"ENSG00000133863:TEX15","2":"564.9442","3":"3.657717","4":"0.7061780","5":"5.179596","6":"2.223669e-07","7":"4.333052e-05"},{"1":"ENSG00000188987:HIST1H4D","2":"2708.9938","3":"3.638894","4":"0.5757547","5":"6.320216","6":"2.611979e-10","7":"1.854727e-07"},{"1":"ENSG00000271167:LINC01109","2":"201.7779","3":"3.619977","4":"0.8803585","5":"4.111935","6":"3.923564e-05","7":"2.026852e-03"},{"1":"ENSG00000269558:AC104057.1","2":"120.9082","3":"3.583236","4":"0.9088553","5":"3.942581","6":"8.060947e-05","7":"3.495998e-03"},{"1":"ENSG00000112183:RBM24","2":"210.3410","3":"3.580412","4":"0.6858640","5":"5.220295","6":"1.786387e-07","7":"3.615491e-05"},{"1":"ENSG00000237807:RP11-400K9.4","2":"179.0465","3":"3.560735","4":"0.9127533","5":"3.901092","6":"9.575959e-05","7":"3.948669e-03"},{"1":"ENSG00000270866:RP11-706J10.2","2":"173.1928","3":"3.555561","4":"0.8957097","5":"3.969546","6":"7.200968e-05","7":"3.219686e-03"},{"1":"ENSG00000136750:GAD2","2":"1477.7334","3":"3.553368","4":"0.9252852","5":"3.840295","6":"1.228867e-04","7":"4.749390e-03"},{"1":"ENSG00000200434:RNA5-8SP2","2":"501.5696","3":"3.546196","4":"0.7302167","5":"4.856362","6":"1.195622e-06","7":"1.456121e-04"},{"1":"ENSG00000221676:RNU6ATAC","2":"267.2400","3":"3.546123","4":"0.6830374","5":"5.191696","6":"2.083873e-07","7":"4.101647e-05"},{"1":"ENSG00000165105:RASEF","2":"241.6959","3":"3.539092","4":"0.7553151","5":"4.685583","6":"2.791639e-06","7":"2.714851e-04"},{"1":"ENSG00000091490:SEL1L3","2":"4947.7448","3":"3.515447","4":"0.5064180","5":"6.941788","6":"3.871668e-12","7":"4.634387e-09"},{"1":"ENSG00000206965:RNU6-5P","2":"205.1958","3":"3.490845","4":"0.8311179","5":"4.200181","6":"2.667020e-05","7":"1.554031e-03"},{"1":"ENSG00000258932:CTD-3006G17.2","2":"186.3873","3":"3.469834","4":"0.6906328","5":"5.024137","6":"5.057022e-07","7":"7.792150e-05"},{"1":"ENSG00000203867:RBM20","2":"454.2318","3":"3.466192","4":"0.6216696","5":"5.575617","6":"2.466540e-08","7":"7.979613e-06"},{"1":"ENSG00000207493:SNORA46","2":"159.0710","3":"3.432386","4":"0.7423551","5":"4.623645","6":"3.770553e-06","7":"3.368173e-04"},{"1":"ENSG00000177283:FZD8","2":"943.0419","3":"3.421624","4":"0.6617359","5":"5.170680","6":"2.332442e-07","7":"4.468485e-05"},{"1":"ENSG00000151882:CCL28","2":"116.2338","3":"3.387484","4":"0.8904234","5":"3.804352","6":"1.421760e-04","7":"5.294636e-03"},{"1":"ENSG00000153404:PLEKHG4B","2":"342.3394","3":"3.383045","4":"0.8919133","5":"3.793020","6":"1.488259e-04","7":"5.505951e-03"},{"1":"ENSG00000138759:FRAS1","2":"1047.5739","3":"3.374007","4":"0.8792806","5":"3.837236","6":"1.244268e-04","7":"4.795639e-03"},{"1":"ENSG00000258548:LINC00645","2":"1135.0018","3":"3.361091","4":"0.7013443","5":"4.792355","6":"1.648351e-06","7":"1.822102e-04"},{"1":"ENSG00000135824:RGS8","2":"990.1227","3":"3.356990","4":"0.8075035","5":"4.157245","6":"3.221091e-05","7":"1.761035e-03"},{"1":"ENSG00000156510:HKDC1","2":"310.8280","3":"3.352334","4":"0.6327340","5":"5.298172","6":"1.169680e-07","7":"2.677801e-05"},{"1":"ENSG00000261185:RP11-317P15.5","2":"642.3084","3":"3.343505","4":"0.8270559","5":"4.042660","6":"5.284833e-05","7":"2.541999e-03"},{"1":"ENSG00000164142:FAM160A1","2":"104.3119","3":"3.241987","4":"0.8559870","5":"3.787425","6":"1.522166e-04","7":"5.584162e-03"},{"1":"ENSG00000197057:DTHD1","2":"489.8414","3":"3.219929","4":"0.7864476","5":"4.094270","6":"4.235011e-05","7":"2.161000e-03"},{"1":"ENSG00000169515:CCDC8","2":"296.5953","3":"3.212904","4":"0.6276918","5":"5.118601","6":"3.078098e-07","7":"5.487528e-05"},{"1":"ENSG00000180592:SKIDA1","2":"454.5368","3":"3.179393","4":"0.4218431","5":"7.536909","6":"4.812435e-14","7":"1.550900e-10"},{"1":"ENSG00000186960:C14orf23","2":"591.6774","3":"3.172566","4":"0.4463771","5":"7.107366","6":"1.182785e-12","7":"1.708716e-09"},{"1":"ENSG00000074219:TEAD2","2":"1370.9023","3":"3.137649","4":"0.4560940","5":"6.879392","6":"6.010858e-12","7":"6.806079e-09"},{"1":"ENSG00000112333:NR2E1","2":"1442.4843","3":"3.121038","4":"0.4315342","5":"7.232424","6":"4.744472e-13","7":"7.950786e-10"},{"1":"ENSG00000144152:FBLN7","2":"660.4919","3":"3.090439","4":"0.4119387","5":"7.502181","6":"6.276467e-14","7":"1.753017e-10"},{"1":"ENSG00000202538:RNU4-2","2":"17222.8981","3":"3.085751","4":"0.6868300","5":"4.492744","6":"7.031129e-06","7":"5.589547e-04"},{"1":"ENSG00000228065:RP11-222A11.1","2":"256.8652","3":"3.076268","4":"0.6453683","5":"4.766686","6":"1.872806e-06","7":"2.011825e-04"},{"1":"ENSG00000177133:LINC00982","2":"484.8924","3":"3.065205","4":"0.6626791","5":"4.625474","6":"3.737420e-06","7":"3.360069e-04"},{"1":"ENSG00000256316:HIST1H3F","2":"2025.7403","3":"3.033389","4":"0.6441904","5":"4.708839","6":"2.491321e-06","7":"2.521108e-04"},{"1":"ENSG00000226684:RP3-418C23.2","2":"187.7430","3":"2.992337","4":"0.4747945","5":"6.302383","6":"2.931033e-10","7":"2.046594e-07"},{"1":"ENSG00000251002:AE000661.37","2":"146.4731","3":"2.990781","4":"0.8113108","5":"3.686356","6":"2.274877e-04","7":"7.405281e-03"},{"1":"ENSG00000149451:ADAM33","2":"538.4450","3":"2.979518","4":"0.7448884","5":"3.999952","6":"6.335523e-05","7":"2.923202e-03"},{"1":"ENSG00000239183:SNORA84","2":"135.1077","3":"2.974613","4":"0.6231283","5":"4.773677","6":"1.808925e-06","7":"1.958267e-04"},{"1":"ENSG00000111186:WNT5B","2":"283.5161","3":"2.963797","4":"0.5037063","5":"5.883979","6":"4.005188e-09","7":"1.806777e-06"},{"1":"ENSG00000162873:KLHDC8A","2":"10099.7079","3":"2.961924","4":"0.6569006","5":"4.508937","6":"6.515315e-06","7":"5.219104e-04"},{"1":"ENSG00000007062:PROM1","2":"1944.5568","3":"2.936248","4":"0.5616905","5":"5.227520","6":"1.717990e-07","7":"3.510985e-05"},{"1":"ENSG00000164649:CDCA7L","2":"3993.3334","3":"2.932198","4":"0.4900807","5":"5.983093","6":"2.189396e-09","7":"1.206904e-06"},{"1":"ENSG00000173040:EVC2","2":"432.5603","3":"2.929803","4":"0.7029366","5":"4.167947","6":"3.073547e-05","7":"1.721474e-03"},{"1":"ENSG00000256128:LINC00944","2":"115.1050","3":"2.921666","4":"0.8119918","5":"3.598147","6":"3.204919e-04","7":"9.526296e-03"},{"1":"ENSG00000272235:RP11-22L13.1","2":"174.8014","3":"2.921209","4":"0.6047062","5":"4.830790","6":"1.359921e-06","7":"1.582608e-04"},{"1":"ENSG00000221420:SNORA81","2":"1284.9420","3":"2.904979","4":"0.7818869","5":"3.715344","6":"2.029275e-04","7":"6.893501e-03"},{"1":"ENSG00000075275:CELSR1","2":"2583.2977","3":"2.890494","4":"0.5907535","5":"4.892893","6":"9.936446e-07","7":"1.280884e-04"},{"1":"ENSG00000117069:ST6GALNAC5","2":"345.3144","3":"2.887394","4":"0.7757033","5":"3.722292","6":"1.974224e-04","7":"6.760821e-03"},{"1":"ENSG00000125869:LAMP5","2":"1018.9278","3":"2.887196","4":"0.5839706","5":"4.944078","6":"7.650500e-07","7":"1.071966e-04"},{"1":"ENSG00000074047:GLI2","2":"1455.7537","3":"2.865226","4":"0.4810124","5":"5.956656","6":"2.574517e-09","7":"1.350443e-06"},{"1":"ENSG00000250740:RP11-710F7.2","2":"114.9966","3":"2.847521","4":"0.6864846","5":"4.147975","6":"3.354297e-05","7":"1.814391e-03"},{"1":"ENSG00000082497:SERTAD4","2":"194.1412","3":"2.834870","4":"0.6035022","5":"4.697365","6":"2.635393e-06","7":"2.610160e-04"},{"1":"ENSG00000232270:INTS4L2","2":"676.0662","3":"2.834013","4":"0.6741882","5":"4.203593","6":"2.627113e-05","7":"1.537192e-03"},{"1":"ENSG00000139445:FOXN4","2":"124.6480","3":"2.826444","4":"0.6387147","5":"4.425206","6":"9.635020e-06","7":"7.169790e-04"},{"1":"ENSG00000272796:RP1-74M1.3","2":"868.9933","3":"2.816671","4":"0.5624946","5":"5.007463","6":"5.515201e-07","7":"8.311488e-05"},{"1":"ENSG00000196966:HIST1H3E","2":"414.2804","3":"2.808804","4":"0.6349976","5":"4.423330","6":"9.719118e-06","7":"7.206769e-04"},{"1":"ENSG00000211459:MT-RNR1","2":"95620.3184","3":"2.801717","4":"0.7206277","5":"3.887883","6":"1.011222e-04","7":"4.089299e-03"},{"1":"ENSG00000015413:DPEP1","2":"163.8005","3":"2.796864","4":"0.7429998","5":"3.764287","6":"1.670250e-04","7":"6.006448e-03"},{"1":"ENSG00000133216:EPHB2","2":"4015.4897","3":"2.740162","4":"0.4691977","5":"5.840102","6":"5.216900e-09","7":"2.163980e-06"},{"1":"ENSG00000201457:SNORA55","2":"166.8783","3":"2.731900","4":"0.7477871","5":"3.653313","6":"2.588785e-04","7":"8.179271e-03"},{"1":"ENSG00000264940:SNORD3C","2":"58041.0533","3":"2.715462","4":"0.7211373","5":"3.765527","6":"1.661983e-04","7":"5.981852e-03"},{"1":"ENSG00000172878:METAP1D","2":"2317.2241","3":"2.709265","4":"0.5623449","5":"4.817800","6":"1.451497e-06","7":"1.656961e-04"},{"1":"ENSG00000135451:TROAP","2":"845.0542","3":"2.703996","4":"0.4670694","5":"5.789280","6":"7.068864e-09","7":"2.875243e-06"},{"1":"ENSG00000149571:KIRREL3","2":"1042.6759","3":"2.696906","4":"0.6199909","5":"4.349912","6":"1.361919e-05","7":"9.211221e-04"},{"1":"ENSG00000237819:AC002454.1","2":"109.3229","3":"2.688832","4":"0.6659513","5":"4.037580","6":"5.400534e-05","7":"2.585776e-03"},{"1":"ENSG00000264014:MIR4253","2":"141.0914","3":"2.682638","4":"0.4849833","5":"5.531403","6":"3.176793e-08","7":"9.506553e-06"},{"1":"ENSG00000257056:RP11-966I7.2","2":"289.2775","3":"2.678814","4":"0.4220325","5":"6.347411","6":"2.189693e-10","7":"1.609424e-07"},{"1":"ENSG00000139155:SLCO1C1","2":"1875.9125","3":"2.676306","4":"0.5542879","5":"4.828368","6":"1.376567e-06","7":"1.597542e-04"},{"1":"ENSG00000135476:ESPL1","2":"985.0336","3":"2.671527","4":"0.4053800","5":"6.590180","6":"4.392936e-11","7":"4.182774e-08"},{"1":"ENSG00000163009:C2orf48","2":"113.2164","3":"2.665428","4":"0.6233025","5":"4.276299","6":"1.900263e-05","7":"1.182935e-03"},{"1":"ENSG00000123201:GUCY1B2","2":"181.0850","3":"2.655277","4":"0.7136428","5":"3.720737","6":"1.986420e-04","7":"6.782484e-03"},{"1":"ENSG00000222000:AC092675.3","2":"241.5844","3":"2.632205","4":"0.4374359","5":"6.017351","6":"1.772946e-09","7":"1.003751e-06"},{"1":"ENSG00000101670:LIPG","2":"1897.8813","3":"2.631680","4":"0.6593185","5":"3.991516","6":"6.565234e-05","7":"2.999460e-03"},{"1":"ENSG00000256879:RP11-284H19.1","2":"150.4771","3":"2.612928","4":"0.3421416","5":"7.636979","6":"2.223780e-14","7":"8.469568e-11"},{"1":"ENSG00000117650:NEK2","2":"504.5604","3":"2.598944","4":"0.5265711","5":"4.935600","6":"7.990460e-07","7":"1.107754e-04"},{"1":"ENSG00000176165:FOXG1","2":"7968.9282","3":"2.594923","4":"0.4381816","5":"5.922027","6":"3.179980e-09","7":"1.549131e-06"},{"1":"ENSG00000210082:MT-RNR2","2":"449763.7534","3":"2.590708","4":"0.6908317","5":"3.750130","6":"1.767430e-04","7":"6.222392e-03"},{"1":"ENSG00000155760:FZD7","2":"2533.7162","3":"2.588263","4":"0.5967620","5":"4.337177","6":"1.443245e-05","7":"9.597576e-04"},{"1":"ENSG00000129173:E2F8","2":"362.5440","3":"2.579272","4":"0.4919936","5":"5.242491","6":"1.584235e-07","7":"3.352098e-05"},{"1":"ENSG00000225091:SNORA71A","2":"208.1908","3":"2.566752","4":"0.5563481","5":"4.613572","6":"3.958068e-06","7":"3.513204e-04"},{"1":"ENSG00000071909:MYO3B","2":"244.6404","3":"2.560916","4":"0.6917873","5":"3.701884","6":"2.140045e-04","7":"7.144000e-03"},{"1":"ENSG00000197459:HIST1H2BH","2":"842.7395","3":"2.549236","4":"0.6846055","5":"3.723658","6":"1.963571e-04","7":"6.737411e-03"},{"1":"ENSG00000111206:FOXM1","2":"2653.8815","3":"2.541667","4":"0.3970819","5":"6.400864","6":"1.545002e-10","7":"1.198664e-07"},{"1":"ENSG00000233639:LINC01158","2":"1913.3494","3":"2.540532","4":"0.5234183","5":"4.853732","6":"1.211594e-06","7":"1.467044e-04"},{"1":"ENSG00000200795:RNU4-1","2":"5301.0159","3":"2.529725","4":"0.6357571","5":"3.979075","6":"6.918394e-05","7":"3.130088e-03"},{"1":"ENSG00000063127:SLC6A16","2":"147.1219","3":"2.522272","4":"0.4489528","5":"5.618123","6":"1.930433e-08","7":"6.522217e-06"},{"1":"ENSG00000272511:RP11-180N14.1","2":"270.5165","3":"2.486100","4":"0.6811700","5":"3.649750","6":"2.624952e-04","7":"8.268599e-03"},{"1":"ENSG00000182968:SOX1","2":"1259.4969","3":"2.476572","4":"0.6819620","5":"3.631540","6":"2.817354e-04","7":"8.646997e-03"},{"1":"ENSG00000117399:CDC20","2":"896.8744","3":"2.475797","4":"0.4177855","5":"5.926001","6":"3.104007e-09","7":"1.540291e-06"},{"1":"ENSG00000165891:E2F7","2":"751.8563","3":"2.459589","4":"0.4500534","5":"5.465105","6":"4.626322e-08","7":"1.275130e-05"},{"1":"ENSG00000197846:HIST1H2BF","2":"2127.7195","3":"2.456072","4":"0.4863934","5":"5.049558","6":"4.428328e-07","7":"7.218864e-05"},{"1":"ENSG00000207165:SNORA70","2":"1204.9984","3":"2.440405","4":"0.5878811","5":"4.151188","6":"3.307538e-05","7":"1.794939e-03"},{"1":"ENSG00000145632:PLK2","2":"2365.4498","3":"2.434710","4":"0.5922637","5":"4.110854","6":"3.941981e-05","7":"2.031356e-03"},{"1":"ENSG00000151789:ZNF385D","2":"2131.3921","3":"2.433049","4":"0.6704732","5":"3.628854","6":"2.846827e-04","7":"8.712040e-03"},{"1":"ENSG00000168658:VWA3B","2":"534.7087","3":"2.410889","4":"0.4743096","5":"5.082944","6":"3.716299e-07","7":"6.329039e-05"},{"1":"ENSG00000182319:SGK223","2":"2574.2730","3":"2.409963","4":"0.3894353","5":"6.188352","6":"6.079650e-10","7":"3.859196e-07"},{"1":"ENSG00000109452:INPP4B","2":"2904.7942","3":"2.402330","4":"0.5127521","5":"4.685169","6":"2.797296e-06","7":"2.714851e-04"},{"1":"ENSG00000138798:EGF","2":"491.2952","3":"2.402068","4":"0.5600660","5":"4.288902","6":"1.795590e-05","7":"1.138067e-03"},{"1":"ENSG00000178999:AURKB","2":"577.4682","3":"2.399356","4":"0.4645804","5":"5.164567","6":"2.409960e-07","7":"4.589330e-05"},{"1":"ENSG00000124693:HIST1H3B","2":"6683.1588","3":"2.395404","4":"0.5810398","5":"4.122616","6":"3.745937e-05","7":"1.956808e-03"},{"1":"ENSG00000143416:SELENBP1","2":"498.0082","3":"2.394126","4":"0.4459566","5":"5.368518","6":"7.938631e-08","7":"1.967982e-05"},{"1":"ENSG00000183856:IQGAP3","2":"1334.6502","3":"2.391487","4":"0.3601194","5":"6.640817","6":"3.119495e-11","7":"3.111696e-08"},{"1":"ENSG00000175906:ARL4D","2":"446.3678","3":"2.387029","4":"0.5164190","5":"4.622271","6":"3.795617e-06","7":"3.383349e-04"},{"1":"ENSG00000101447:FAM83D","2":"560.1051","3":"2.385019","4":"0.4002427","5":"5.958933","6":"2.538903e-09","7":"1.350443e-06"},{"1":"ENSG00000237096:RP4-782G3.1","2":"115.1429","3":"2.379306","4":"0.5926216","5":"4.014882","6":"5.947554e-05","7":"2.780946e-03"},{"1":"ENSG00000198353:HOXC4","2":"111.9040","3":"2.369292","4":"0.6549730","5":"3.617389","6":"2.975898e-04","7":"9.027895e-03"},{"1":"ENSG00000118898:PPL","2":"645.5584","3":"2.368226","4":"0.5420627","5":"4.368916","6":"1.248649e-05","7":"8.646639e-04"},{"1":"ENSG00000172572:PDE3A","2":"2536.2739","3":"2.368066","4":"0.3162728","5":"7.487415","6":"7.024341e-14","7":"1.839280e-10"},{"1":"ENSG00000075218:GTSE1","2":"739.4786","3":"2.367954","4":"0.4164212","5":"5.686440","6":"1.297147e-08","7":"4.725563e-06"},{"1":"ENSG00000124575:HIST1H1D","2":"3199.9219","3":"2.342945","4":"0.5429464","5":"4.315243","6":"1.594275e-05","7":"1.035537e-03"},{"1":"ENSG00000143061:IGSF3","2":"3944.8545","3":"2.342093","4":"0.5518642","5":"4.243966","6":"2.196037e-05","7":"1.332038e-03"},{"1":"ENSG00000106571:GLI3","2":"3038.3629","3":"2.325014","4":"0.4727414","5":"4.918153","6":"8.736471e-07","7":"1.183078e-04"},{"1":"ENSG00000117724:CENPF","2":"12380.0236","3":"2.323928","4":"0.4961347","5":"4.684067","6":"2.812377e-06","7":"2.714851e-04"},{"1":"ENSG00000013810:TACC3","2":"1574.6715","3":"2.319503","4":"0.4252180","5":"5.454855","6":"4.901277e-08","7":"1.324768e-05"},{"1":"ENSG00000224825:RP11-171A24.3","2":"104.3360","3":"2.318627","4":"0.6363791","5":"3.643468","6":"2.689891e-04","7":"8.422127e-03"},{"1":"ENSG00000179603:GRM8","2":"347.3115","3":"2.314213","4":"0.6241698","5":"3.707665","6":"2.091789e-04","7":"7.027705e-03"},{"1":"ENSG00000271860:RP11-436D23.1","2":"208.3132","3":"2.311613","4":"0.5987390","5":"3.860803","6":"1.130151e-04","7":"4.462553e-03"},{"1":"ENSG00000109113:RAB34","2":"1042.8573","3":"2.299290","4":"0.4467341","5":"5.146886","6":"2.648455e-07","7":"4.953118e-05"},{"1":"ENSG00000123485:HJURP","2":"684.3001","3":"2.287870","4":"0.4590694","5":"4.983713","6":"6.237562e-07","7":"9.105319e-05"},{"1":"ENSG00000177602:GSG2","2":"286.6823","3":"2.274176","4":"0.3877142","5":"5.865600","6":"4.475108e-09","7":"1.932831e-06"},{"1":"ENSG00000235363:SNRPGP10","2":"624.9416","3":"2.270428","4":"0.5179555","5":"4.383442","6":"1.168186e-05","7":"8.211814e-04"},{"1":"ENSG00000121621:KIF18A","2":"784.5794","3":"2.255233","4":"0.4551883","5":"4.954507","6":"7.251405e-07","7":"1.022888e-04"},{"1":"ENSG00000101057:MYBL2","2":"956.1312","3":"2.252465","4":"0.4430940","5":"5.083493","6":"3.705564e-07","7":"6.329039e-05"},{"1":"ENSG00000142611:PRDM16","2":"1677.4332","3":"2.245185","4":"0.5661020","5":"3.966043","6":"7.307560e-05","7":"3.253456e-03"},{"1":"ENSG00000075702:WDR62","2":"837.0920","3":"2.242546","4":"0.3733124","5":"6.007156","6":"1.888063e-09","7":"1.054672e-06"},{"1":"ENSG00000270124:RP11-118F19.1","2":"205.3633","3":"2.235852","4":"0.4813858","5":"4.644617","6":"3.407081e-06","7":"3.137136e-04"},{"1":"ENSG00000118495:PLAGL1","2":"740.2544","3":"2.235358","4":"0.5672976","5":"3.940362","6":"8.135888e-05","7":"3.524850e-03"},{"1":"ENSG00000185008:ROBO2","2":"9680.3144","3":"2.230040","4":"0.5990384","5":"3.722700","6":"1.971037e-04","7":"6.757496e-03"},{"1":"ENSG00000105877:DNAH11","2":"2028.2853","3":"2.225012","4":"0.5700196","5":"3.903396","6":"9.485245e-05","7":"3.936206e-03"},{"1":"ENSG00000180340:FZD2","2":"285.4706","3":"2.220342","4":"0.4456423","5":"4.982342","6":"6.281919e-07","7":"9.118888e-05"},{"1":"ENSG00000178814:OPLAH","2":"443.2600","3":"2.214703","4":"0.5485783","5":"4.037168","6":"5.410041e-05","7":"2.587371e-03"},{"1":"ENSG00000237649:KIFC1","2":"1171.6061","3":"2.209918","4":"0.4324104","5":"5.110696","6":"3.209733e-07","7":"5.650075e-05"},{"1":"ENSG00000203811:HIST2H3C","2":"4786.3682","3":"2.191062","4":"0.4984363","5":"4.395873","6":"1.103287e-05","7":"7.887746e-04"},{"1":"ENSG00000186185:KIF18B","2":"1006.7593","3":"2.165951","4":"0.3913607","5":"5.534411","6":"3.122761e-08","7":"9.412093e-06"},{"1":"ENSG00000147082:CCNB3","2":"278.9834","3":"2.164786","4":"0.4409359","5":"4.909524","6":"9.129754e-07","7":"1.216449e-04"},{"1":"ENSG00000197697:HIST1H2BE","2":"1454.6458","3":"2.163603","4":"0.5765034","5":"3.752975","6":"1.747484e-04","7":"6.167391e-03"},{"1":"ENSG00000197496:SLC2A10","2":"1119.8447","3":"2.161328","4":"0.5438495","5":"3.974129","6":"7.063722e-05","7":"3.171861e-03"},{"1":"ENSG00000207405:SNORA64","2":"345.4786","3":"2.158629","4":"0.5761786","5":"3.746457","6":"1.793497e-04","7":"6.303570e-03"},{"1":"ENSG00000142494:SLC47A1","2":"661.7807","3":"2.142567","4":"0.4620526","5":"4.637062","6":"3.533959e-06","7":"3.224018e-04"},{"1":"ENSG00000164185:ZNF474","2":"130.5472","3":"2.133233","4":"0.4677533","5":"4.560595","6":"5.100887e-06","7":"4.291198e-04"},{"1":"ENSG00000224093:RP5-1033H22.2","2":"217.6884","3":"2.125057","4":"0.4281116","5":"4.963792","6":"6.912988e-07","7":"9.784448e-05"},{"1":"ENSG00000066279:ASPM","2":"6341.3693","3":"2.124992","4":"0.5197622","5":"4.088394","6":"4.343705e-05","7":"2.195169e-03"},{"1":"ENSG00000161888:SPC24","2":"725.1165","3":"2.119755","4":"0.4177674","5":"5.074008","6":"3.895229e-07","7":"6.553840e-05"},{"1":"ENSG00000024526:DEPDC1","2":"903.8783","3":"2.119058","4":"0.4731215","5":"4.478888","6":"7.503274e-06","7":"5.864732e-04"},{"1":"ENSG00000214944:ARHGEF28","2":"1075.7588","3":"2.115314","4":"0.5275519","5":"4.009680","6":"6.080109e-05","7":"2.831819e-03"},{"1":"ENSG00000076555:ACACB","2":"1612.5906","3":"2.106603","4":"0.3762334","5":"5.599193","6":"2.153523e-08","7":"7.160463e-06"},{"1":"ENSG00000157399:ARSE","2":"138.1362","3":"2.106285","4":"0.4616819","5":"4.562199","6":"5.062057e-06","7":"4.274138e-04"},{"1":"ENSG00000218422:AC016773.1","2":"346.6943","3":"2.103953","4":"0.4134296","5":"5.089024","6":"3.599116e-07","7":"6.205142e-05"},{"1":"ENSG00000147485:PXDNL","2":"393.6382","3":"2.102594","4":"0.4925774","5":"4.268556","6":"1.967426e-05","7":"1.217509e-03"},{"1":"ENSG00000072840:EVC","2":"2887.9632","3":"2.098898","4":"0.3064199","5":"6.849745","6":"7.398172e-12","7":"7.947344e-09"},{"1":"ENSG00000197409:HIST1H3D","2":"5316.5701","3":"2.094610","4":"0.4431227","5":"4.726931","6":"2.279382e-06","7":"2.363730e-04"},{"1":"ENSG00000090776:EFNB1","2":"1806.5364","3":"2.081736","4":"0.3954568","5":"5.264130","6":"1.408549e-07","7":"3.089588e-05"},{"1":"ENSG00000140009:ESR2","2":"618.8283","3":"2.079108","4":"0.5686048","5":"3.656507","6":"2.556754e-04","7":"8.101349e-03"},{"1":"ENSG00000143195:ILDR2","2":"11864.8644","3":"2.078850","4":"0.5145871","5":"4.039841","6":"5.348746e-05","7":"2.563910e-03"},{"1":"ENSG00000214717:ZBED1","2":"250.0744","3":"2.078359","4":"0.5285287","5":"3.932348","6":"8.412031e-05","7":"3.610882e-03"},{"1":"ENSG00000267396:RP11-845C23.3","2":"377.3202","3":"2.072581","4":"0.5148924","5":"4.025270","6":"5.690986e-05","7":"2.687981e-03"},{"1":"ENSG00000166803:KIAA0101","2":"373.8873","3":"2.071067","4":"0.4563221","5":"4.538608","6":"5.662692e-06","7":"4.697792e-04"},{"1":"ENSG00000157456:CCNB2","2":"959.6384","3":"2.061142","4":"0.4395012","5":"4.689731","6":"2.735648e-06","7":"2.684075e-04"},{"1":"ENSG00000049759:NEDD4L","2":"6173.1219","3":"2.060064","4":"0.4448486","5":"4.630932","6":"3.640231e-06","7":"3.308189e-04"},{"1":"ENSG00000267491:CTD-2373H9.5","2":"151.8970","3":"2.050915","4":"0.5648090","5":"3.631165","6":"2.821445e-04","7":"8.646997e-03"},{"1":"ENSG00000199785:SNORA52","2":"231.3372","3":"2.050557","4":"0.5242611","5":"3.911329","6":"9.178971e-05","7":"3.857101e-03"},{"1":"ENSG00000173894:CBX2","2":"1272.7939","3":"2.049145","4":"0.4183155","5":"4.898564","6":"9.653942e-07","7":"1.260741e-04"},{"1":"ENSG00000134247:PTGFRN","2":"5114.7330","3":"2.048236","4":"0.3819060","5":"5.363193","6":"8.176350e-08","7":"2.014989e-05"},{"1":"ENSG00000007372:PAX6","2":"3216.4001","3":"2.043788","4":"0.3847829","5":"5.311535","6":"1.087056e-07","7":"2.547732e-05"},{"1":"ENSG00000148773:MKI67","2":"8637.3103","3":"2.041775","4":"0.4676627","5":"4.365914","6":"1.265924e-05","7":"8.712383e-04"},{"1":"ENSG00000137259:HIST1H2AB","2":"794.1360","3":"2.040968","4":"0.5169845","5":"3.947833","6":"7.886167e-05","7":"3.434417e-03"},{"1":"ENSG00000109220:CHIC2","2":"1536.1712","3":"2.039253","4":"0.3845247","5":"5.303308","6":"1.137227e-07","7":"2.646896e-05"},{"1":"ENSG00000100427:MLC1","2":"6813.3026","3":"2.030066","4":"0.5132725","5":"3.955143","6":"7.648897e-05","7":"3.373163e-03"},{"1":"ENSG00000116194:ANGPTL1","2":"1179.1336","3":"2.027712","4":"0.5565612","5":"3.643287","6":"2.691784e-04","7":"8.422127e-03"},{"1":"ENSG00000142945:KIF2C","2":"1205.6714","3":"2.019816","4":"0.4311040","5":"4.685218","6":"2.796626e-06","7":"2.714851e-04"},{"1":"ENSG00000175063:UBE2C","2":"1117.0161","3":"2.017552","4":"0.5104957","5":"3.952143","6":"7.745458e-05","7":"3.412155e-03"},{"1":"ENSG00000129195:FAM64A","2":"782.8360","3":"2.015519","4":"0.4892564","5":"4.119557","6":"3.796017e-05","7":"1.978037e-03"},{"1":"ENSG00000121152:NCAPH","2":"759.1061","3":"2.014728","4":"0.4156879","5":"4.846733","6":"1.255114e-06","7":"1.495444e-04"},{"1":"ENSG00000143476:DTL","2":"1132.3556","3":"2.014417","4":"0.4206624","5":"4.788678","6":"1.678837e-06","7":"1.850917e-04"},{"1":"ENSG00000140451:PIF1","2":"326.1837","3":"2.014201","4":"0.4487841","5":"4.488130","6":"7.185107e-06","7":"5.668928e-04"},{"1":"ENSG00000270408:JAG1","2":"902.2993","3":"2.006882","4":"0.4210900","5":"4.765922","6":"1.879916e-06","7":"2.014298e-04"},{"1":"ENSG00000131747:TOP2A","2":"11638.6622","3":"2.005287","4":"0.5200997","5":"3.855583","6":"1.154543e-04","7":"4.533232e-03"},{"1":"ENSG00000101384:JAG1","2":"13601.9452","3":"2.002077","4":"0.4009741","5":"4.993033","6":"5.943842e-07","7":"8.768214e-05"},{"1":"ENSG00000239268:RP11-384F7.2","2":"321.2299","3":"-2.002005","4":"0.5426356","5":"-3.689408","6":"2.247761e-04","7":"7.357027e-03"},{"1":"ENSG00000107742:SPOCK2","2":"4703.6774","3":"-2.007148","4":"0.5314826","5":"-3.776508","6":"1.590423e-04","7":"5.788945e-03"},{"1":"ENSG00000260249:RP11-401P9.5","2":"104.5389","3":"-2.018564","4":"0.4160113","5":"-4.852185","6":"1.221086e-06","7":"1.470710e-04"},{"1":"ENSG00000110090:CPT1A","2":"1140.7729","3":"-2.046873","4":"0.4005745","5":"-5.109845","6":"3.224228e-07","7":"5.651843e-05"},{"1":"ENSG00000142549:IGLON5","2":"131.5586","3":"-2.079034","4":"0.4622853","5":"-4.497297","6":"6.882292e-06","7":"5.481628e-04"},{"1":"ENSG00000228214:LINC00693","2":"142.4468","3":"-2.087727","4":"0.4994492","5":"-4.180059","6":"2.914338e-05","7":"1.654420e-03"},{"1":"ENSG00000143473:KCNH1","2":"358.4190","3":"-2.099403","4":"0.5667980","5":"-3.703970","6":"2.122515e-04","7":"7.102457e-03"},{"1":"ENSG00000077063:CTTNBP2","2":"2608.0500","3":"-2.104605","4":"0.5280722","5":"-3.985449","6":"6.735266e-05","7":"3.066729e-03"},{"1":"ENSG00000113231:PDE8B","2":"3370.7352","3":"-2.119170","4":"0.4729245","5":"-4.480990","6":"7.429776e-06","7":"5.829035e-04"},{"1":"ENSG00000228213:NLGN1-AS1","2":"338.5690","3":"-2.127358","4":"0.5117034","5":"-4.157404","6":"3.218841e-05","7":"1.761035e-03"},{"1":"ENSG00000134042:MRO","2":"263.7161","3":"-2.127464","4":"0.5720980","5":"-3.718706","6":"2.002459e-04","7":"6.826120e-03"},{"1":"ENSG00000253553:RP11-586K2.1","2":"282.3713","3":"-2.134119","4":"0.5921286","5":"-3.604147","6":"3.131799e-04","7":"9.358538e-03"},{"1":"ENSG00000182247:UBE2E2","2":"768.5224","3":"-2.166608","4":"0.4488326","5":"-4.827207","6":"1.384615e-06","7":"1.602443e-04"},{"1":"ENSG00000152402:GUCY1A2","2":"5040.8641","3":"-2.167224","4":"0.5183114","5":"-4.181316","6":"2.898268e-05","7":"1.647530e-03"},{"1":"ENSG00000163053:SLC16A14","2":"563.3286","3":"-2.167714","4":"0.3896597","5":"-5.563096","6":"2.650307e-08","7":"8.448246e-06"},{"1":"ENSG00000205682:RP11-798L4.1","2":"154.0503","3":"-2.177579","4":"0.5748672","5":"-3.787970","6":"1.518834e-04","7":"5.581717e-03"},{"1":"ENSG00000145808:ADAMTS19","2":"102.1412","3":"-2.194216","4":"0.5954648","5":"-3.684878","6":"2.288118e-04","7":"7.431062e-03"},{"1":"ENSG00000269891:ARHGAP19-SLIT1","2":"2598.7579","3":"-2.196315","4":"0.5304835","5":"-4.140214","6":"3.469825e-05","7":"1.861310e-03"},{"1":"ENSG00000167680:SEMA6B","2":"663.3282","3":"-2.213553","4":"0.4132754","5":"-5.356122","6":"8.502723e-08","7":"2.071056e-05"},{"1":"ENSG00000243197:TUSC7","2":"447.2388","3":"-2.213850","4":"0.6086306","5":"-3.637427","6":"2.753750e-04","7":"8.531117e-03"},{"1":"ENSG00000072163:LIMS2","2":"334.1529","3":"-2.218899","4":"0.5328479","5":"-4.164226","6":"3.124103e-05","7":"1.740483e-03"},{"1":"ENSG00000262223:RP11-1055B8.3","2":"122.3041","3":"-2.222318","4":"0.5865746","5":"-3.788636","6":"1.514764e-04","7":"5.576540e-03"},{"1":"ENSG00000090659:CD209","2":"136.9067","3":"-2.228125","4":"0.6215569","5":"-3.584748","6":"3.374035e-04","7":"9.885670e-03"},{"1":"ENSG00000135905:DOCK10","2":"5627.7742","3":"-2.239887","4":"0.5579795","5":"-4.014282","6":"5.962710e-05","7":"2.781824e-03"},{"1":"ENSG00000204136:GGTA1P","2":"360.3983","3":"-2.253815","4":"0.5380319","5":"-4.188999","6":"2.801873e-05","7":"1.610212e-03"},{"1":"ENSG00000111962:UST","2":"2877.1464","3":"-2.254353","4":"0.4539230","5":"-4.966378","6":"6.821495e-07","7":"9.687679e-05"},{"1":"ENSG00000182118:FAM89A","2":"159.4671","3":"-2.260478","4":"0.4199588","5":"-5.382619","6":"7.340978e-08","7":"1.852712e-05"},{"1":"ENSG00000006016:CRLF1","2":"308.5955","3":"-2.283814","4":"0.5268421","5":"-4.334912","6":"1.458182e-05","7":"9.666227e-04"},{"1":"ENSG00000224668:IPO8P1","2":"156.3992","3":"-2.286803","4":"0.4148303","5":"-5.512624","6":"3.535234e-08","7":"1.035725e-05"},{"1":"ENSG00000163618:CADPS","2":"4159.2833","3":"-2.290772","4":"0.6176412","5":"-3.708905","6":"2.081578e-04","7":"7.015905e-03"},{"1":"ENSG00000255701:AC090953.1","2":"554.0609","3":"-2.291764","4":"0.5453847","5":"-4.202105","6":"2.644448e-05","7":"1.544866e-03"},{"1":"ENSG00000159263:SIM2","2":"331.1521","3":"-2.291868","4":"0.6205013","5":"-3.693575","6":"2.211229e-04","7":"7.294445e-03"},{"1":"ENSG00000149403:GRIK4","2":"638.4344","3":"-2.298124","4":"0.5324123","5":"-4.316436","6":"1.585684e-05","7":"1.034620e-03"},{"1":"ENSG00000223764:RP11-54O7.3","2":"110.3486","3":"-2.299326","4":"0.4549544","5":"-5.053971","6":"4.327165e-07","7":"7.081507e-05"},{"1":"ENSG00000173210:ABLIM3","2":"1003.7289","3":"-2.299375","4":"0.4952212","5":"-4.643126","6":"3.431769e-06","7":"3.152938e-04"},{"1":"ENSG00000204176:SYT15","2":"215.6487","3":"-2.300154","4":"0.4362839","5":"-5.272149","6":"1.348354e-07","7":"2.973120e-05"},{"1":"ENSG00000111728:ST8SIA1","2":"1418.4662","3":"-2.314427","4":"0.4724819","5":"-4.898445","6":"9.659811e-07","7":"1.260741e-04"},{"1":"ENSG00000131459:GFPT2","2":"1474.6724","3":"-2.321219","4":"0.5616391","5":"-4.132938","6":"3.581559e-05","7":"1.904180e-03"},{"1":"ENSG00000186642:PDE2A","2":"898.3078","3":"-2.350245","4":"0.5951775","5":"-3.948813","6":"7.853969e-05","7":"3.431095e-03"},{"1":"ENSG00000227660:RP11-162J8.2","2":"163.0496","3":"-2.351720","4":"0.4836557","5":"-4.862383","6":"1.159807e-06","7":"1.424930e-04"},{"1":"ENSG00000116675:DNAJC6","2":"2158.1809","3":"-2.352600","4":"0.5543649","5":"-4.243775","6":"2.197910e-05","7":"1.332038e-03"},{"1":"ENSG00000157890:MEGF11","2":"1221.9709","3":"-2.363347","4":"0.3713128","5":"-6.364842","6":"1.954908e-10","7":"1.462516e-07"},{"1":"ENSG00000242173:ARHGDIG","2":"246.3597","3":"-2.383980","4":"0.6084273","5":"-3.918266","6":"8.918843e-05","7":"3.785764e-03"},{"1":"ENSG00000269975:RP11-58B17.2","2":"206.9985","3":"-2.385142","4":"0.5195340","5":"-4.590926","6":"4.412841e-06","7":"3.867698e-04"},{"1":"ENSG00000272143:FGF14-AS2","2":"141.0177","3":"-2.397162","4":"0.4673912","5":"-5.128812","6":"2.915765e-07","7":"5.265344e-05"},{"1":"ENSG00000147724:FAM135B","2":"886.0731","3":"-2.401817","4":"0.5991824","5":"-4.008490","6":"6.110817e-05","7":"2.838278e-03"},{"1":"ENSG00000207955:MIR219-2","2":"179.8909","3":"-2.416961","4":"0.4614297","5":"-5.237983","6":"1.623414e-07","7":"3.383727e-05"},{"1":"ENSG00000174348:PODN","2":"370.2636","3":"-2.417685","4":"0.4351914","5":"-5.555451","6":"2.768952e-08","7":"8.620745e-06"},{"1":"ENSG00000055118:KCNH2","2":"597.6804","3":"-2.420757","4":"0.6566355","5":"-3.686607","6":"2.272637e-04","7":"7.405281e-03"},{"1":"ENSG00000151090:THRB","2":"536.4620","3":"-2.438472","4":"0.6577440","5":"-3.707327","6":"2.094581e-04","7":"7.031447e-03"},{"1":"ENSG00000158445:KCNB1","2":"1001.8979","3":"-2.440834","4":"0.6575817","5":"-3.711834","6":"2.057631e-04","7":"6.957584e-03"},{"1":"ENSG00000125845:BMP2","2":"870.7819","3":"-2.446861","4":"0.5725220","5":"-4.273828","6":"1.921450e-05","7":"1.192580e-03"},{"1":"ENSG00000168874:ATOH8","2":"1007.9344","3":"-2.453530","4":"0.4862556","5":"-5.045761","6":"4.517207e-07","7":"7.306888e-05"},{"1":"ENSG00000186231:KLHL32","2":"459.1889","3":"-2.454302","4":"0.6848724","5":"-3.583591","6":"3.389030e-04","7":"9.908124e-03"},{"1":"ENSG00000241261:RPL17P19","2":"184.6122","3":"-2.456294","4":"0.5568026","5":"-4.411426","6":"1.026919e-05","7":"7.508338e-04"},{"1":"ENSG00000197134:ZNF257","2":"124.4993","3":"-2.478556","4":"0.6025351","5":"-4.113546","6":"3.896269e-05","7":"2.017728e-03"},{"1":"ENSG00000171044:XKR6","2":"237.8072","3":"-2.486525","4":"0.6356276","5":"-3.911921","6":"9.156491e-05","7":"3.856575e-03"},{"1":"ENSG00000064989:CALCRL","2":"4378.1276","3":"-2.497571","4":"0.5306757","5":"-4.706398","6":"2.521318e-06","7":"2.539198e-04"},{"1":"ENSG00000172995:ARPP21","2":"914.9728","3":"-2.499863","4":"0.6581098","5":"-3.798549","6":"1.455454e-04","7":"5.400909e-03"},{"1":"ENSG00000215346:AF131215.5","2":"188.6992","3":"-2.502182","4":"0.6553450","5":"-3.818114","6":"1.344759e-04","7":"5.084717e-03"},{"1":"ENSG00000254556:AF131215.4","2":"104.1756","3":"-2.510030","4":"0.6996340","5":"-3.587632","6":"3.336943e-04","7":"9.808914e-03"},{"1":"ENSG00000156103:MMP16","2":"7087.2680","3":"-2.511874","4":"0.6091476","5":"-4.123588","6":"3.730157e-05","7":"1.952376e-03"},{"1":"ENSG00000128285:MCHR1","2":"215.7845","3":"-2.536791","4":"0.6601354","5":"-3.842835","6":"1.216215e-04","7":"4.719033e-03"},{"1":"ENSG00000165495:PKNOX2","2":"576.7823","3":"-2.555153","4":"0.3677487","5":"-6.948096","6":"3.702490e-12","7":"4.562231e-09"},{"1":"ENSG00000161681:SHANK1","2":"114.2295","3":"-2.556236","4":"0.6888323","5":"-3.710970","6":"2.064667e-04","7":"6.970124e-03"},{"1":"ENSG00000197959:DNM3","2":"4132.7729","3":"-2.563224","4":"0.4593002","5":"-5.580717","6":"2.395296e-08","7":"7.839915e-06"},{"1":"ENSG00000169783:LINGO1","2":"612.5171","3":"-2.568782","4":"0.6347840","5":"-4.046703","6":"5.194413e-05","7":"2.510034e-03"},{"1":"ENSG00000171246:NPTX1","2":"731.4768","3":"-2.575098","4":"0.6476487","5":"-3.976072","6":"7.006288e-05","7":"3.154608e-03"},{"1":"ENSG00000117318:ID3","2":"2620.9090","3":"-2.581006","4":"0.4991880","5":"-5.170408","6":"2.335836e-07","7":"4.468485e-05"},{"1":"ENSG00000125968:ID1","2":"1183.1748","3":"-2.591917","4":"0.6049001","5":"-4.284868","6":"1.828474e-05","7":"1.150209e-03"},{"1":"ENSG00000254545:RP11-84A19.3","2":"121.9254","3":"-2.601303","4":"0.6726323","5":"-3.867348","6":"1.100253e-04","7":"4.360935e-03"},{"1":"ENSG00000116544:DLGAP3","2":"128.9873","3":"-2.604280","4":"0.6677496","5":"-3.900085","6":"9.615881e-05","7":"3.961232e-03"},{"1":"ENSG00000173705:SUSD5","2":"904.0435","3":"-2.621503","4":"0.7197789","5":"-3.642095","6":"2.704282e-04","7":"8.454917e-03"},{"1":"ENSG00000117154:IGSF21","2":"335.9436","3":"-2.635867","4":"0.7157939","5":"-3.682439","6":"2.310136e-04","7":"7.479378e-03"},{"1":"ENSG00000131773:KHDRBS3","2":"1616.2595","3":"-2.639585","4":"0.6084860","5":"-4.337956","6":"1.438142e-05","7":"9.594101e-04"},{"1":"ENSG00000041353:RAB27B","2":"463.1501","3":"-2.643890","4":"0.6541554","5":"-4.041684","6":"5.306862e-05","7":"2.549667e-03"},{"1":"ENSG00000187672:ERC2","2":"620.6908","3":"-2.666702","4":"0.6068292","5":"-4.394485","6":"1.110359e-05","7":"7.911305e-04"},{"1":"ENSG00000152583:SPARCL1","2":"49371.4794","3":"-2.668528","4":"0.6123674","5":"-4.357723","6":"1.314227e-05","7":"8.967352e-04"},{"1":"ENSG00000101463:SYNDIG1","2":"463.6472","3":"-2.679046","4":"0.6291339","5":"-4.258308","6":"2.059801e-05","7":"1.267565e-03"},{"1":"ENSG00000227502:RP1-249H1.4","2":"496.5241","3":"-2.711129","4":"0.5178941","5":"-5.234910","6":"1.650651e-07","7":"3.423467e-05"},{"1":"ENSG00000174938:SEZ6L2","2":"620.8363","3":"-2.713944","4":"0.7230086","5":"-3.753681","6":"1.742566e-04","7":"6.165948e-03"},{"1":"ENSG00000001626:CFTR","2":"132.2769","3":"-2.715223","4":"0.6280869","5":"-4.323006","6":"1.539178e-05","7":"1.012306e-03"},{"1":"ENSG00000133958:UNC79","2":"891.0644","3":"-2.730413","4":"0.6915984","5":"-3.947975","6":"7.881500e-05","7":"3.434417e-03"},{"1":"ENSG00000146147:MLIP","2":"210.4398","3":"-2.739113","4":"0.6935638","5":"-3.949332","6":"7.836972e-05","7":"3.427244e-03"},{"1":"ENSG00000157680:DGKI","2":"1946.8562","3":"-2.739433","4":"0.5575816","5":"-4.913063","6":"8.966456e-07","7":"1.204005e-04"},{"1":"ENSG00000143632:ACTA1","2":"126.1544","3":"-2.743567","4":"0.6069417","5":"-4.520315","6":"6.174780e-06","7":"5.013418e-04"},{"1":"ENSG00000102452:NALCN","2":"1342.8871","3":"-2.748542","4":"0.6918406","5":"-3.972797","6":"7.103360e-05","7":"3.181460e-03"},{"1":"ENSG00000185666:SYN3","2":"161.1279","3":"-2.753651","4":"0.6196669","5":"-4.443760","6":"8.840037e-06","7":"6.661032e-04"},{"1":"ENSG00000245648:RP11-277P12.20","2":"199.4006","3":"-2.776714","4":"0.5347100","5":"-5.192935","6":"2.070047e-07","7":"4.101647e-05"},{"1":"ENSG00000107295:SH3GL2","2":"918.4078","3":"-2.792773","4":"0.6754617","5":"-4.134614","6":"3.555519e-05","7":"1.892738e-03"},{"1":"ENSG00000269918:AF131215.9","2":"186.1558","3":"-2.798206","4":"0.6690473","5":"-4.182374","6":"2.884812e-05","7":"1.644343e-03"},{"1":"ENSG00000171798:KNDC1","2":"340.0299","3":"-2.804108","4":"0.7625325","5":"-3.677363","6":"2.356580e-04","7":"7.600377e-03"},{"1":"ENSG00000054356:PTPRN","2":"547.6781","3":"-2.817514","4":"0.7507531","5":"-3.752917","6":"1.747887e-04","7":"6.167391e-03"},{"1":"ENSG00000077080:ACTL6B","2":"120.2058","3":"-2.833143","4":"0.7451508","5":"-3.802106","6":"1.434711e-04","7":"5.338118e-03"},{"1":"ENSG00000102466:FGF14","2":"484.4211","3":"-2.840315","4":"0.7773177","5":"-3.653995","6":"2.581914e-04","7":"8.169886e-03"},{"1":"ENSG00000132561:MATN2","2":"2610.9001","3":"-2.841428","4":"0.5216082","5":"-5.447438","6":"5.110049e-08","7":"1.363602e-05"},{"1":"ENSG00000124194:GDAP1L1","2":"202.9103","3":"-2.849977","4":"0.5649926","5":"-5.044274","6":"4.552474e-07","7":"7.335612e-05"},{"1":"ENSG00000174521:TTC9B","2":"102.7447","3":"-2.860662","4":"0.5630359","5":"-5.080781","6":"3.758865e-07","7":"6.375614e-05"},{"1":"ENSG00000187527:ATP13A5","2":"177.7859","3":"-2.861829","4":"0.5587753","5":"-5.121611","6":"3.029371e-07","7":"5.423739e-05"},{"1":"ENSG00000255310:AF131215.2","2":"278.6602","3":"-2.871609","4":"0.6689191","5":"-4.292910","6":"1.763468e-05","7":"1.122804e-03"},{"1":"ENSG00000155265:GOLGA7B","2":"270.4963","3":"-2.874450","4":"0.6958063","5":"-4.131107","6":"3.610198e-05","7":"1.909713e-03"},{"1":"ENSG00000206579:XKR4","2":"881.8107","3":"-2.875204","4":"0.7706133","5":"-3.731059","6":"1.906766e-04","7":"6.606756e-03"},{"1":"ENSG00000141837:CACNA1A","2":"858.5855","3":"-2.893767","4":"0.6496908","5":"-4.454068","6":"8.425835e-06","7":"6.441612e-04"},{"1":"ENSG00000226869:LHFPL3-AS1","2":"238.0479","3":"-2.910678","4":"0.7596822","5":"-3.831441","6":"1.273948e-04","7":"4.892032e-03"},{"1":"ENSG00000144406:UNC80","2":"3641.1625","3":"-2.914675","4":"0.6699270","5":"-4.350735","6":"1.356819e-05","7":"9.198045e-04"},{"1":"ENSG00000170017:ALCAM","2":"5169.8838","3":"-2.925467","4":"0.5635039","5":"-5.191565","6":"2.085334e-07","7":"4.101647e-05"},{"1":"ENSG00000116983:HPCAL4","2":"438.8479","3":"-2.978138","4":"0.7963998","5":"-3.739501","6":"1.843857e-04","7":"6.432006e-03"},{"1":"ENSG00000100433:KCNK10","2":"206.4647","3":"-2.978419","4":"0.7936400","5":"-3.752859","6":"1.748293e-04","7":"6.167391e-03"},{"1":"ENSG00000134594:RAB33A","2":"121.7694","3":"-2.990533","4":"0.4792868","5":"-6.239548","6":"4.388384e-10","7":"2.918275e-07"},{"1":"ENSG00000131771:PPP1R1B","2":"429.3309","3":"-3.010443","4":"0.6821053","5":"-4.413457","6":"1.017329e-05","7":"7.464274e-04"},{"1":"ENSG00000148408:CACNA1B","2":"246.5342","3":"-3.012500","4":"0.8065386","5":"-3.735097","6":"1.876431e-04","7":"6.525236e-03"},{"1":"ENSG00000154118:JPH3","2":"442.0551","3":"-3.028357","4":"0.5965771","5":"-5.076221","6":"3.850152e-07","7":"6.504118e-05"},{"1":"ENSG00000148053:NTRK2","2":"16974.0458","3":"-3.030012","4":"0.6638594","5":"-4.564237","6":"5.013130e-06","7":"4.251520e-04"},{"1":"ENSG00000242441:GTF2A1L","2":"202.2602","3":"-3.047908","4":"0.6199705","5":"-4.916214","6":"8.823382e-07","7":"1.188603e-04"},{"1":"ENSG00000196581:AJAP1","2":"668.5800","3":"-3.048376","4":"0.7574846","5":"-4.024340","6":"5.713530e-05","7":"2.693123e-03"},{"1":"ENSG00000249307:LINC01088","2":"448.0789","3":"-3.054286","4":"0.7337318","5":"-4.162674","6":"3.145417e-05","7":"1.743566e-03"},{"1":"ENSG00000214870:AC004540.5","2":"152.2819","3":"-3.062440","4":"0.7288341","5":"-4.201835","6":"2.647605e-05","7":"1.544866e-03"},{"1":"ENSG00000213809:KLRK1","2":"238.2692","3":"-3.072058","4":"0.4962505","5":"-6.190538","6":"5.995907e-10","7":"3.859196e-07"},{"1":"ENSG00000260990:RP3-518E13.2","2":"160.7752","3":"-3.078324","4":"0.7873820","5":"-3.909569","6":"9.246094e-05","7":"3.877298e-03"},{"1":"ENSG00000116147:TNR","2":"8309.7342","3":"-3.082269","4":"0.7869927","5":"-3.916515","6":"8.983826e-05","7":"3.805636e-03"},{"1":"ENSG00000176406:RIMS2","2":"1071.9619","3":"-3.098535","4":"0.8006074","5":"-3.870230","6":"1.087327e-04","7":"4.330187e-03"},{"1":"ENSG00000137731:FXYD2","2":"189.3553","3":"-3.104217","4":"0.6875652","5":"-4.514796","6":"6.337780e-06","7":"5.125894e-04"},{"1":"ENSG00000152578:GRIA4","2":"5199.4899","3":"-3.139367","4":"0.7701087","5":"-4.076524","6":"4.571391e-05","7":"2.285423e-03"},{"1":"ENSG00000248587:GDNF-AS1","2":"174.3623","3":"-3.157592","4":"0.6518002","5":"-4.844417","6":"1.269838e-06","7":"1.504667e-04"},{"1":"ENSG00000136531:SCN2A","2":"1411.4307","3":"-3.173059","4":"0.6640491","5":"-4.778349","6":"1.767401e-06","7":"1.933297e-04"},{"1":"ENSG00000227757:AP000282.2","2":"110.6457","3":"-3.225460","4":"0.7315397","5":"-4.409139","6":"1.037825e-05","7":"7.548556e-04"},{"1":"ENSG00000196542:SPTSSB","2":"219.6625","3":"-3.236583","4":"0.6163081","5":"-5.251566","6":"1.508111e-07","7":"3.256820e-05"},{"1":"ENSG00000136928:GABBR2","2":"839.7770","3":"-3.241039","4":"0.6919061","5":"-4.684218","6":"2.810309e-06","7":"2.714851e-04"},{"1":"ENSG00000166897:ELFN2","2":"548.1635","3":"-3.243355","4":"0.8025518","5":"-4.041303","6":"5.315492e-05","7":"2.550888e-03"},{"1":"ENSG00000101542:CDH20","2":"831.4483","3":"-3.245234","4":"0.7106947","5":"-4.566284","6":"4.964465e-06","7":"4.241987e-04"},{"1":"ENSG00000116254:CHD5","2":"413.7435","3":"-3.246795","4":"0.8508340","5":"-3.816014","6":"1.356248e-04","7":"5.118919e-03"},{"1":"ENSG00000120729:MYOT","2":"138.5330","3":"-3.264262","4":"0.6498995","5":"-5.022718","6":"5.094544e-07","7":"7.792150e-05"},{"1":"ENSG00000237289:CKMT1B","2":"126.2853","3":"-3.286741","4":"0.7838088","5":"-4.193295","6":"2.749315e-05","7":"1.588725e-03"},{"1":"ENSG00000165388:ZNF488","2":"535.7848","3":"-3.303469","4":"0.8718830","5":"-3.788890","6":"1.513217e-04","7":"5.575747e-03"},{"1":"ENSG00000143502:SUSD4","2":"459.6828","3":"-3.330238","4":"0.6315730","5":"-5.272927","6":"1.342649e-07","7":"2.973120e-05"},{"1":"ENSG00000132329:RAMP1","2":"501.8569","3":"-3.338860","4":"0.8021070","5":"-4.162611","6":"3.146285e-05","7":"1.743566e-03"},{"1":"ENSG00000092051:JPH4","2":"684.0404","3":"-3.339482","4":"0.6253136","5":"-5.340491","6":"9.269506e-08","7":"2.206511e-05"},{"1":"ENSG00000245832:RP11-179A16.1","2":"297.8936","3":"-3.403065","4":"0.7598922","5":"-4.478353","6":"7.522127e-06","7":"5.868520e-04"},{"1":"ENSG00000204442:FAM155A","2":"382.7206","3":"-3.426716","4":"0.7579985","5":"-4.520743","6":"6.162308e-06","7":"5.013007e-04"},{"1":"ENSG00000154654:NCAM2","2":"3508.4355","3":"-3.432532","4":"0.6815835","5":"-5.036115","6":"4.750758e-07","7":"7.488647e-05"},{"1":"ENSG00000263371:AP000673.1","2":"110.1923","3":"-3.460184","4":"0.8989464","5":"-3.849155","6":"1.185261e-04","7":"4.629019e-03"},{"1":"ENSG00000002746:HECW1","2":"304.6331","3":"-3.461973","4":"0.8273125","5":"-4.184602","6":"2.856662e-05","7":"1.632796e-03"},{"1":"ENSG00000161082:CELF5","2":"303.7389","3":"-3.468807","4":"0.6898597","5":"-5.028279","6":"4.949013e-07","7":"7.679218e-05"},{"1":"ENSG00000168959:GRM5","2":"505.0082","3":"-3.495385","4":"0.9301318","5":"-3.757946","6":"1.713136e-04","7":"6.092687e-03"},{"1":"ENSG00000176204:LRRTM4","2":"1110.5983","3":"-3.504747","4":"0.8588214","5":"-4.080880","6":"4.486542e-05","7":"2.251062e-03"},{"1":"ENSG00000182870:GALNT9","2":"136.4146","3":"-3.511724","4":"0.7183904","5":"-4.888322","6":"1.016991e-06","7":"1.306921e-04"},{"1":"ENSG00000145681:HAPLN1","2":"434.0725","3":"-3.520961","4":"0.8576578","5":"-4.105321","6":"4.037530e-05","7":"2.067877e-03"},{"1":"ENSG00000173404:INSM1","2":"102.4034","3":"-3.523776","4":"0.9623123","5":"-3.661781","6":"2.504684e-04","7":"7.973688e-03"},{"1":"ENSG00000153820:SPHKAP","2":"301.2180","3":"-3.544131","4":"0.8667206","5":"-4.089128","6":"4.329984e-05","7":"2.190878e-03"},{"1":"ENSG00000186369:LINC00643","2":"102.4955","3":"-3.548255","4":"0.7584821","5":"-4.678099","6":"2.895464e-06","7":"2.744468e-04"},{"1":"ENSG00000260903:XKR7","2":"102.3749","3":"-3.549876","4":"0.8235583","5":"-4.310413","6":"1.629500e-05","7":"1.051894e-03"},{"1":"ENSG00000100181:TPTEP1","2":"552.4446","3":"-3.551519","4":"0.7254689","5":"-4.895480","6":"9.806588e-07","7":"1.271972e-04"},{"1":"ENSG00000255819:KLRC4-KLRK1","2":"375.0960","3":"-3.559474","4":"0.5450476","5":"-6.530576","6":"6.551727e-11","7":"5.718429e-08"},{"1":"ENSG00000095713:CRTAC1","2":"372.8592","3":"-3.591630","4":"0.6380819","5":"-5.628792","6":"1.814766e-08","7":"6.389042e-06"},{"1":"ENSG00000162728:KCNJ9","2":"714.1787","3":"-3.596284","4":"0.5560455","5":"-6.467606","6":"9.956740e-11","7":"8.404794e-08"},{"1":"ENSG00000170091:NSG2","2":"797.4938","3":"-3.598025","4":"0.9267246","5":"-3.882519","6":"1.033801e-04","7":"4.164199e-03"},{"1":"ENSG00000148798:INA","2":"605.1599","3":"-3.621242","4":"0.8141428","5":"-4.447920","6":"8.670590e-06","7":"6.559706e-04"},{"1":"ENSG00000187902:SHISA7","2":"299.4193","3":"-3.651465","4":"0.6499074","5":"-5.618439","6":"1.926907e-08","7":"6.522217e-06"},{"1":"ENSG00000243319:FGF14-IT1","2":"272.5157","3":"-3.654343","4":"0.7998671","5":"-4.568688","6":"4.907875e-06","7":"4.222083e-04"},{"1":"ENSG00000272781:MDGA2","2":"625.1012","3":"-3.661496","4":"0.6847589","5":"-5.347131","6":"8.935910e-08","7":"2.139257e-05"},{"1":"ENSG00000154914:USP43","2":"114.4435","3":"-3.667432","4":"0.6169835","5":"-5.944132","6":"2.779254e-09","7":"1.402854e-06"},{"1":"ENSG00000139915:MDGA2","2":"997.3626","3":"-3.677501","4":"0.6863565","5":"-5.358004","6":"8.414632e-08","7":"2.061585e-05"},{"1":"ENSG00000144278:GALNT13","2":"2928.2088","3":"-3.679934","4":"0.7951921","5":"-4.627730","6":"3.696957e-06","7":"3.330839e-04"},{"1":"ENSG00000144230:GPR17","2":"338.3136","3":"-3.704356","4":"0.8420853","5":"-4.399027","6":"1.087374e-05","7":"7.800602e-04"},{"1":"ENSG00000018236:CNTN1","2":"6475.5491","3":"-3.710244","4":"0.8303696","5":"-4.468184","6":"7.888627e-06","7":"6.131615e-04"},{"1":"ENSG00000203805:PPAPDC1A","2":"198.3986","3":"-3.713692","4":"0.8120728","5":"-4.573103","6":"4.805533e-06","7":"4.151089e-04"},{"1":"ENSG00000255082:GRM5-AS1","2":"241.5898","3":"-3.764424","4":"0.9515178","5":"-3.956231","6":"7.614151e-05","7":"3.361379e-03"},{"1":"ENSG00000167654:ATCAY","2":"1538.9353","3":"-3.772948","4":"0.7136297","5":"-5.286984","6":"1.243497e-07","7":"2.812489e-05"},{"1":"ENSG00000184408:KCND2","2":"1305.0875","3":"-3.835282","4":"0.9293865","5":"-4.126681","6":"3.680360e-05","7":"1.939480e-03"},{"1":"ENSG00000243915:PRKRIRP2","2":"136.0819","3":"-3.837596","4":"0.6497391","5":"-5.906364","6":"3.497417e-09","7":"1.636386e-06"},{"1":"ENSG00000150051:MKX","2":"109.5541","3":"-3.900304","4":"0.7071680","5":"-5.515385","6":"3.480163e-08","7":"1.034053e-05"},{"1":"ENSG00000187122:SLIT1","2":"1981.9035","3":"-3.905272","4":"0.7616159","5":"-5.127614","6":"2.934378e-07","7":"5.276214e-05"},{"1":"ENSG00000155897:ADCY8","2":"641.9610","3":"-3.919321","4":"0.7137886","5":"-5.490871","6":"3.999565e-08","7":"1.147683e-05"},{"1":"ENSG00000188730:VWC2","2":"142.9117","3":"-3.952357","4":"0.8762287","5":"-4.510646","6":"6.463061e-06","7":"5.187163e-04"},{"1":"ENSG00000205810:KLRC3","2":"255.1010","3":"-3.976524","4":"0.6962522","5":"-5.711328","6":"1.120980e-08","7":"4.230943e-06"},{"1":"ENSG00000257842:NOVA1-AS1","2":"173.7307","3":"-3.992321","4":"0.6230043","5":"-6.408177","6":"1.472701e-10","7":"1.164129e-07"},{"1":"ENSG00000148357:HMCN2","2":"414.9299","3":"-4.093493","4":"0.8664833","5":"-4.724261","6":"2.309540e-06","7":"2.383140e-04"},{"1":"ENSG00000089199:CHGB","2":"545.4277","3":"-4.140743","4":"0.7936902","5":"-5.217077","6":"1.817691e-07","7":"3.661161e-05"},{"1":"ENSG00000114279:FGF12","2":"1904.6784","3":"-4.193712","4":"0.8583679","5":"-4.885682","6":"1.030716e-06","7":"1.312518e-04"},{"1":"ENSG00000164418:GRIK2","2":"3292.2787","3":"-4.204442","4":"0.8131319","5":"-5.170677","6":"2.332476e-07","7":"4.468485e-05"},{"1":"ENSG00000184221:OLIG1","2":"1148.5156","3":"-4.205053","4":"1.1100313","5":"-3.788229","6":"1.517249e-04","7":"5.580786e-03"},{"1":"ENSG00000104313:EYA1","2":"245.1158","3":"-4.226314","4":"0.8090343","5":"-5.223900","6":"1.751938e-07","7":"3.562984e-05"},{"1":"ENSG00000255641:NKG2-E","2":"345.8587","3":"-4.233324","4":"0.7239259","5":"-5.847732","6":"4.983214e-09","7":"2.108805e-06"},{"1":"ENSG00000196132:MYT1","2":"597.2625","3":"-4.254082","4":"0.7680167","5":"-5.539048","6":"3.041196e-08","7":"9.232674e-06"},{"1":"ENSG00000218475:RP11-93K7.1","2":"156.0799","3":"-4.254427","4":"0.8756861","5":"-4.858393","6":"1.183423e-06","7":"1.445466e-04"},{"1":"ENSG00000239519:CADM2-AS1","2":"170.4031","3":"-4.331102","4":"0.7471211","5":"-5.797054","6":"6.748991e-09","7":"2.772049e-06"},{"1":"ENSG00000183542:KLRC4","2":"171.7693","3":"-4.331346","4":"0.7122458","5":"-6.081251","6":"1.192483e-09","7":"7.036489e-07"},{"1":"ENSG00000224675:AC009227.2","2":"174.4691","3":"-4.370484","4":"0.9181592","5":"-4.760050","6":"1.935451e-06","7":"2.058013e-04"},{"1":"ENSG00000101204:CHRNA4","2":"112.0347","3":"-4.458768","4":"0.8788954","5":"-5.073150","6":"3.912835e-07","7":"6.557128e-05"},{"1":"ENSG00000205809:KLRC2","2":"373.8816","3":"-4.511848","4":"0.8111422","5":"-5.562340","6":"2.661817e-08","7":"8.448246e-06"},{"1":"ENSG00000104112:SCG3","2":"4646.7612","3":"-4.661433","4":"0.7624710","5":"-6.113587","6":"9.741625e-10","7":"5.914861e-07"},{"1":"ENSG00000101311:FERMT1","2":"1295.9922","3":"-4.740991","4":"0.8322592","5":"-5.696532","6":"1.222691e-08","7":"4.533153e-06"},{"1":"ENSG00000183090:FREM3","2":"202.0580","3":"-4.890584","4":"0.9229112","5":"-5.299084","6":"1.163850e-07","7":"2.677801e-05"},{"1":"ENSG00000160862:AZGP1","2":"305.0874","3":"-4.910138","4":"0.7593208","5":"-6.466486","6":"1.003078e-10","7":"8.404794e-08"},{"1":"ENSG00000145526:CDH18","2":"325.0296","3":"-4.965195","4":"1.3000447","5":"-3.819249","6":"1.338583e-04","7":"5.079706e-03"},{"1":"ENSG00000173157:ADAMTS20","2":"311.9936","3":"-5.049167","4":"0.9506486","5":"-5.311286","6":"1.088541e-07","7":"2.547732e-05"},{"1":"ENSG00000187398:LUZP2","2":"3093.3705","3":"-5.101987","4":"0.9637524","5":"-5.293878","6":"1.197497e-07","7":"2.726583e-05"},{"1":"ENSG00000164796:CSMD3","2":"2555.2303","3":"-5.143318","4":"0.9270974","5":"-5.547764","6":"2.893461e-08","7":"8.848290e-06"},{"1":"ENSG00000145451:GLRA3","2":"156.7315","3":"-5.155038","4":"0.9349241","5":"-5.513857","6":"3.510532e-08","7":"1.035725e-05"},{"1":"ENSG00000154975:CA10","2":"367.2088","3":"-5.246239","4":"1.2471240","5":"-4.206670","6":"2.591615e-05","7":"1.520668e-03"},{"1":"ENSG00000187416:LHFPL3","2":"2353.6227","3":"-5.448684","4":"0.7579027","5":"-7.189160","6":"6.519109e-13","7":"1.050454e-09"},{"1":"ENSG00000153993:SEMA3D","2":"566.7408","3":"-5.505097","4":"0.8283070","5":"-6.646203","6":"3.007508e-11","7":"3.073160e-08"},{"1":"ENSG00000122584:NXPH1","2":"431.2498","3":"-5.580754","4":"0.7583931","5":"-7.358656","6":"1.857705e-13","7":"3.891429e-10"},{"1":"ENSG00000267339:LINC00906","2":"108.8366","3":"-5.585157","4":"0.9932278","5":"-5.623239","6":"1.874097e-08","7":"6.435679e-06"},{"1":"ENSG00000205927:OLIG2","2":"2586.8728","3":"-5.596454","4":"0.7919706","5":"-7.066492","6":"1.588993e-12","7":"2.080340e-09"},{"1":"ENSG00000215612:HMX1","2":"106.4733","3":"-5.641662","4":"0.7887072","5":"-7.153051","6":"8.487021e-13","7":"1.269871e-09"},{"1":"ENSG00000131094:C1QL1","2":"2000.1352","3":"-5.778284","4":"0.5881267","5":"-9.824897","6":"8.796396e-23","7":"3.685250e-18"},{"1":"ENSG00000038295:TLL1","2":"253.7274","3":"-5.814804","4":"0.9382736","5":"-6.197343","6":"5.742415e-10","7":"3.759039e-07"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-09-28 14:01:28
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [5457cb7] 2020-09-28: Add plot for virtual 4c update
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
##  [1] cowplot_0.9.4 ggrepel_0.8.0 scales_1.1.1  ggplot2_3.1.0 purrr_0.3.4  
##  [6] glue_1.4.2    magrittr_1.5  dplyr_0.8.0   readr_1.3.1   tidyr_0.8.2  
## [11] here_0.1     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.5         git2r_0.27.1       pillar_1.4.6       compiler_3.5.1    
##  [5] RColorBrewer_1.1-2 plyr_1.8.6         tools_3.5.1        digest_0.6.25     
##  [9] evaluate_0.14      lifecycle_0.2.0    tibble_3.0.3       gtable_0.3.0      
## [13] pkgconfig_2.0.3    rlang_0.4.7        yaml_2.2.1         xfun_0.17         
## [17] withr_2.2.0        stringr_1.4.0      knitr_1.29         vctrs_0.3.4       
## [21] hms_0.5.3          rprojroot_1.3-2    grid_3.5.1         tidyselect_1.1.0  
## [25] R6_2.4.1           rmarkdown_1.11     farver_2.0.3       codetools_0.2-15  
## [29] backports_1.1.9    ellipsis_0.3.1     htmltools_0.5.0    assertthat_0.2.1  
## [33] colorspace_1.4-1   stringi_1.5.3      lazyeval_0.2.2     munsell_0.5.0     
## [37] crayon_1.3.4
```

</details>


***

<!-- END OF END MATTER -->
