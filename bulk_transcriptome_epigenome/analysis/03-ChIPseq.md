---
title: "03 - Epigenome / ChIP-seq analysis"
date: "17 September, 2020"
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
## Cache: ~/tmp/G34-gliomas/bulk_transcriptome_epigenome/03/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/bulk_transcriptome_epigenome/output/03
```

```
## G34-gliomas/bulk_transcriptome_epigenome/figures/03//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Here we take the G34 patient tumor epigenomic data (genic promoter scores for H3K27me3 and H3K27ac) and cross these values with the differentially-expressed genes by RNAseq, to identify genes with concordant differences
at the transcriptomic and epigenomic level.


# Libraries


```r
# General
library(tidyr)
library(readr)
library(dplyr)
library(magrittr)
library(glue)

# For plotting
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
ggplot2::theme_set(theme_min(base_size = 13))

# Custom
source(here::here(subdir, "analysis/functions.R"))
```

# Analysis

## Load data

Load tumor data:


```r
dge <- read_tsv(here(subdir, "output/01/DGE_HGG-IDHvsHGG-G34R.V_padj<0.01_baseMean>100_absLFC>2.tsv")) %>% 
  separate(ID, into = c("ENSID", "symbol"), sep = ":")
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

Load and tidy epigenomics data:


```r
epi <- read_tsv(here(subdir, "data/20191112.RefSeq.K27me3K27ac.txt"))
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   chr = col_character(),
##   strand = col_character(),
##   name = col_character(),
##   ID = col_character(),
##   X71 = col_logical()
## )
```

```
## See spec(...) for full column specifications.
```


Clean up gene names converted to dates:


```r
epi %>% distinct(name) %>% pull(name) %>% sort() %>% head(25)
```

```
##  [1] "01-Dec" "01-Mar" "01-Sep" "02-Mar" "02-Sep" "03-Mar" "03-Sep" "04-Mar"
##  [9] "04-Sep" "05-Mar" "05-Sep" "06-Mar" "06-Sep" "07-Mar" "07-Sep" "08-Mar"
## [17] "08-Sep" "09-Mar" "09-Sep" "10-Mar" "10-Sep" "11-Mar" "11-Sep" "12-Sep"
## [25] "14-Sep"
```

```r
epi <- epi %>%
  mutate(name = case_when(
    name == "01-Mar" ~ "MARCH1",
    name == "02-Mar" ~ "MARCH2",
    name == "03-Mar" ~ "MARCH3",
    name == "04-Mar" ~ "MARCH4",
    name == "05-Mar" ~ "MARCH5",
    name == "06-Mar" ~ "MARCH6",
    name == "07-Mar" ~ "MARCH7",
    name == "08-Mar" ~ "MARCH8",
    name == "09-Mar" ~ "MARCH9",
    name == "10-Mar" ~ "MARCH10",
    name == "11-Mar" ~ "MARCH11",
    name == "01-Sep" ~ "SEPT1",
    name == "02-Sep" ~ "SEPT2",
    name == "03-Sep" ~ "SEPT3",
    name == "04-Sep" ~ "SEPT4",
    name == "05-Sep" ~ "SEPT5",
    name == "06-Sep" ~ "SEPT6",
    name == "07-Sep" ~ "SEPT7",
    name == "08-Sep" ~ "SEPT8",
    name == "09-Sep" ~ "SEPT9",
    name == "10-Sep" ~ "SEPT10",
    name == "11-Sep" ~ "SEPT11",
    name == "12-Sep" ~ "SEPT12",
    name == "14-Sep" ~ "SEPT14",
    name == "01-Dec" ~ "DEC1",
    TRUE ~ name))
```


```r
epi_mean <- epi %>%
  select(1:7, G34RV.K27me3.median, nonG34RV.K27me3.median, G34RV.K27ac.median, nonG34RV.K27ac.median, G34R.zK27m3, G34R.zK27ac) %>% 
  gather(Stat, Score, 8:ncol(.)) %>% 
  group_by(name, Stat) %>% 
  # Take the mean over different transcripts (although note that for each transcript,
  # the score is the same anyway)
  summarize(Score = mean(Score)) %>% 
  spread(Stat, Score) %>% 
  ungroup()

epi_indiv <- epi %>% 
  select(name, 8:40)

epi_meta <- data.frame(Sample = colnames(epi)[8:42],
                       Mark = c(rep("K27me3", 18), rep("K27Ac", 17)),
                       Group = c(rep("G34R/V", 5), rep("IDH/SETD2", 7), rep("WT", 6),
                                 rep("G34R/V", 6), rep("IDH/SETD2", 6), rep("WT", 5)))

rr_saveRDS(desc = "Median and z-scored values of promoter H3K27me3/ac for G34 and non-G34 entities",
           object = epi_mean,
           file = glue("{out}/genic_promoters_histone_summary.Rds"))
```

```
## ...writing description of genic_promoters_histone_summary.Rds to G34-gliomas/bulk_transcriptome_epigenome/output/03/genic_promoters_histone_summary.desc
```

## Crossing epigenomics with other data

### DGE


```r
# Join the differentially expressed genes table and the epigenomics data
# table based on gene symbol
dge_and_epi <- dge %>%
  left_join(epi_mean, by = c("symbol" = "name")) %>% 
  arrange(abs(log2FoldChange)) %T>% 
  rr_write_tsv(path = glue("{out}/DGE_and_histone_marks.tsv"),
               desc = "Integrated RNA-seq differential expression and histone mark ChIPseq values for patient tumor samples")
```

```
## ...writing description of DGE_and_histone_marks.tsv to G34-gliomas/bulk_transcriptome_epigenome/output/03/DGE_and_histone_marks.desc
```

```r
dge_and_epi %>% 
  ggplot(aes(x = G34R.zK27ac, y = G34R.zK27m3)) +
  geom_hline(yintercept = 0.5, size = 0.5, colour = "gray90") +
  geom_vline(xintercept = 0.5, size = 0.5, colour = "gray90") +
  geom_point(aes(colour = log2FoldChange, size = -log10(padj)), alpha = 0.8) +
  scale_colour_gradientn(colours = colorRampPalette(c("navy", "white", "red"))(n = 200)) +
  cowplot::theme_cowplot()
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/03//dge_and_epi-1.png)<!-- -->

```r
goi <- dge_and_epi %>% filter(G34R.zK27ac > 0.9 | G34R.zK27m3 > 0.9) %>%
  filter(abs(log2FoldChange) > 3) %>%
  pull(symbol)

dge_and_epi %>% 
  rr_ggplot(aes(x = G34R.zK27ac, y = G34R.zK27m3), plot_num = 2) +
  geom_hline(yintercept = 0.5, size = 0.5, colour = "gray90") +
  geom_vline(xintercept = 0.5, size = 0.5, colour = "gray90") +
  geom_point(aes(colour = log2FoldChange, size = -log10(padj)), alpha = 0.8) +
  scale_colour_gradientn(colours = colorRampPalette(c("navy", "white", "red"))(n = 200)) +
  ggrepel::geom_label_repel(data = dge_and_epi %>% filter(symbol %in% goi), aes(label = symbol))
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/03//dge_and_epi-2.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/03//dge_and_epi...*]~</span>


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-09-17 13:17:46
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [918687c] 2020-09-17: Update TOC links
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
## [1] ggrepel_0.8.0      RColorBrewer_1.1-2 ggplot2_3.1.0      glue_1.4.2        
## [5] magrittr_1.5       dplyr_0.8.0        readr_1.3.1        tidyr_0.8.2       
## [9] here_0.1          
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.5       git2r_0.27.1     pillar_1.4.6     compiler_3.5.1  
##  [5] plyr_1.8.6       tools_3.5.1      digest_0.6.25    evaluate_0.14   
##  [9] lifecycle_0.2.0  tibble_3.0.3     gtable_0.3.0     pkgconfig_2.0.3 
## [13] rlang_0.4.7      yaml_2.2.1       xfun_0.17        withr_2.2.0     
## [17] stringr_1.4.0    knitr_1.29       vctrs_0.3.4      hms_0.5.3       
## [21] cowplot_0.9.4    rprojroot_1.3-2  grid_3.5.1       tidyselect_1.1.0
## [25] R6_2.4.1         rmarkdown_1.11   farver_2.0.3     purrr_0.3.4     
## [29] codetools_0.2-15 backports_1.1.9  scales_1.1.1     ellipsis_0.3.1  
## [33] htmltools_0.5.0  assertthat_0.2.1 colorspace_1.4-1 labeling_0.3    
## [37] stringi_1.5.3    lazyeval_0.2.2   munsell_0.5.0    crayon_1.3.4
```

</details>


***

<!-- END OF END MATTER -->
