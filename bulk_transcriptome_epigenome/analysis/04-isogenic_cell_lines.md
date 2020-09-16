---
title: "01 - Bulk RNAseq"
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
## Document index: 04
```

```r
# Specify where to save outputs
out        <- here(subdir, "output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here(subdir, "figures", doc_id, "/"); dir.create(figout, recursive = TRUE)
cache      <- paste0("~/tmp/", basename(here()), "/", subdir, "/", doc_id, "/")

message("Cache: ", cache)
```

```
## Cache: ~/tmp/G34-gliomas/bulk_transcriptome_epigenome/04/
```

</details>

The root directory of this project is:

```
## /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
```

Outputs and figures will be saved at these paths, relative to project root:

```
## G34-gliomas/bulk_transcriptome_epigenome/output/04
```

```
## G34-gliomas/bulk_transcriptome_epigenome/figures/04//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

In this document, we analyze bulk RNA-seq data from the G34-mutant patient-derived cell line GMB002 and isogenic matches, in which CRISPR editing was used to remove the mutation. Cell lines were also subjected to a differentiation protocol.


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
library(cowplot)
ggplot2::theme_set(theme_min(base_size = 13))

# Custom
source(here::here(subdir, "analysis/functions.R"))
source(here::here(subdir, "analysis/ssgsea.R"))
```

# Analysis

## Prep bulk RNAseq analysis

### Sample metadata

Load in the metadata table for all bulk RNAseq samples to be included
in the analysis:


```r
(cell_line_meta <- readxl::read_xlsx(here(subdir, "data/20200219-cell_line_metadata.xlsx")) %>% 
   mutate(Group_broad = paste0(Genotype_broad, "_", Media)) %>% 
   mutate(Group_broad = gsub("/", "", Group_broad)))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Nickname"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Lab"],"name":[3],"type":["chr"],"align":["left"]},{"label":["Library Strategy"],"name":[4],"type":["chr"],"align":["left"]},{"label":["Run Type"],"name":[5],"type":["chr"],"align":["left"]},{"label":["Read Length"],"name":[6],"type":["chr"],"align":["left"]},{"label":["Cell Line"],"name":[7],"type":["chr"],"align":["left"]},{"label":["Clone"],"name":[8],"type":["chr"],"align":["left"]},{"label":["Passage"],"name":[9],"type":["chr"],"align":["left"]},{"label":["Media"],"name":[10],"type":["chr"],"align":["left"]},{"label":["Genotype"],"name":[11],"type":["chr"],"align":["left"]},{"label":["Genotype_broad"],"name":[12],"type":["chr"],"align":["left"]},{"label":["Replicate"],"name":[13],"type":["chr"],"align":["left"]},{"label":["R1.fastq"],"name":[14],"type":["chr"],"align":["left"]},{"label":["R2.fastq"],"name":[15],"type":["chr"],"align":["left"]},{"label":["Group"],"name":[16],"type":["chr"],"align":["left"]},{"label":["Batch"],"name":[17],"type":["chr"],"align":["left"]},{"label":["Group_broad"],"name":[18],"type":["chr"],"align":["left"]}],"data":[{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/CMC1118G8_p11","2":"CMC1118G8_p11","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"CMC1118G8","8":"Parental","9":"p11","10":"ser","11":"H3F3A(+/G34R)","12":"G34R","13":"1","14":"HI.5043.008.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.CMC1118G8_p11.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.CMC1118G8_p11.pair2.fastq.gz","16":"H3F3A(+/G34R) ser","17":"New","18":"G34R_ser"},{"1":"njabado/RNAseq/2017-01-09-PROJECT13824_RUN3892_NJ-RNAseq-clinical-Nov_2016_hg/HSJD-GBM002-p42","2":"429-HSJD-GBM002-p42","3":"Jabado","4":"RiboZero","5":"__NA__","6":"__NA__","7":"HSJD-GBM002","8":"Parental","9":"__NA__","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"__NA__","14":"__NA__","15":"__NA__","16":"H3F3A(+/G34R) nsm","17":"Old","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_A10_p58","2":"HSJD-GBM002_A10_p58","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A10","9":"p58","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"1","14":"HI.5043.007.NEBNext_dual_i7_C11---NEBNext_dual_i5_C11.HSJD-GBM002_A10_p58.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_C11---NEBNext_dual_i5_C11.HSJD-GBM002_A10_p58.pair2.fastq.gz","16":"H3F3A(+/G34R) nsm","17":"New","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_A10_p59","2":"HSJD-GBM002_A10_p59","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A10","9":"p59","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"2","14":"HI.5043.008.NEBNext_dual_i7_B12---NEBNext_dual_i5_B12.HSJD-GBM002_A10_p59.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_B12---NEBNext_dual_i5_B12.HSJD-GBM002_A10_p59.pair2.fastq.gz","16":"H3F3A(+/G34R) nsm","17":"New","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cA10_p58","2":"GBM002_cA10_p58","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A10","9":"p58","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"3","14":"NS.1209.002.NEBNext_dual_i7_B10---NEBNext_dual_i5_B10.GBM002_cA10_p58_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_B10---NEBNext_dual_i5_B10.GBM002_cA10_p58_R2.fastq.gz","16":"H3F3A(+/G34R) nsm","17":"New","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_F06_p55","2":"HSJD-GBM002_F06_p55","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"F06","9":"p55","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"1","14":"HI.5043.007.NEBNext_dual_i7_B11---NEBNext_dual_i5_B11.HSJD-GBM002_F06_p55.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_B11---NEBNext_dual_i5_B11.HSJD-GBM002_F06_p55.pair2.fastq.gz","16":"H3F3A(+/G34R) nsm","17":"New","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_F06_p56","2":"HSJD-GBM002_F06_p56","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"F06","9":"p56","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"2","14":"HI.5043.008.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.HSJD-GBM002_F06_p56.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.HSJD-GBM002_F06_p56.pair2.fastq.gz","16":"H3F3A(+/G34R) nsm","17":"New","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cF06_p64","2":"GBM002_cF06_p64","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"F06","9":"p64","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"3","14":"NS.1209.002.NEBNext_dual_i7_A10---NEBNext_dual_i5_A10.GBM002_cF06_p64_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_A10---NEBNext_dual_i5_A10.GBM002_cF06_p64_R2.fastq.gz","16":"H3F3A(+/G34R) nsm","17":"New","18":"G34R_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cA10_p58_SerDif","2":"GBM002_cA10_p58_SerDif","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A10","9":"p58","10":"ser","11":"H3F3A(+/G34R)","12":"G34R","13":"1","14":"NS.1209.002.NEBNext_dual_i7_F10---NEBNext_dual_i5_F10.GBM002_cA10_p58_SerDif_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_F10---NEBNext_dual_i5_F10.GBM002_cA10_p58_SerDif_R2.fastq.gz","16":"H3F3A(+/G34R) ser","17":"New","18":"G34R_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cF06_p64_SerDif","2":"GBM002_cF06_p64_SerDif","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"F06","9":"p64","10":"ser","11":"H3F3A(+/G34R)","12":"G34R","13":"1","14":"NS.1209.002.NEBNext_dual_i7_E10---NEBNext_dual_i5_E10.GBM002_cF06_p64_SerDif_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_E10---NEBNext_dual_i5_E10.GBM002_cF06_p64_SerDif_R2.fastq.gz","16":"H3F3A(+/G34R) ser","17":"New","18":"G34R_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_A09_p58","2":"HSJD-GBM002_A09_p58","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A09","9":"p58","10":"nsm","11":"H3F3A(+/+)","12":"KO/Repair","13":"1","14":"HI.5043.007.NEBNext_dual_i7_E11---NEBNext_dual_i5_E11.HSJD-GBM002_A09_p58.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_E11---NEBNext_dual_i5_E11.HSJD-GBM002_A09_p58.pair2.fastq.gz","16":"H3F3A(+/+) nsm","17":"New","18":"KORepair_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_A09_p59","2":"HSJD-GBM002_A09_p59","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A09","9":"p59","10":"nsm","11":"H3F3A(+/+)","12":"KO/Repair","13":"2","14":"HI.5043.008.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.HSJD-GBM002_A09_p59.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.HSJD-GBM002_A09_p59.pair2.fastq.gz","16":"H3F3A(+/+) nsm","17":"New","18":"KORepair_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cC08_p54","2":"GBM002_cC08_p54","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A09","9":"p54","10":"nsm","11":"H3F3A(+/+)","12":"KO/Repair","13":"3","14":"NS.1209.002.NEBNext_dual_i7_D10---NEBNext_dual_i5_D10.GBM002_cC08_p54_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_D10---NEBNext_dual_i5_D10.GBM002_cC08_p54_R2.fastq.gz","16":"H3F3A(+/+) nsm","17":"New","18":"KORepair_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_C08_p54","2":"HSJD-GBM002_C08_p54","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"C08","9":"p54","10":"nsm","11":"H3F3A(+/-)","12":"KO/Repair","13":"1","14":"HI.5043.007.NEBNext_dual_i7_D11---NEBNext_dual_i5_D11.HSJD-GBM002_C08_p54.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_D11---NEBNext_dual_i5_D11.HSJD-GBM002_C08_p54.pair2.fastq.gz","16":"H3F3A(+/-) nsm","17":"New","18":"KORepair_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/HSJD-GBM002_C08_p56","2":"HSJD-GBM002_C08_p56","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"C08","9":"p56","10":"nsm","11":"H3F3A(+/-)","12":"KO/Repair","13":"2","14":"HI.5043.008.NEBNext_dual_i7_C12---NEBNext_dual_i5_C12.HSJD-GBM002_C08_p56.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_C12---NEBNext_dual_i5_C12.HSJD-GBM002_C08_p56.pair2.fastq.gz","16":"H3F3A(+/-) nsm","17":"New","18":"KORepair_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cA09_p58","2":"GBM002_cA09_p58","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"C08","9":"p58","10":"nsm","11":"H3F3A(+/-)","12":"KO/Repair","13":"3","14":"NS.1209.002.NEBNext_dual_i7_C10---NEBNext_dual_i5_C10.GBM002_cA09_p58_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_C10---NEBNext_dual_i5_C10.GBM002_cA09_p58_R2.fastq.gz","16":"H3F3A(+/-) nsm","17":"New","18":"KORepair_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cA09_p58_SerDif","2":"GBM002_cA09_p58_SerDif","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"A09","9":"p58","10":"ser","11":"H3F3A(+/+)","12":"KO/Repair","13":"1","14":"NS.1209.002.NEBNext_dual_i7_G10---NEBNext_dual_i5_G10.GBM002_cA09_p58_SerDif_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_G10---NEBNext_dual_i5_G10.GBM002_cA09_p58_SerDif_R2.fastq.gz","16":"H3F3A(+/+) ser","17":"New","18":"KORepair_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/GBM002_cC08_p58_SerDif","2":"GBM002_cC08_p58_SerDif","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"HSJD-GBM002","8":"C08","9":"p58","10":"ser","11":"H3F3A(+/-)","12":"KO/Repair","13":"1","14":"NS.1209.002.NEBNext_dual_i7_H10---NEBNext_dual_i5_H10.GBM002_cC08_p58_SerDif_R1.fastq.gz","15":"NS.1209.002.NEBNext_dual_i7_H10---NEBNext_dual_i5_H10.GBM002_cC08_p58_SerDif_R2.fastq.gz","16":"H3F3A(+/-) ser","17":"New","18":"KORepair_ser"},{"1":"njabado/RNAseq/2015-08-18-cell_lines/KNS42-C18-H3_3G34V","2":"238-KNS42-C18-H3_3G34V","3":"Jabado","4":"RiboZero","5":"__NA__","6":"__NA__","7":"KNS-42","8":"Parental","9":"__NA__","10":"nsm","11":"H3F3A(+/G34V)","12":"G34V","13":"__NA__","14":"__NA__","15":"__NA__","16":"H3F3A(+/G34R) nsm","17":"Old","18":"G34V_nsm"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/KNS42_parent_p34","2":"KNS42_parent_p34","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"KNS-42","8":"Parental","9":"p34","10":"ser","11":"H3F3A(+/G34V)","12":"G34V","13":"1","14":"HI.5043.007.NEBNext_dual_i7_F11---NEBNext_dual_i5_F11.KNS42_parent_p34.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_F11---NEBNext_dual_i5_F11.KNS42_parent_p34.pair2.fastq.gz","16":"H3F3A(+/G34V) ser","17":"New","18":"G34V_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/KNS42_G34V-KO_c1-16_p42","2":"KNS42_G34V-KO_c1-16_p42","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"KNS-42","8":"c1-16","9":"p42","10":"ser","11":"H3F3A(+/-)","12":"KO/Repair","13":"1","14":"HI.5043.008.NEBNext_dual_i7_E12---NEBNext_dual_i5_E12.KNS42_G34V-KO_c1-16_p42.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_E12---NEBNext_dual_i5_E12.KNS42_G34V-KO_c1-16_p42.pair2.fastq.gz","16":"H3F3A(+/-) ser","17":"New","18":"KORepair_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/KNS42_G34V-KO_c1-2_p41","2":"KNS42_G34V-KO_c1-2_p41","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"KNS-42","8":"c1-2","9":"p41","10":"ser","11":"H3F3A(+/-)","12":"KO/Repair","13":"1","14":"HI.5043.007.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.KNS42_G34V-KO_c1-2_p41.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.KNS42_G34V-KO_c1-2_p41.pair2.fastq.gz","16":"H3F3A(+/-) ser","17":"New","18":"KORepair_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/KNS42_G34V-KO_c1-9_p42","2":"KNS42_G34V-KO_c1-9_p42","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"KNS-42","8":"c1-9","9":"p42","10":"ser","11":"H3F3A(+/-)","12":"KO/Repair","13":"1","14":"HI.5043.007.NEBNext_dual_i7_H11---NEBNext_dual_i5_H11.KNS42_G34V-KO_c1-9_p42.pair1.fastq.gz","15":"HI.5043.007.NEBNext_dual_i7_H11---NEBNext_dual_i5_H11.KNS42_G34V-KO_c1-9_p42.pair2.fastq.gz","16":"H3F3A(+/-) ser","17":"New","18":"KORepair_ser"},{"1":"njabado/RNAseq/2019-11-23-G34-lines-RNAseq-@sjessa/KNS42_G34V-KO_c2-2_p40","2":"KNS42_G34V-KO_c2-2_p40","3":"Jabado","4":"riboMinus","5":"PAIRED_END","6":"run_100","7":"KNS-42","8":"c2-2","9":"p40","10":"ser","11":"H3F3A(+/-)","12":"KO/Repair","13":"1","14":"HI.5043.008.NEBNext_dual_i7_F12---NEBNext_dual_i5_F12.KNS42_G34V-KO_c2-2_p40.pair1.fastq.gz","15":"HI.5043.008.NEBNext_dual_i7_F12---NEBNext_dual_i5_F12.KNS42_G34V-KO_c2-2_p40.pair2.fastq.gz","16":"H3F3A(+/-) ser","17":"New","18":"KORepair_ser"},{"1":"njabado/RNAseq/2015-09-07-batch2/PS10-801_CL","2":"244-PS10-801_CL","3":"Jabado","4":"RiboZero","5":"__NA__","6":"__NA__","7":"PS10-801","8":"Parental","9":"__NA__","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"__NA__","14":"__NA__","15":"__NA__","16":"H3F3A(+/G34R) nsm","17":"Old","18":"G34R_nsm"},{"1":"njabado/RNAseq/2016-02-16-PROJECT12674_RUN3489_NJ-RNAseq-GBM-liver-Jan-11-2016/PS10-801_C13","2":"319-PS10-801_C13","3":"Jabado","4":"RiboZero","5":"__NA__","6":"__NA__","7":"PS10-801","8":"Parental","9":"__NA__","10":"nsm","11":"H3F3A(+/G34R)","12":"G34R","13":"__NA__","14":"__NA__","15":"__NA__","16":"H3F3A(+/G34R) nsm","17":"Old","18":"G34R_nsm"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# Number of samples per group
table(cell_line_meta$Group)
```

```
## 
##    H3F3A(+/-) nsm    H3F3A(+/-) ser    H3F3A(+/+) nsm    H3F3A(+/+) ser 
##                 3                 5                 3                 1 
## H3F3A(+/G34R) nsm H3F3A(+/G34R) ser H3F3A(+/G34V) ser 
##                10                 3                 1
```

### Prepare bulk RNA-seq pipeline input

Palettes for groups:


```r
palette_groups <- c("G34R_ser" = "cyan4",
                    "G34R_nsm" = "cyan1",
                    "G34V_ser" = "cyan4",
                    "G34V_nsm" = "cyan1",
                    "KORepair_ser" = "navy",
                    "KORepair_nsm" = "blue")

save(palette_groups, file = glue("{out}/cl_palettes.Rda"))
```

We will need one version with all samples together, since normalization
is dataset-dependent, and we may want to compare expression levels across groups.


```r
dir.create(file.path(out, "All_cell_line_data"), showWarnings = FALSE)

cell_line_meta %>%
  select(ID, Nickname, Group = Group_broad) %>%
  write_tsv(file.path(out, "All_cell_line_data", "info.samples.tsv"))

data.frame(Group = unique(cell_line_meta$Group_broad),
           Label = unique(cell_line_meta$Group_broad),
           Order = seq_along(unique(cell_line_meta$Group_broad)),
           Color = palette_groups[unique(cell_line_meta$Group_broad)],
           stringsAsFactors = FALSE) %>%
  write_tsv(file.path(out, "All_cell_line_data", "info.groups.tsv"))
```


## Run pipeline on Beluga

In order to run the in-house DESeq2/expression analysis pipeline, the inputs
created above are synced to the pipeline folder, and then
the following outputs are loaded for further analysis below and in the subsequent analyses:

- Alignment stats
- Differential gene expession analysis output
- Counts


```r
pipeline_path <- "../../../2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/Cell_lines3/"
```

## Alignment stats


```r
align_stats <- read_tsv(file.path(pipeline_path, "All_cell_line_data", "alignment.statistics.tsv")) %>%
  separate(mitochondrialPercentage, sep = "%", into = "mitochondrialPercentage") %>%
  mutate(mitochondrialPercentage = as.numeric(mitochondrialPercentage))
```

```
## Parsed with column specification:
## cols(
##   .default = col_character(),
##   rawReadCount = col_double(),
##   cleanReadCount = col_double(),
##   unmappedReadCount = col_double(),
##   mappedReadCount = col_double(),
##   wholeGeneReadCount = col_double(),
##   exonReadCount = col_double(),
##   intronReadCount = col_double(),
##   utr5ReadCount = col_double(),
##   cdsReadCount = col_double(),
##   utr3ReadCount = col_double(),
##   mitochondrialReadCount = col_double(),
##   ribosomalReadCount = col_double()
## )
```

```
## See spec(...) for full column specifications.
```

```r
# DT::datatable(align_stats, filter = "top", class = 'white-space: nowrap')

left_join(cell_line_meta, align_stats, by = c("Nickname" = "name")) %>%
  write_tsv(glue("{out}/cell_line_metadata_with_align_stats.tsv"))
```

## G34R vs KO / repair in serum

### Load DGE


```r
dge_serum <- read_tsv(file.path(pipeline_path, "GBM002_ser/diff/Ensembl.ensGene.exon/G34R_servsKORepair_ser.tsv")) %>% 
  separate(ID, into = c("ENSID", "gene_symbol"), sep = ":") %>% 
  mutate(log10p = -log10(padj))
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
dge_serum %>% 
  filter(baseMean > 100 & padj < 0.05) %>% 
  rr_ggplot(aes(x = log2FoldChange, y = log10p, label = gene_symbol), alpha = 0.8, plot_num = 1) +
  geom_point()
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/04//dge_serum-1.png)<!-- -->

```r
dge_stem <- read_tsv(file.path(pipeline_path, "GBM002_nsm/diff/Ensembl.ensGene.exon/G34R_nsmvsKORepair_nsm.tsv")) %>% 
  separate(ID, into = c("ENSID", "gene_symbol"), sep = ":") %>% 
  mutate(log10p = -log10(padj))
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


### Targeted DGE, for serum condition - GBM002

Note that I have taken care, here, to use the counts normalized
for the serum samples together only.


```r
goi <- c("DLX1", "DLX2", "DLX5", "DLX6", "PDGFRA", "GSX2")

gbm002_counts <- extract_pipeline_counts(file.path(pipeline_path, "GBM002_ser/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                       goi = goi) %>%
  left_join(cell_line_meta, by = c("sample" = "Nickname")) %>% 
  left_join(dge_serum, by = c("gene_symbol" = "gene_symbol", "gene_ensg" = "ENSID"))
```

```
## Joining, by = "gene_symbol"
```

```r
dge_serum %>% filter(gene_symbol %in% goi)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ENSID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene_symbol"],"name":[2],"type":["chr"],"align":["left"]},{"label":["baseMean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["log2FoldChange"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["lfcSE"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["stat"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["pvalue"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["padj"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["log10p"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"ENSG00000105880","2":"DLX5","3":"658.09433","4":"1.5905839","5":"0.3301889","6":"4.8171945","7":"1.455907e-06","8":"0.0001727677","9":"3.76253756"},{"1":"ENSG00000144355","2":"DLX1","3":"15755.71196","4":"0.5199826","5":"0.1370992","6":"3.7927471","7":"1.489899e-04","8":"0.0066446211","9":"2.17752978"},{"1":"ENSG00000134853","2":"PDGFRA","3":"15148.73844","4":"-0.6204999","5":"0.1868425","6":"-3.3209778","7":"8.970267e-04","8":"0.0246858090","9":"1.60755264"},{"1":"ENSG00000115844","2":"DLX2","3":"18568.20544","4":"0.3318288","5":"0.1283568","6":"2.5852057","7":"9.732091e-03","8":"0.1303521990","9":"0.88488164"},{"1":"ENSG00000006377","2":"DLX6","3":"713.25248","4":"1.3350848","5":"0.6551704","6":"2.0377672","7":"4.157322e-02","8":"0.3330453073","9":"0.47749668"},{"1":"ENSG00000180613","2":"GSX2","3":"79.72791","4":"0.4182405","5":"0.5458985","6":"0.7661508","7":"4.435866e-01","8":"0.8930264147","9":"0.04913569"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

For the DLXs only


```r
dlx_limits <- list("DLX1"   = c(0, 21000),
                   "DLX2"   = c(0, 21000),
                   "DLX5"   = c(0, 1600),
                   "DLX6"   = c(0, 1600),
                   "PDGFRA" = c(0, 22000),
                   "GSX2"   = c(0, 200))

cl_dotplot <- function(gene, lims, counts, levels = c("G34R_nsm", "KORepair_nsm",
                                                      "G34R_ser", "KORepair_ser")) {
  
  counts %>%
    filter(gene_symbol == gene) %>% 
    mutate(Group_broad = factor(Group_broad,
                                levels = levels)) %>% 
    ggplot(aes(x = Group_broad, y = gene_expression)) +
    geom_jitter(aes(colour = Group_broad), size = 5, alpha = 0.87, width = 0.1) +
    scale_colour_manual(values = c("cyan4", "midnightblue")) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.2) +
    theme_min(border_colour = "black") +
    ylim(lims) +
    ggtitle(gene) +
    rotateX() +
    noLegend()
  
}

gbm002_counts %>% 
  select(sample, gene_symbol, gene_expression, Group_broad) %>% 
  write_tsv(glue("{figout}/dlx_serum_boxplot-1.source_data.tsv"))

imap(dlx_limits, ~ cl_dotplot(.y, .x, counts = gbm002_counts)) %>% 
  {plot_grid(plotlist = ., ncol = 6)}
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/04//dlx_serum_boxplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/04//dlx_serum_boxplot...*]~</span>


### Targeted DGE, for stem condition - GBM002


```r
gbm002_counts_nsm <- extract_pipeline_counts(file.path(pipeline_path, "GBM002_nsm/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                           goi = goi) %>%
  left_join(cell_line_meta, by = c("sample" = "Nickname")) %>% 
  left_join(dge_stem, by = c("gene_symbol" = "gene_symbol", "gene_ensg" = "ENSID"))
```

```
## Joining, by = "gene_symbol"
```

```r
dge_stem %>% filter(gene_symbol %in% goi)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ENSID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene_symbol"],"name":[2],"type":["chr"],"align":["left"]},{"label":["baseMean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["log2FoldChange"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["lfcSE"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["stat"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["pvalue"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["padj"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["log10p"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"ENSG00000006377","2":"DLX6","3":"284.06987","4":"0.542042546","5":"0.16101295","6":"3.36645312","7":"0.0007614152","8":"0.04816356","9":"1.3172813866"},{"1":"ENSG00000105880","2":"DLX5","3":"247.82526","4":"0.514908241","5":"0.25068774","6":"2.05398255","7":"0.0399773806","8":"0.45437972","9":"0.3425810605"},{"1":"ENSG00000134853","2":"PDGFRA","3":"20045.91584","4":"0.214689417","5":"0.16504960","6":"1.30075700","7":"0.1933416441","8":"0.80481776","9":"0.0943024495"},{"1":"ENSG00000115844","2":"DLX2","3":"4129.20301","4":"0.085125191","5":"0.12504680","6":"0.68074669","7":"0.4960317901","8":"0.97508622","9":"0.0109569832"},{"1":"ENSG00000180613","2":"GSX2","3":"16.66248","4":"0.278731316","5":"0.41677302","6":"0.66878446","7":"0.5036329829","8":"0.97670982","9":"0.0102344449"},{"1":"ENSG00000144355","2":"DLX1","3":"3380.09662","4":"0.004416948","5":"0.09802634","6":"0.04505879","7":"0.9640604498","8":"0.99902440","9":"0.0004239044"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
dlx_limits <- list("DLX1"   = c(0, 4000),
                   "DLX2"   = c(0, 5300),
                   "DLX5"   = c(0, 400),
                   "DLX6"   = c(0, 500),
                   "PDGFRA" = c(0, 30000),
                   "GSX2"   = c(0, 30))

gbm002_counts_nsm %>% 
  select(sample, gene_symbol, gene_expression, Group_broad) %>% 
  write_tsv(glue("{figout}/dlx_stem_boxplot-1.source_data.tsv"))

imap(dlx_limits, ~ cl_dotplot(.y, .x, counts = gbm002_counts_nsm)) %>% 
  {plot_grid(plotlist = ., ncol = 6)}
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/04//dlx_stem_boxplot-1.png)<!-- --><br><span style="color:#0d00ff">~[figure/source data @ *G34-gliomas/bulk_transcriptome_epigenome/figures/04//dlx_stem_boxplot...*]~</span>


## ssGSEA

### Prep counts


```r
raw_counts <- file.path(pipeline_path, "All_GBM002_new/counts/Ensembl.ensGene.exon.raw.tsv.gz") %>%
  read.table(header = T, row.names = 1, check.names = F, sep = "\t")

head(rownames(raw_counts) %<>% strsplit(":") %>% sapply(getElement, 1))
```

```
## [1] "ENSG00000118473" "ENSG00000162426" "ENSG00000157191" "ENSG00000169504"
## [5] "ENSG00000142920" "ENSG00000186094"
```

```r
raw_counts <- as.matrix(raw_counts)

raw_counts[1:5, 1:5]
```

```
##                 HSJD-GBM002_A10_p58 HSJD-GBM002_A10_p59 GBM002_cA10_p58
## ENSG00000118473                  37                  30              80
## ENSG00000162426                  90                  81             174
## ENSG00000157191                 475                 473            1017
## ENSG00000169504                7542                6570           19732
## ENSG00000142920                  55                  48              89
##                 HSJD-GBM002_F06_p55 HSJD-GBM002_F06_p56
## ENSG00000118473                 114                  62
## ENSG00000162426                  80                  40
## ENSG00000157191                 531                 364
## ENSG00000169504               11481                5515
## ENSG00000142920                  66                  20
```

```r
save(raw_counts, file = glue("{out}/cl_raw_counts.Rda"))

dim(raw_counts)
```

```
## [1] 60234    16
```


### Load signatures

Loading the cell type gene signatures prepared previously for GSEA, and selecting
some representative signatures for relevant cell types:


```r
signatures <- readRDS(here(subdir, "output/02/signatures_ens.Rds"))

# Select pathways of interest
poi <- list("newborn_inhib_neuron" = signatures$`HF nIN1`,
            "mature_inhib_neuron" = signatures$`HP IN-PV`,
            "opc" = signatures$`F-p6 OPC`,
            "astrocyte" = signatures$`F-p0 ASTR2`,
            "rgc" = signatures$`HF RG-early`,
            "fetal_excit_neurons" = signatures$`HF EN-PFC1`)
```

### Run ssGSEA


```r
cl_ssgsea <- ssgsea_le(expr_mat = raw_counts,
                       gene_sets = poi,
                       # Following exactly the same parameters as we used
                       # in the NG bulk analysis
                       alpha = 0.75,
                       normalize = FALSE)

cl_ssgsea$enrichment_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Signature") %>%
  write_tsv(glue("{out}/bulk_ssgsea_scores.tsv"))

leading_edge <- cl_ssgsea$leading_edge
save(leading_edge, file = glue("{out}/cl_ssgsea_le.Rda"))
```

### Visualize results


```r
cl_ssgsea_tidy <- cl_ssgsea$enrichment_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Signature") %>%
  gather(Sample, Score, 2:ncol(.)) %>%
  left_join(cell_line_meta, by = c("Sample" = "Nickname"))

cl_ssgsea_tidy %>%
  filter(Signature %in% c("newborn_inhib_neuron", "rgc")) %>%
  mutate(Signature = factor(Signature, levels = c("rgc", "newborn_inhib_neuron"))) %>%
  filter(Media == "ser") %>%
  mutate(Group_broad = factor(Group_broad,
                              levels = c("G34R_nsm", "KORepair_nsm",
                                         "G34R_ser", "KORepair_ser"))) %>%
  rr_ggplot(aes(x = Group_broad, y = Score), plot_num = 1) +
  geom_jitter(aes(colour = Group_broad), size = 5, alpha = 0.87, width = 0.1) +
  scale_colour_manual(values = c("cyan4", "midnightblue")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.2) +
  theme_min(border_colour = "black") +
  facet_wrap(~ Signature) +
  rotateX() +
  noLegend() +
  ylim(c(22400, 24000))
```

![](/lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas/bulk_transcriptome_epigenome/figures/04//viz_ssgsea-1.png)<!-- -->



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2020-09-16 18:25:59
```

The git repository and last commit:

```
## Local:    master /lustre03/project/6004736/sjessa/from_beluga/HGG-G34/G34-gliomas
## Remote:   master @ origin (git@github.com:fungenomics/G34-gliomas.git)
## Head:     [250bd7d] 2020-09-16: Actually update lockfile, with monocle this time...
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
##  [1] pbapply_1.4-3 cowplot_0.9.4 ggrepel_0.8.0 scales_1.1.1  ggplot2_3.1.0
##  [6] purrr_0.3.4   glue_1.4.2    magrittr_1.5  dplyr_0.8.0   readr_1.3.1  
## [11] tidyr_0.8.2   here_0.1     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.5          git2r_0.27.1        plyr_1.8.6         
##  [4] pillar_1.4.6        compiler_3.5.1      RColorBrewer_1.1-2 
##  [7] BiocManager_1.30.10 tools_3.5.1         digest_0.6.25      
## [10] evaluate_0.14       lifecycle_0.2.0     tibble_3.0.3       
## [13] gtable_0.3.0        pkgconfig_2.0.3     rlang_0.4.7        
## [16] parallel_3.5.1      yaml_2.2.1          xfun_0.17          
## [19] withr_2.2.0         stringr_1.4.0       knitr_1.29         
## [22] vctrs_0.3.4         hms_0.5.3           rprojroot_1.3-2    
## [25] grid_3.5.1          tidyselect_1.1.0    R6_2.4.1           
## [28] rmarkdown_1.11      farver_2.0.3        codetools_0.2-15   
## [31] backports_1.1.9     ellipsis_0.3.1      htmltools_0.5.0    
## [34] assertthat_0.2.1    colorspace_1.4-1    renv_0.10.0        
## [37] labeling_0.3        stringi_1.5.3       lazyeval_0.2.2     
## [40] munsell_0.5.0       crayon_1.3.4
```

</details>


***

<!-- END OF END MATTER -->
