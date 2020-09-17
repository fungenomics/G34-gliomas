# 2020-03-02
# Prepare tables to accompany the paper for bulk RNAseq

# ---------------------------------------------------------------------------- #

doc_id <- "06"
out <- paste0("output/", doc_id); dir.create(out, recursive = TRUE)

# 0. Set-up ----

library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(glue)
library(here)

# 1.1. Sequencing stats and QC ----

# First assemble the metadata & stats for all samples included

# This includes all cell line RNAseq used in any analysis
cl_meta    <- read_tsv("output/08/cell_line_metadata_with_align_stats.tsv") %>% 
  filter(!(`Cell Line` == "KNS-42" & Genotype_broad == "KO/Repair")) %>% 
  mutate(mitochondrialPercentage = paste0(mitochondrialPercentage, "%"))

# This includes all tumor samples used in any analysis
tumor_meta <- read_xlsx("../metadata/2020-02-25_bulkRNAseq_metadata.xlsx")

# Get alignment stats for tumors
tumor_align_stats <- read_tsv("../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/All_tumor_samples/alignment.statistics.tsv")

tumor_meta <- left_join(tumor_meta, tumor_align_stats, by = c("Sample" = "name"))

seq_stats_to_keep <- quos(rawReadCount,
                          cleanReadRetentionPercentage,
                          mappedPercentage,
                          wholeGenePercentage,
                          mitochondrialPercentage,
                          ribosomalPercentage)

# Harmonize the columns
tumor_meta2 <- tumor_meta %>% 
  rename(Internal_sample = Sample) %>% 
  mutate(Source = case_when(
    is.na(IC_Patient) ~ "Jessa et al, Nature Genetics, 2019",
    TRUE ~ "This study"
  )) %>% 
  mutate(Sample = case_when(
    is.na(IC_Patient) ~ Sample_anon_ID_NatGenet,
    TRUE ~ paste0(IC_Patient, "_", IC_Sample)
  )) %>% 
  # Include the Internal_sample for tracking purposes, but DELETE from final table
  select(Internal_sample, Sample, Source, Group, Genotype_G34, Genotype_PDGFRA, Location_broad = Location, Location = Location2,
         !!! seq_stats_to_keep) %>% 
  mutate(Tissue = "Tumor")

cl_meta2 <- cl_meta %>% 
  select(Sample = `Cell Line`, Clone, Condition = Group, Genotype, Media,
         !!! seq_stats_to_keep) %>% 
  # Tidy the nomenclature
  mutate(Tissue = "Patient-derived cell line",
         Media = recode(Media, "nsm" = "Stem cell media", "ser" = "Serum"),
         Condition = case_when(
           grepl("nsm", Condition) ~ gsub("nsm", "stem", Condition),
           grepl("ser", Condition) ~ gsub("ser", "serum", Condition)
         ),
         Source = "This study")

# Combine tumor and cell lines
table_bulk_rnaseq_stats <- bind_rows(tumor_meta2, cl_meta2) %>%
  select(Sample, Source, Tissue, Group, Genotype_G34, Genotype_PDGFRA, Location_broad, Location,
         Clone, Condition, Genotype, Media,
         everything())

write_tsv(table_bulk_rnaseq_stats, glue("{out}/table_bulk_rnaseq_stats.tsv"))


# 1.2. Counts ----
tumor_counts <- "../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/All_tumor_samples/counts/Ensembl.ensGene.exon.raw.tsv.gz" %>% 
  read.table(header = T, sep = "\t", check.names = FALSE)

table_bulk_rnaseq_counts <- tumor_counts
colnames(table_bulk_rnaseq_counts) <- tumor_meta2 %>% data.frame() %>% tibble::column_to_rownames(var = "Internal_sample") %>% .[colnames(tumor_counts), "Sample"]
table_bulk_rnaseq_counts <- table_bulk_rnaseq_counts %>% tibble::rownames_to_column(var = "ID")
write_tsv(table_bulk_rnaseq_counts, glue("{out}/table_bulk_rnaseq_counts.tsv"))

# 1.3. DGE ----
# Tumors
tumor_raptor_path <- "../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/"
comps_t <- list.files(tumor_raptor_path, full.names = TRUE)[grepl("^HGG", list.files(tumor_raptor_path)) & grepl("vs", list.files(tumor_raptor_path)) & grepl("batch_covariate", list.files(tumor_raptor_path))]
names(comps_t) <- gsub("_batch_covariate", "", basename(comps_t))
comps_t <- comps_t[c("HGG-G34R.V_vs_HGG-IDH",
                     "HGG-G34R.V_vs_HGG-H3.3-K27M-pons",
                     "HGG-G34R.V_vs_HGG-H3.3-K27M-thalamic",
                     "HGG-G34R.V_vs_HGG-WT",
                     "HGG-G34R.V_vs_PNET-HGNET-BCOR")]

# Cell lines
cl_raptor_path <- "../../../from_beluga/HGG/2019-09_bulk_RNAseq/2020-01_G34_submission1_add_samples/Cell_lines3"
comps_cl <- list.files(cl_raptor_path, full.names = TRUE)[grepl("^GBM002", list.files(cl_raptor_path)) & !grepl("vs", list.files(cl_raptor_path))]
names(comps_cl) <- basename(comps_cl)

comps <- c(comps_t, comps_cl)

# > names(comps)
# [1] "HGG-G34R.V_vs_HGG-H3.3-K27M-pons"     "HGG-G34R.V_vs_HGG-H3.3-K27M-thalamic" "HGG-G34R.V_vs_HGG-IDH"               
# [4] "HGG-G34R.V_vs_HGG-WT"                 "HGG-G34R.V_vs_PNET-HGNET-BCOR"        "GBM002_nsm"                          
# [7] "GBM002_ser"    

table_dge <- imap_dfc(comps, ~ read_tsv(Sys.glob(file.path(.x, "diff/Ensembl.ensGene.exon/*.tsv"))) %>% 
                        mutate(Comparison = .y) %>% 
                        arrange(ID)) %>% 
  arrange(pvalue)

write_tsv(table_dge, glue("{out}/table_dge.tsv"))

# 1.4. GSEA ----
load("output/02/fgsea_df.Rda")

table_gsea <- fgsea_df %>% 
  filter(padj < 0.01) %>% 
  # Clean up the source of signatures
  mutate(Reference_dataset = 
           case_when(
             !is.na(Sample) ~ "Jessa et al 2019",
             Age == "Human fetal" ~ "Nowakowski et al 2017",
             Age == "Human ped/adult" ~ "Velmeshev et al 2019"
           )) %>% 
  select(Comparison, Cell_type, Age, Reference_dataset, pval, padj, ES, NES, leadingEdge) %>% 
  mutate(Comparison = factor(Comparison, levels = c(
    "HGG-G34R.V_vs_HGG-IDH",
    "HGG-G34R.V_vs_HGG-WT",
    "HGG-G34R.V_vs_HGG-H3.3-K27M-pons",
    "HGG-G34R.V_vs_HGG-H3.3-K27M-thalamic",
    "HGG-G34R.V_vs_PNET-HGNET-BCOR"
  ))) %>% 
  arrange(Comparison, desc(NES))

write_tsv(table_gsea, glue("{out}/table_gsea.tsv"))
