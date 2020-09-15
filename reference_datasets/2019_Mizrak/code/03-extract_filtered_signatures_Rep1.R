# Selin Jessa
# 2020-07-07
#
# Here I will get 100-gene signatures for the clusters, filtered 
# based on gene stats

library(tidyr)
library(dplyr)
library(glue)
library(readr)
library(magrittr)
library(ggplot2)
library(purrr)

out <- "../processed_data/03"

# Load gene annotation
anno <- read_tsv("../../../sjessa/from_hydra/single_cell/2019-02_revision1/input/gene_annotation_with_homologous_mice_genes.txt")

# Get markers from paper supplementary, and add cell type labels
markers <- read_tsv("../processed_data/01-cluster_markers.tsv") %>% 
  dplyr::rename(mmusculus_external_gene_id = gene,
                Cluster = cluster)

# Annotate the genes with mouse/human homologs and ENSEMBL IDs
markers_anno <- markers %>% 
  left_join(select(anno,
                   hsapiens_external_gene_id,
                   mmusculus_external_gene_id,
                   hsapiens_ensembl_gene_id,
                   mmusculus_ensembl_gene_id = mmusculus_homolog_ensembl_gene),
            by = c("mmusculus_external_gene_id")) %>%
  distinct(Cluster, mmusculus_external_gene_id, .keep_all = TRUE)

# Save the table of annotated markers
write_tsv(markers_anno, glue("{out}-markers.annotated.tsv"))

# Filter to positive LFC markers; rank by adj pval, and take the top 100
# as we did for the Nat Genet 2019 paper
filterMarkers <- function(markers, gene_col, sp = "mm", n_top = 100) {
  
  gene_col <- rlang::enquo(gene_col)
  
  markers %>% 
    {
      if (sp == "mm") dplyr::filter(., !grepl("Rps|Rpl|Mrps|Mrpl", !!gene_col))
      else dplyr::filter(., !grepl("RPS|RPL|MRPS|MRPL", !!gene_col))
    } %>% 
    group_by(Cluster) %>% 
    dplyr::filter(avg_logFC > 0) %>% 
    arrange(p_val_adj) %>% 
    dplyr::slice(1:n_top) %>% 
    ungroup()
  
}

# Do the filtering
sigs <- markers_anno %>% 
  filterMarkers(mmusculus_external_gene_id, sp = "mm")

# Check # markers per cluster
sort(table(sigs$Cluster))

# Filter out signatures with < 80 genes and clusters with < 100 cells
(enough_genes <- sigs %>% count(Cluster) %>% filter(n >= 70) %>% pull(Cluster))

sigs2 <- sigs %>% filter(Cluster %in% enough_genes)

# Sanity check
sort(table(sigs2$Cluster))

sigs2_qc <- sigs2 %>%
  group_by(Cluster) %>% 
  summarise(median_logFC = median(avg_logFC),
            n_gene_FC_above_1.5 = sum(avg_logFC > log2(1.5))) %>% 
  arrange(median_logFC) %T>% 
  readr::write_tsv(glue("{out}-sigs_qc.tsv"))

sigs2_qc %>% 
  mutate(Cluster = factor(Cluster, levels = unique(.$Cluster))) %>% 
  gather(key = "stat", value = "value", 2:length(.)) %>%
  ggplot(aes(x = Cluster, y = value)) +
  geom_bar(aes(fill = Cluster), stat = "identity", width = 0.8) +
  facet_wrap(~ stat, ncol = 2, scales = "free_x") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

# Don't filter anything for now, although the clustering resolution seems a bit low
sigs3 <- sigs2

# Extract the signatures, filter NA values
signatures <- list("hg_sym" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, hsapiens_external_gene_id) %>% .[!is.na(.)]),
                   "mm_sym" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, mmusculus_external_gene_id) %>% .[!is.na(.)]),
                   "hg_ens" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, hsapiens_ensembl_gene_id) %>% .[!is.na(.)]),
                   "mm_ens" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, mmusculus_ensembl_gene_id) %>% .[!is.na(.)]))

saveRDS(signatures, file = glue("{out}-cluster_gene_signatures.Rds"))