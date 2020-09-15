# Selin Jessa
# 2020-01-12
#
# Here I will get 100-gene signatures for the clusters, filtered 
# based on gene stats

library(tidyverse)
library(readr)

out <- "../processed_data/01"

anno <- readr:: read_tsv("../../../sjessa/from_hydra/single_cell/2019-02_revision1/input/gene_annotation_with_homologous_mice_genes.txt")
markers <- readxl::read_xls("../data/Data-S3_cell_type_markers.xls") %>% 
  rename(Cluster = `Cell type`,
         hsapiens_ensembl_gene_id = `Gene ID`,
         hsapiens_external_gene_id = `Gene name`,
         avg_logFC = FC,
         p_val_adj = FDR)

markers_anno <- markers %>% 
  left_join(select(anno,
                   hsapiens_external_gene_id,
                   mmusculus_external_gene_id,
                   hsapiens_ensembl_gene_id,
                   mmusculus_ensembl_gene_id = mmusculus_homolog_ensembl_gene),
            by = c("hsapiens_external_gene_id", "hsapiens_ensembl_gene_id")) %>%
  distinct(Cluster, hsapiens_external_gene_id, .keep_all = TRUE)

write_tsv(markers_anno, glue("{out}-markers.annotated.tsv"))

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

sigs <- markers_anno %>% 
  filterMarkers(hsapiens_external_gene_id, sp = "hg")

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
  theme_min() +
  theme(legend.position = "none")

# Not filtering any signatures for now
sigs3 <- sigs2

signatures <- list("hg_sym" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, hsapiens_external_gene_id)),
                   "mm_sym" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, mmusculus_external_gene_id)),
                   "hg_ens" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, hsapiens_ensembl_gene_id)),
                   "mm_ens" = split(sigs3, f = sigs3$Cluster) %>% map(~ pull(.x, mmusculus_ensembl_gene_id)))

saveRDS(signatures, file = glue("{out}-cluster_gene_signatures.Rds"))

data.frame(Cluster = names(signatures$hg_sym)) %>% 
  write_tsv(glue("{out}-cluster_names.tsv"))
