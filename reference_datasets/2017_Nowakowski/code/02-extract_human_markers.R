# Selin Jessa
# 2019-11-26
#
# Here I will get 100-gene signatures for the clusters, filtered the same
# way we did our atlas

library(tidyverse)

anno <- readr:: read_tsv("../../../sjessa/from_hydra/single_cell/2019-02_revision1/input/gene_annotation_with_homologous_mice_genes.txt")

# 0. Update filt ----
# Filter by fold change since we don't have adj p-val
filterMarkers <- function(markers, n_top = 100) {
  
  markers %>%
    dplyr::filter(., !grepl("RPS|RPL|MRPS|MRPL", hsapiens_external_gene_id)) %>% 
                    group_by(Cluster) %>% 
    dplyr::filter(avg_logFC > 0) %>% 
    arrange(desc(avg_logFC)) %>% 
    dplyr::slice(1:n_top) %>% 
    ungroup()
  
}


# 1. markers ----
nowakowski_markers <- read_csv("../data/nowakowski2017_tableS5_cluster_markers.csv") %>% 
  filter(avg_diff > 0) %>% 
  select(hsapiens_external_gene_id = gene,
         p_val_adj = p_val,
         avg_logFC = avg_diff,
         cluster,
         pct.1,
         pct.2)

# Populate additional columns. Since we don't have an adjusted p-value, we copied
# the p-value to that column, but won't use it for ordering genes.

nowa_anno <- nowakowski_markers %>% 
  left_join(select(anno, hsapiens_external_gene_id, mmusculus_external_gene_id,
                   hsapiens_ensembl_gene_id, mmusculus_ensembl_gene_id = mmusculus_homolog_ensembl_gene),
                         by = c("hsapiens_external_gene_id")) %>%
  distinct(cluster, mmusculus_external_gene_id, .keep_all = TRUE)

write_tsv(nowa_anno, "../processed_data/nowakowski_markers.annotated.tsv")

nowa_sigs <- nowa_anno %>% 
  rename(Cluster = cluster) %>% 
  filterMarkers()

# Double check we get 100 markers per cluster
all(table(nowa_sigs$Cluster) == 100)
# [1] TRUE

nowakowski_signatures <- list("hg_sym" = split(nowa_sigs, f = nowa_sigs$Cluster) %>% map(~ pull(.x, hsapiens_external_gene_id)),
                              "mm_sym" = split(nowa_sigs, f = nowa_sigs$Cluster) %>% map(~ pull(.x, mmusculus_external_gene_id)),
                              "hg_ens" = split(nowa_sigs, f = nowa_sigs$Cluster) %>% map(~ pull(.x, hsapiens_ensembl_gene_id)),
                              "mm_ens" = split(nowa_sigs, f = nowa_sigs$Cluster) %>% map(~ pull(.x, mmusculus_ensembl_gene_id)))

save(nowakowski_signatures, file = "../processed_data/nowakowski_signatures.Rda")
