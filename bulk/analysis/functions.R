# Helper functions for the bulk RNAseq/epigenome analysis


#' Based off of code from Nicolas De Jay. Extract counts for genes of interest from
#' the raptor bulk RNA-seq pipeline and wrangle them into tidy format
#'
#' @param path Character, path to counts file
#' @param goi Character, vector of one or more gene symbols for genes of interest
#' for which to retrieve counts
#'
#' @return Data frame with four columns: sample, gene_expression, gene_ensg, gene_symbol
extract_pipeline_counts <- function(path, goi) {
    
    counts <- read.table(path, header = T, sep = "\t", check.names = FALSE)
    
    # Our pipeline labels genes using #{ENSG}:#{SYMBOL}
    ensg_to_symbol <- data.frame("gene_id" = rownames(counts) %>% as.character) %>%
        dplyr::mutate("gene_ensg"   = gene_id %>% as.character %>% strsplit(":") %>% sapply(`[[`, 1)) %>%
        dplyr::mutate("gene_symbol" = gene_id %>% as.character %>% strsplit(":") %>% sapply(`[[`, 2))
    
    # Convert all genes from #{SYMBOL} to #{ENSG}:#{SYMBOL}
    genes.symbol <- data.frame(genes = goi)
    colnames(genes.symbol) = "gene_symbol"
    
    genes.ensg <- genes.symbol %>% left_join(ensg_to_symbol) %>% dplyr::select(gene_id)
    
    # Subset counts table
    counts.subset <- counts[genes.ensg %>% pull(gene_id) %>% as.character, , drop = F] %>%
        as.matrix() %>%
        reshape2::melt() %>%
        setNames(c("gene_id", "sample", "gene_expression")) %>%
        left_join(ensg_to_symbol, by = "gene_id") %>%
        dplyr::select(-gene_id)
    
    return(counts.subset)
    
}


#' Convert Ensembl gene IDs to gene symbols
ensembl2symbols_safe <- function(genes, ...) {
    
    sym <- lapply(genes, ensembl2symbols, ...)
    empty <- which(unlist(map(sym, ~ identical(unlist(.x), character(0)))))
    sym[empty] <- genes[empty]
    
    return(unlist(sym))
    
}


# Helpers for GSEA analysis ----------------------------------------------------

get_gene_ranks <- function(pathway, stats, direction = "enriched") {
    
    if (direction == "enriched") rnk <- rank(-stats)
    else if (direction == "depleted") rnk <- rank(stats)
    
    all_genes_ranked <- sort(rnk)
    pathway_genes_ranked <- all_genes_ranked[pathway[pathway %in% names(stats)]]
    
    x <- data.frame(gene = names(pathway_genes_ranked),
                    rank = unname(pathway_genes_ranked)) %>%
        arrange(rank)
    
    x
    
}