
library(pbapply)
library(magrittr)

#' integrate_cdfs
#' 
#' This code was adapted from GSVA and the GSEA code from the Broad Institute
#' 
#' @param gene_set Character vector of genes
#' @param gene_ranking Numeric vector specifying gene ranking
#' @param R Matrix of gene ranks
#' @param j Integer, specifying sample (column of expression matrix)
#' @param alpha Numeric, exponent giving the weight in the weighted ECDF $P^w_{in}$
integrate_cdfs <- function(gene_set, gene_ranking, expr_mat, R, j, alpha) {
  
  # Binary 1/0 vector indicating whether each gene (in ranked order)
  # is in the gene set or not
  indicator_in_gene_set <- match(rownames(expr_mat)[gene_ranking], gene_set)
  indicator_in_gene_set[!is.na(indicator_in_gene_set)] <- 1
  indicator_in_gene_set[ is.na(indicator_in_gene_set)] <- 0
  indicator_out <- 1 - indicator_in_gene_set
  
  # P_in (the weighted ECDF of genes in the set)
  # i.e. it will take a weighted step whenever a gene is in the set
  P_in <- cumsum( (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha ) /
    sum( (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha )
  
  # P_out (the ECDF of genes not in the set)
  # i.e. it will be length N - Ng, with a step whenever the gene is not in the set
  P_out <- cumsum( !indicator_in_gene_set ) /
    sum( !indicator_in_gene_set )
  
  # The difference in CDFs
  cdf_diffs <- P_in - P_out
  es <- sum(cdf_diffs)
  
  # Calculate running sum statistic, and leading edge subset
  n_g <- length(gene_set)
  n <- nrow(R)
  steps_in <- (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha
  step_norm_in <- 1/sum(steps_in)
  step_norm_out <- 1/(n - n_g)
  
  running_sum <- cumsum(indicator_in_gene_set * steps_in * step_norm_in -  # Increment the score by a weighted step if a gene is a hit
                          indicator_out * step_norm_out)                   # Decrement by an fixed size step if a gene is a miss
  
  # The leading edge consists of all the genes in the gene set, which come up
  # in the ranked list prior to the max enrichment score (diff. between ECDFs)
  leading_edge <- names(running_sum)[1:which.max(running_sum)]
  
  # We only care about the ones in the gene set
  leading_edge_indicator <- 1 * ((rownames(R) %in% leading_edge) & (rownames(R) %in% gene_set))
  names(leading_edge_indicator) <- rownames(R)
  
  return(list(
    leading = leading_edge_indicator,
    es = sum(cdf_diffs),
    running = running_sum))
  
}

#' This code was adapted from GSVA
ssgsea_le <- function(expr_mat, gene_sets, alpha = 0.25, normalize = TRUE,
                      save_le = TRUE, n_cores = 1, verbose = FALSE) {
  
  require(pbapply)
  
  if (verbose) message("1. Converting gene expr to ranks...")
  
  # Convert gene expression data in X to ranks within each sample
  R <- apply(expr_mat, 2, function(x) as.integer(rank(x, na.last = TRUE)))
  rownames(R) <- rownames(expr_mat)
  
  # # Collect the differences between CDFs
  # diff_list <- vector("list", length = n)
  
  if (verbose) message(glue("2. Calculating enrichment, parallelizing over ", n_cores, " cores..."))
  
  # For each sample S_1 to S_n, calculate the score for each gene set
  # Parallelize over samples
  ssgsea_out <- pblapply(1:ncol(expr_mat), function(j) {
    
    # The order (by indices) by which to retrieve genes from the matrix
    # to get them from high -> low
    
    # BUG: I think that subsetting the dataframe in this way does not properly
    # deal with any NAs for that sample
    gene_ranking <- order(R[, j], decreasing = TRUE, na.last = TRUE)
    sample_out <- lapply(gene_sets, integrate_cdfs, gene_ranking, expr_mat, R, j, alpha)
    
    return(sample_out)
    
  }, cl = n_cores)
  
  # TODO: this is super unwieldy and unneeded!
  # Would be better to simply return the genes in the set which are in the leading
  # edge as a char vector, that could be saved.
  # But, that may not be the solution, since the leading edge of interest,
  # in our case, is defined per-sample, and so we'll only really care about it
  # with respect to biological groups of samples
  filter_zeros <- function(df) {
    
    df[rowSums(df) != 0, ]
    
  }
  
  if (verbose) message("3. Tidying ouput...")
  
  if (save_le) {
    
    # Tidy the binary data indicating if genes are within the leading edge subset
    leading_mat <- lapply(ssgsea_out,
                          function(sample_out) purrr::transpose(sample_out)$leading) %>% 
      magrittr::set_names(colnames(expr_mat)) %>%
      purrr::transpose() %>%
      lapply(bind_rows) %>% 
      lapply(as.data.frame) %>% 
      lapply(magrittr::set_rownames, rownames(expr_mat)) %>% 
      lapply(filter_zeros)
    
    # Tidy the running sum data
    # Can't coerce to a dataframe because then order would be lost
    running_sum <- lapply(ssgsea_out,
                          function(sample_out) purrr::transpose(sample_out)$running) %>% 
      magrittr::set_names(colnames(expr_mat)) %>%
      purrr::transpose()
    
  } else {
    
    leading_mat <- NULL
    
  }
  
  # Tidy the enrichment scores for each set for each sample
  es <- lapply(ssgsea_out, function(sample_out) unlist(map(sample_out, "es"))) %>%
    data.frame %>%
    as.matrix()
  
  # Normalize the whole dataset
  normES <- function(mat) apply(mat, 2, function(x, mat) x /
                                  (range(mat)[2] - range(mat)[1]), mat)
  
  if (normalize) es <- normES(es)
  
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(expr_mat)
  
  return(list("enrichment_scores" = es,
              "leading_edge" = leading_mat))
  
  if (verbose) message("4. Done.")
  
}


#' Compute ssGSEA scores on bulk data given SC-derived signatures,
#' using custom function which returns leading edge
#' 
#' Given path to cluster marker input files,
#' reformat the matrix, and compute the scores with the
#' bulk RNA-seq data. Code adapted from Nicolas De Jay.
#'
#' @param marker_path character Path to marker incidence matrix
#' @param counts list Matrix of counts from the RNA-seq data
#' @param out_prefix character Prefix for path/file name where scores should be saved as tsv
#' will be appended with ".ssgsea_scores.txt"
compute_scores_le <- function(marker_path, counts, out_prefix,
                              gene_col = "hsapiens_ensembl_gene_id",
                              test_mode = FALSE,
                              save_le = TRUE,
                              n_cores = 1) {
  
  require(pbapply)
  
  if (is.character(marker_path)) {
    
    markers_mat <- read_tsv(marker_path)
    
    markers <- markers_mat %>%
      .[, 6:ncol(.)] %>%
      sapply(function (i) markers_mat[[gene_col]][i], 
             simplify = FALSE)
    
  } else if (is.list(marker_path)) markers <- marker_path
  
  
  if (test_mode) { # Just run for two samples, to make sure it works
    
    out <- ssgsea_le(expr_mat  = counts[,c(1, 2)],
                     gene_sets = markers,
                     alpha     = 0.75,
                     normalize = FALSE,
                     save_le   = save_le,
                     n_cores   = 1)
    
    return(out)
    
  }
  
  out <- ssgsea_le(expr_mat  = counts,
                   gene_sets = markers,
                   alpha     = 0.75,
                   normalize = FALSE,
                   save_le   = save_le,
                   n_cores   = n_cores)
  
  out$enrichment_scores %>%
    data.frame %>% 
    magrittr::set_colnames(colnames(counts)) %>% 
    tibble::rownames_to_column(var = "Signature") %>%
    {write.table(file = glue("{out_prefix}.ssgsea_scores.txt"),
                 x = .,
                 sep = "\t",
                 quote = TRUE,
                 row.names = FALSE)}
  
  if (save_le) {
    
    leading_edge <- out$leading_edge
    save(leading_edge, file = glue("{out_prefix}.ssgsea_le.Rda"))  
    
  }
  
  
  invisible(out$enrichment_scores)
  
}



#' @param genelists A named list, where each element is a character vector of
#' genes. Expects gene symbols.
#' @param summary_type Character specifying how to treat single cells for 
#' the enrichment. If "cluster_mean", computes cluster mean expression and
#' evaluates enrichment per cluster. If "metacells", assumes the seurat object
#' contains metacells (by \code{seurat2metacells}) and evaluates enrichment
#' per metacell.
#' @param gene_space Character vector with gene symbols, used to restrict
#' the space of genes considered in ssGSEA. Note, this is *not* the list of genes
#' for which enrichment will be assessed. Rather, the enrichment of gene set
#' will be computed among these genes.
#' 
#' @return A dataframe, with the first column "Cell" (where the order matches
#' the order in \code{seurat}), and following columns corresponding to signatures,
#' named as in the \code{genelists} input.
compute_scores_sc_le <- function(seurat, genelists,
                                 out_prefix = "./",
                                 summary_type = "cluster_mean",
                                 n_cores = 1,
                                 save_le = FALSE,
                                 gene_space = NULL,
                                 return_counts = FALSE) {
  
  require(pbapply)
  
  # Get counts
  if (summary_type == "cluster_mean") counts <- seurat %>% cytokit::meanClusterExpression()
  else if (summary_type %in% c("metacells", "single_cells")) {
    
    counts <- as.matrix(seurat@data)
    
  }
  
  # Reduce counts to the space of genes to consider
  if (!is.null(gene_space)) {
    gene_space_detected <- unique(gene_space[gene_space %in% rownames(counts)])
    counts <- counts[gene_space_detected, ]
    message("0. Reducing the gene space to the ", nrow(counts), " genes detectd from the provided genes")
  }
  
  if (return_counts) return(counts)
  
  out <- ssgsea_le(expr_mat  = counts,
                   gene_sets = genelists,
                   alpha     = 0.75,
                   normalize = FALSE,
                   save_le   = save_le,
                   n_cores   = n_cores)
  
  # Tidy scores and write to file
  out$enrichment_scores %>%
    data.frame %>% 
    magrittr::set_rownames(NULL) %>% 
    magrittr::set_rownames(rownames(out$enrichment_scores)) %>% 
    magrittr::set_colnames(colnames(counts)) %>% 
    tibble::rownames_to_column(var = "Signature") %>%
    {write.table(file = glue("{out_prefix}.{summary_type}.ssgsea_scores.txt"),
                 x = .,
                 sep = "\t",
                 quote = TRUE,
                 row.names = FALSE)}
  
  if (isTRUE(save_le)) {
    
    leading_edge <- out$leading_edge
    save(leading_edge, file = glue("{out_prefix}.{summary_type}.ssgsea_le.Rda"))  
    
  }
  
  if (length(genelists) == 1) return(unlist(out$enrichment_scores)[1, ])
  
}


#' @describeIn compute_scores_le
compute_scores_sc_split <- function(seurat,
                                    clusters = levels(seurat@ident),
                                    out_prefix,
                                    genelists,
                                    n_cores = 1,
                                    save_le = FALSE,
                                    gene_space = NULL) {
  for (i in seq_along(clusters)) {
    
    seurat_sub <- Seurat::SubsetData(seurat, ident.use = clusters[i])
    compute_scores_sc_le(seurat_sub,
                         out_prefix = paste0(out_prefix, ".cluster", clusters[i]),
                         genelists = genelists,
                         n_cores = n_cores,
                         save_le = save_le,
                         summary_type = "single_cells",
                         gene_space = gene_space)
  }
  
  # Gather output
  
  cluster_files <- list.files(dirname(out_prefix),
                              pattern = paste0(basename(out_prefix), "\\..*.single_cells.ssgsea_scores.txt"), full.names = TRUE)
  
  # purrr::map(cluster_files,
  #            ~ .x %>% read_tsv %>% as.data.frame %>% tibble::column_to_rownames(var = "Signature")) %>% 
  #   {do.call(cbind, .)} %>%
  #   tibble::rownames_to_column(var = "Signature") %>%
  #   {write.table(file = glue("{out_prefix}.single_cells.ssgsea_scores.txt"),
  #                x = .,
  #                sep = "\t",
  #                quote = TRUE,
  #                row.names = FALSE)}
  
  purrr::map(cluster_files,
                        ~ .x %>% read_tsv %>% as.data.frame %>% tibble::column_to_rownames(var = "Signature")) %>% 
    {do.call(cbind, .)} %>%
    t() %>% 
    # Ensure the order matches the Seurat object
    .[seurat@cell.names, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Cell") %>% 
    {write.table(file = glue("{out_prefix}.single_cells.ssgsea_scores.txt"),
                 x = .,
                 sep = "\t",
                 quote = TRUE,
                 row.names = FALSE)}
  
  file.remove(cluster_files)
  
}


#' @describeIn compute_scores_le
compute_scores_sc_bin <- function(seurat,
                                  n_bins = 10,
                                  out_prefix,
                                  genelists,
                                  n_cores = 1,
                                  save_le = FALSE,
                                  gene_space = NULL) {
  
  seurat@meta.data$Bin <- cut(1:length(seurat@cell.names), labels = FALSE, n_bins)
  seurat <- SetAllIdent(seurat, "Bin")
  compute_scores_sc_split(seurat,
                          out_prefix = out_prefix,
                          genelists = genelists,
                          n_cores = 1,
                          save_le = save_le,
                          gene_space = gene_space)
  
}





#' Get leading edge genes
#'
#' For a given signature (cluster), return the genes which are
#' in the leading edge subset in some proportion of the samples. Since more, rather
#' than fewer, genes in the gene set tend to be in the leading edge, the default
#' proportion is 1, i.e. the most stringent.
#' 
#' @param leading_edge The leading edge list object produced by \code{compute_scores_le},
#' either a character vector with the path, or the list itself.
get_leading_edge <- function(leading_edge,
                             cluster,
                             samples = NULL,
                             proportion = 1,
                             return_mat = FALSE,
                             return_symbols = TRUE,
                             leading_edge_only = TRUE) {
  
  if (is.character(leading_edge)) load(leading_edge)
  
  if (!(cluster %in% names(leading_edge))) stop("Cluster not found in signatures...")
  
  if (is.null(samples)) le_clust <- leading_edge[[cluster]]
  else if (length(samples) == 1) {
    
    # Simply return the leading edge genes, by definition, for one sample  
    return(rownames(leading_edge[[cluster]])[which(leading_edge[[cluster]][, samples] == 1)])
    
  } else {
    
    # Subset to the leading edge info for all relevant samples
    le_clust <- leading_edge[[cluster]][, samples]
    message("...Identifying leading edge genes for ", ncol(le_clust), " samples")
    
  }
  
  if (return_mat) return(le_clust)
  
  if (leading_edge_only) le_genes <- le_clust[rowSums(le_clust) >= proportion * ncol(le_clust), ] %>%
      rownames()
  else le_genes <- le_clust %>% rownames()
  
  if (return_symbols & grepl("^ENS", head(le_genes[[1]]))) return(ensembl2symbols(le_genes))
  else return(le_genes)
  
  
}



ensembl2symbols_safe <- function(genes, ...) {
  
  sym <- lapply(genes, ensembl2symbols, ...)
  empty <- which(unlist(map(sym, ~ identical(unlist(.x), character(0)))))
  sym[empty] <- genes[empty]
  
  return(unlist(sym))
  
}
