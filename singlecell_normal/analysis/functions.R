



filter_NA <- function(x) {

    x[!is.na(x)]

}


#' A distance function for Spearman rank correlation
#' Copied from https://davetang.org/muse/2010/11/26/hierarchical-clustering-with-p-values-using-spearman/
spearman <- function(x, ...) {
  x <- as.matrix(x)
  res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
  res <- as.dist(res)
  attr(res, "method") <- "spearman"
  return(res)
}




prop <- function(x) {
    
    sum(x)/length(x)
    
}

#' Calculate proprotion of cells in each cluster where a gene is detected
#'
#' @param seurat Seurat object
#' @param cluster_col Character, name of column in \code{seurat@meta.data} where
#' cluster identities are stored
#' @param genes Character, gene or genes of interest
#'
#' @return
#' @export
calc_pct1 <- function(seurat,
                      cluster_col,
                      genes) {
    
    require(data.table)
    
    message("## Processing ", seurat@project.name)
    
    seurat <- SetAllIdent(seurat, cluster_col)
    
    x <- as.matrix(seurat@data)[genes, ]
    # Binarize
    x[x > 0] <- 1
    x <- as.data.table(t(x))
    # Add cluster info for cells
    x[, Cluster := as.character(seurat@ident)]
    # Get prop of cells expressing a gene, within each cluster
    x[, lapply(.SD, prop), by = Cluster]
    
}



#' This function is adapted from Samantha Worme, adapted from the monocle package
ordered_pseudotime_heatmap <- function(cds_subset,
                                       order_by_pseudotime = TRUE,
                                       
                                       cluster_rows = FALSE,
                                       hclust_method = "ward.D2",
                                       num_clusters = 6,
                                       hmcols = NULL,
                                       add_annotation_row = NULL,
                                       add_annotation_col = NULL,
                                       show_rownames = FALSE,
                                       use_gene_short_name = TRUE,
                                       norm_method = c("log", "vstExprs"),
                                       scale_max = 3,
                                       scale_min = -3,
                                       trend_formula = '~sm.ns(Pseudotime, df=3)',
                                       prefix = "./",
                                       cores=1,
                                       return_curves = FALSE,
                                       return_maxes = FALSE) {
    
    num_clusters <- min(num_clusters, nrow(cds_subset))
    pseudocount <- 1
    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), max(pData(cds_subset)$Pseudotime),length.out = 100))
    
    # Get smoothed expression curves using Monocle
    m <- monocle::genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,
                                  relative_expr = T, new_data = newdata)
    
    if (return_curves) return(m)
    
    # Remove genes with no expression in any condition
    m=m[!apply(m,1,sum)==0,]
    
    norm_method <- match.arg(norm_method)
    
    if(norm_method == 'vstExprs' && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
        m = vstExprs(cds_subset, expr_matrix=m)
    }
    else if(norm_method == 'log') {
        m = log10(m + pseudocount)
    }
    
    # Row-center the data i.e. make mean zero
    m=m[!apply(m,1,sd)==0,]
    m=Matrix::t(scale(Matrix::t(m),center=TRUE))
    m=m[is.na(row.names(m)) == FALSE,]
    m[is.nan(m)] = 0
    m[m>scale_max] = scale_max
    m[m<scale_min] = scale_min
    
    # Order matrix by pseudotime position of each gene's expression maximum
    time_max <- max.col(m)
    m <- cbind(m, time_max)
    m <- m[order(m[ , 101]), ]
    
    if (return_maxes) return(m)
    
    # Remove the column with value of the maxes
    m <- m[ , -101]
    
    # If we don't order by pseudotime, then order by input
    if (!order_by_pseudotime) m <- m[rownames(cds_subset), ]
    
    heatmap_matrix <- m
    
    # Calculate distance matrix
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    
    if(is.null(hmcols)) {
        bks <- seq(-3.1,3.1, by = 0.1)
        hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    }
    else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    
    hm_fun <- purrr::partial(pheatmap,
                             heatmap_matrix,
                             useRaster = T,
                             cluster_cols = FALSE,
                             cluster_rows = cluster_rows,
                             show_rownames = show_rownames,
                             show_colnames = FALSE,
                             clustering_distance_rows = row_dist,
                             clustering_method = hclust_method,
                             cutree_rows = num_clusters,
                             silent = TRUE,
                             breaks = bks,
                             cellheight = 13,
                             fontsize = 6,
                             border_color = NA,
                             color = hmcols)
    
    hm_fun(filename = paste0(prefix, ".png"))
    hm_fun(filename = paste0(prefix, ".pdf"))
    
    
}

