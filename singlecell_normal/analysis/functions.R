


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
