#' @title Extract Differentially Expressed Genes (DEGs) Lists
#' @description Extract significantly upregulated and downregulated genes from multiple DEG analysis results.
#' @param DEGs_lists A list of dataframes containing DEG analysis results.
#' @param logFC_cut A numeric value specifying the log fold change cutoff for significant DEGs. Default is 1.
#' @param p_cut A numeric value specifying the p-value cutoff for significant DEGs. Default is 0.05.
#' @return A list containing two lists: one for upregulated genes and one for downregulated genes across the provided datasets.
#' @examples
#' \dontrun{
#' df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
#' results1 <- DEGs_analysis(df1)
#' df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
#' results2 <- DEGs_analysis(df2)
#' DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
#' results <- get_DEGs_list(DEGs_lists)
#' }
#' @export
get_DEGs_list <- function(DEGs_lists, logFC_cut = 1, p_cut = 0.05) {
  # Helper function to extract upregulated and downregulated genes based on cutoffs
  extract_genes <- function(DEG_data, logFC_cut = 1, p_cut = 0.05) {
    # Remove rows with NA values
    DEG_data <- DEG_data[complete.cases(DEG_data), ]

    # Extract significantly upregulated genes
    up_genes <- DEG_data[DEG_data$logFC > logFC_cut & DEG_data$P.Value < p_cut, ]
    up_genes <- up_genes[order(up_genes$logFC, decreasing = TRUE), "gene"]

    # Extract significantly downregulated genes
    down_genes <- DEG_data[DEG_data$logFC < -logFC_cut & DEG_data$P.Value < p_cut, ]
    down_genes <- down_genes[order(down_genes$logFC, decreasing = FALSE), "gene"]

    # Return a list containing upregulated and downregulated genes
    return(list(up = up_genes, down = down_genes))
  }

  # Initialize lists to store upregulated and downregulated genes
  glist_up <- list()
  glist_down <- list()

  # Names of the datasets in DEGs_lists
  gse.list <- names(DEGs_lists)

  # Iterate over each dataset in DEGs_lists
  for (i in seq_along(DEGs_lists)) {
    genes <- extract_genes(DEGs_lists[[i]], logFC_cut, p_cut)
    glist_up[[gse.list[i]]] <- genes$up
    glist_down[[gse.list[i]]] <- genes$down
  }

  # Return a list containing all upregulated and downregulated genes
  return(list(DEG_up = glist_up, DEG_down = glist_down))
}
