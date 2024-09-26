#' @title Perform GSEA analysis for single gene
#' @description This function performs Gene Set Enrichment Analysis (GSEA) based on either correlation results or limma differential analysis results.
#' @param data A data frame containing gene names and corresponding values. For correlation results, the columns should be named `gene` and `r`. For limma results, the columns should be named `gene` and `logFC`.
#' @param gmt_file Path to the GMT file containing gene sets, or directly pass GO/KEGG/Reactome datasets.
#' @param pvalue_cutoff Numeric, the p-value threshold for significance. Default is 0.05.
#' @param data_type Character, type of the input data. Either "correlation" for correlation analysis results or "limma" for limma differential analysis results.
#' @return A GSEA analysis result object.
#' @importFrom dplyr arrange distinct select
#' @importFrom clusterProfiler GSEA
#' @examples
#' \dontrun{
#' df <- get_OSF_data(table = "GSE74706", action = "geo_data")
#' results <- DEGs_analysis(df,tumor_subtype =c("NSCLC"))
#' gsea_result <- GSEA_analysis(results,  gmt_file = BP_GMT_7.5.1, data_type = "limma")
#'
#' results <- coexpression_analysis(df,"RPN1")
#' gsea_result <- GSEA_analysis(results,  gmt_file = BP_GMT_7.5.1)
#' }
#' @export
GSEA_analysis <- function(data, gmt_file, pvalue_cutoff = 0.05, data_type = "correlation") {

  # Check the input data type and process data accordingly
  if (data_type == "correlation") {

    # Process correlation analysis results
    ranked_genes <- data %>%
      arrange(desc(r)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      dplyr::select(gene, r)

    # Ensure the data frame has no row names
    rownames(ranked_genes) <- NULL

    # Prepare data: Correlation values as ranking vector
    geneList <- ranked_genes$r
    names(geneList) <- ranked_genes$gene

  } else if (data_type == "limma") {

    # Process limma differential analysis results
    ranked_genes <- data %>%
      arrange(desc(logFC)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      dplyr::select(gene, logFC)

    # Ensure the data frame has no row names
    rownames(ranked_genes) <- NULL

    # Prepare data: logFC values as ranking vector
    geneList <- ranked_genes$logFC
    names(geneList) <- ranked_genes$gene

  } else {
    stop("Invalid data_type. Choose 'correlation' or 'limma'.")
  }

  # Run GSEA analysis
  gsea_result <- GSEA(
    geneList = geneList,
    TERM2GENE = gmt_file,  # Gene set file, can be GMT file or GO/KEGG/Reactome datasets
    pvalueCutoff = pvalue_cutoff,
    verbose = FALSE
  )

  # Return result
  return(gsea_result)
}
