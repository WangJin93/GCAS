#' @title Correlation analysis of immune cells infiltration
#' @description Calculate the correlation between target gene expression and immune cells infiltration in multiple datasets.
#' @import dplyr psych
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param cor_method Method for correlation analysis, default "spearman".
#' @param TIL_type Algorithm for calculating immune cell infiltration, default "TIMER".
#' @examples
#' \dontrun{
#' dataset <- c("GSE27262", "GSE7670", "GSE19188", "GSE19804", "GSE30219",
#'              "GSE31210", "GSE32665", "GSE32863", "GSE43458", "GSE46539",
#'              "GSE75037", "GSE10072", "GSE74706", "GSE18842", "GSE62113")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' result <- cor_gcas_TIL(df, cor_method = "spearman", TIL_type = "TIMER")
#' }
#' @export
cor_gcas_TIL <- function(df,
                         cor_method = "spearman",
                         TIL_type = "TIMER") {

  # Return NULL if input dataframe is NULL
  if (is.null(df)) {
    return(NULL)
  }

  # Filter and join the target gene expression data with immune cell infiltration data
  df <- df %>%
    dplyr::filter(!subtype %in% c("Normal", "Adjacent")) %>%
    inner_join(., GCAS_TIL, by = "ID")

  # Split the dataframe by dataset
  sss <- split(df, df$dataset)
  dataset <- names(sss)
  nrow <- length(dataset)

  # Get the immune cell types based on the specified algorithm
  sig <- TIL_map[TIL_map$algorithm %in% TIL_type,]$cell_type
  ncol <- length(sig)

  # Initialize matrices to store correlation coefficients and p-values
  rvalue <- matrix(nrow = nrow, ncol = ncol)
  rownames(rvalue) <- dataset
  colnames(rvalue) <- sig
  pvalue <- matrix(nrow = nrow, ncol = ncol)
  rownames(pvalue) <- dataset
  colnames(pvalue) <- sig

  # Calculate correlation for each dataset
  for (i in seq_along(dataset)) {
    sss_can <- sss[[i]]
    if (nrow(sss_can) < 4) {
      next
    } else {
      corr_result <- psych::corr.test(x = sss_can[, colnames(df)[3]],
                                      y = sss_can[, sig],
                                      method = cor_method)
      rvalue[i, ] <- corr_result[["r"]]
      pvalue[i, ] <- corr_result[["p"]]
    }
  }

  # Transpose correlation and p-value matrices, and omit any NA values
  rvalue_T <- t(na.omit(rvalue))
  pvalue_T <- t(na.omit(pvalue))

  # Create a list to hold the results
  plist <- list(r = rvalue_T, p = pvalue_T, sss = sss)

  return(plist)
}
