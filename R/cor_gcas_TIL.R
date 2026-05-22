#' @title Correlation analysis of immune cells infiltration
#' @description Calculate the correlation between target gene expression and immune cells infiltration in multiple datasets.
#' @import dplyr
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param cor_method Method for correlation analysis, default "spearman".
#' @param TIL_type Algorithm for calculating immune cell infiltration, default "TIMER".
#' @param adjust_method Method for adjusting p-values for multiple testing. Options include "none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param conf_level The confidence level for confidence interval calculation (default: 0.95).
#' @return A list containing the correlation results (including n, r, t, p, p.adj, ci_lower, ci_upper) and the merged data.
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
                         TIL_type = "TIMER",
                         adjust_method = "BH",
                         conf_level = 0.95) {

  # Input validation
  if (is.null(df)) {
    warning("Input dataframe is NULL")
    return(NULL)
  }
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  if (!"ID" %in% colnames(df)) {
    stop("df must contain an 'ID' column")
  }
  
  if (!adjust_method %in% c("none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")) {
    stop("adjust_method must be one of: 'none', 'BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni'")
  }

  # Filter and join the target gene expression data with immune cell infiltration data
  df <- df %>%
    dplyr::filter(!subtype %in% c("Normal", "Adjacent")) %>%
    dplyr::inner_join(GCAS_TIL, by = "ID")

  # Get the immune cell types based on the specified algorithm
  sig <- TIL_map[TIL_map$algorithm %in% TIL_type, ]$cell_type
  
  if (length(sig) == 0) {
    warning("No immune cell types found for the specified TIL_type")
    return(NULL)
  }

  # Get the target gene column (third column, after ID and dataset)
  target_col <- colnames(df)[3]
  
  # Use the utility function for correlation calculation
  plist <- .calculate_correlation(df, target_col = target_col, 
                                   sig_cols = sig, cor_method = cor_method,
                                   adjust_method = adjust_method, conf_level = conf_level)
  
  return(plist)
}
