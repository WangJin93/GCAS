#' @title Correlation Analysis of Drug Sensitivity
#' @description Calculate the correlation between target gene expression and anti-tumor drug sensitivity in multiple datasets.
#' @import dplyr
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param cor_method Method for correlation analysis, default "pearson".
#' @param Target.pathway The signaling pathways of anti-tumor drug targets, default "Cell cycle". Use "drug_info" to get the detailed information of these drugs.
#' @param adjust_method Method for adjusting p-values for multiple testing. Options include "none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param conf_level The confidence level for confidence interval calculation (default: 0.95).
#' @return A list containing the correlation results (including n, r, t, p, p.adj, ci_lower, ci_upper) and the merged data.
#' @examples
#' \dontrun{
#' dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210",
#'              "GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072",
#'              "GSE74706","GSE18842","GSE62113")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' result <- cor_gcas_drug(df, Target.pathway = c("Cell cycle"))
#' }
#' @export
cor_gcas_drug <- function(df,
                          cor_method = "pearson",
                          Target.pathway = c("Cell cycle"),
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

  # Filter out normal and adjacent subtypes and join with GCAS drug data
  df <- df %>%
    dplyr::filter(!subtype %in% c("Normal", "Adjacent")) %>%
    dplyr::inner_join(GCAS_drug, by = "ID")

  # Extract target drugs from drug_info based on Target.pathway
  sig <- drug_info[drug_info$Target.pathway %in% Target.pathway, ][c("ID", "Name")]
  sig$drug <- paste(sig$Name, sig$ID, sep = "_")
  sig <- sig$drug
  
  if (length(sig) == 0) {
    warning("No drugs found for the specified Target.pathway")
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
