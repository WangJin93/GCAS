#' @title Perform correlation analysis
#' @description Perform correlation analysis of the mRNA/protein expression data in CPTAC database.
#' @import dplyr stringr psych
#' @param dataset Dataset name. Use 'dataset$Abbre' to get all datasets.
#' @param id1 Gene symbol, you can input one gene symbol.
#' @param id2 Gene symbols, you can input one or multiple symbols.
#' @param tumor_subtype Tumor subtype used for correlation analysis, default is NULL.
#' @param sample_type Sample type used for correlation analysis, default all types: c("Tumor", "Normal").
#' @param cor_method Method for correlation analysis, default "pearson".
#' @param adjust_method Method for adjusting p-values for multiple testing. Options include "none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param conf_level The confidence level for confidence interval calculation (default: 0.95).
#' @return A list containing the correlation results (including n, r, t, p, p.adj, ci_lower, ci_upper) and the merged data.
#' @examples
#' \dontrun{
#' results <- cor_cancer_genelist(dataset = "GSE62113",
#'                                id1 = "STAT3",tumor_subtype = "LC",
#'                                id2 = c("TNS1", "TP53"),
#'                                sample_type = c("Tumor", "Normal"),
#'                                cor_method = "pearson",
#'                                adjust_method = "BH")
#' }
#' @export
cor_cancer_genelist <- function(dataset = "GSE62113",
                                id1 = "STAT3",
                                id2 = c("TNS1", "TP53"),
                                tumor_subtype = NULL,
                                sample_type = c("Tumor", "Normal"),
                                cor_method = "pearson",
                                adjust_method = "BH",
                                conf_level = 0.95) {

  # Input validation
  if (!is.character(dataset) || length(dataset) != 1) {
    stop("dataset must be a single character string")
  }
  
  if (!is.character(id1) || length(id1) != 1) {
    stop("id1 must be a single character string")
  }
  
  if (!is.character(id2)) {
    stop("id2 must be a character vector")
  }
  
  if (!adjust_method %in% c("none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")) {
    stop("adjust_method must be one of: 'none', 'BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni'")
  }

  # Retrieve expression data for the first gene
  data1 <- get_expr_data(dataset, id1)

  # Check if data1 is NULL, return NULL if so
  if (is.null(data1)) {
    warning(paste("No data found for gene", id1, "in dataset", dataset))
    return(NULL)
  }

  # Retrieve expression data for the second gene(s)
  data2 <- get_expr_data(dataset, id2)
  
  if (is.null(data2)) {
    warning(paste("No data found for genes", paste(id2, collapse = ", "), "in dataset", dataset))
    return(NULL)
  }

  # Define common columns for merging
  common_cols <- c("ID", "subtype", "dataset", "tissue", "Patient.ID")
  
  # Merge the two datasets on common columns
  data_merge <- merge(data1, data2, by = common_cols)

  # Determine sample type using utility function
  data_merge <- .determine_sample_type(data_merge, tumor_subtype)

  # Filter data based on sample_type
  data_merge <- data_merge %>% 
    dplyr::filter(type %in% sample_type)

  # Get the gene columns (excluding non-expression columns)
  non_expr_cols <- c(common_cols, "type")
  gene_cols <- setdiff(colnames(data_merge), non_expr_cols)
  
  # Convert expression columns to numeric using dplyr
  data_merge <- data_merge %>% 
    dplyr::mutate(dplyr::across(all_of(gene_cols), as.numeric))

  # Get target and gene set columns
  target_col <- gene_cols[1]
  sig_cols <- gene_cols[-1]
  
  # Calculate sample size
  n <- nrow(data_merge)
  
  # Perform correlation test using psych::corr.test
  result <- psych::corr.test(data_merge[[target_col]], 
                             data_merge[, sig_cols, drop = FALSE], 
                             method = cor_method)

  # Extract correlation coefficients
  r <- as.numeric(result$r)
  
  # Calculate adjusted p-values
  p_values <- as.numeric(result$p)
  adj_p_values <- if (adjust_method == "none") {
    p_values
  } else {
    stats::p.adjust(p_values, method = adjust_method)
  }
  
  # Calculate t-statistic (only for pearson correlation)
  t_values <- if (cor_method == "pearson") {
    r * sqrt((n - 2) / (1 - r^2))
  } else {
    rep(NA, length(r))
  }
  
  # Calculate confidence intervals
  ci_lower <- numeric(length(r))
  ci_upper <- numeric(length(r))
  for (j in seq_along(r)) {
    ci <- .calculate_correlation_ci(r[j], n, conf_level)
    ci_lower[j] <- ci["lower"]
    ci_upper[j] <- ci["upper"]
  }

  # Extract correlation coefficients and p-values using dplyr
  result_df <- data.frame(
    Symbol = colnames(result$r),
    n = n,
    r = r,
    t = t_values,
    p = p_values,
    p.adj = adj_p_values,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

  # Ensure symbol names are correct
  result_df$Symbol <- if (length(id2) == 1) id2 else result_df$Symbol

  # Create a list to hold the results and merged data
  p_df <- list("cor_result" = result_df, "cor_data" = data_merge)

  return(p_df)
}
