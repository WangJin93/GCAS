#' @title Perform correlation analysis
#' @description Perform correlation analysis of the expression data in multiple datasets.
#' @import dplyr
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param geneset_data The expression data of a genelist in multiple datasets, obtained by the get_expr_data() function.
#' @param tumor_subtype Tumor subtype used for correlation analysis, default is NULL.
#' @param sample_type Sample type used for correlation analysis, default all types: c("Tumor", "Normal").
#' @param cor_method Method for correlation analysis, default "pearson".
#' @param adjust_method Method for adjusting p-values for multiple testing. Options include "none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param conf_level The confidence level for confidence interval calculation (default: 0.95).
#' @return A list containing the correlation results (including n, r, t, p, p.adj, ci_lower, ci_upper) and the merged data.
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
#' df <- get_expr_data(genes = "TNS1",datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
#' result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))
#' }
#' @export
cor_gcas_genelist <- function(df,
                              geneset_data,
                              tumor_subtype = NULL,
                              sample_type = c("Tumor", "Normal"),
                              cor_method = "pearson",
                              adjust_method = "BH",
                              conf_level = 0.95) {

  # Input validation
  if (is.null(df)) {
    warning("Input dataframe df is NULL")
    return(NULL)
  }
  
  if (is.null(geneset_data)) {
    warning("Input dataframe geneset_data is NULL")
    return(NULL)
  }
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  if (!is.data.frame(geneset_data)) {
    stop("geneset_data must be a data frame")
  }
  
  if (!adjust_method %in% c("none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")) {
    stop("adjust_method must be one of: 'none', 'BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni'")
  }
  
  required_cols <- c("ID", "subtype", "dataset", "tissue", "Patient.ID")
  if (!all(required_cols %in% colnames(df))) {
    stop(paste("df must contain columns:", paste(required_cols, collapse = ", ")))
  }
  
  if (!all(required_cols %in% colnames(geneset_data))) {
    stop(paste("geneset_data must contain columns:", paste(required_cols, collapse = ", ")))
  }

  # Merge the target gene expression data with the genelist expression data
  df <- merge(df, geneset_data,
              by = c("ID", "subtype", "dataset", "tissue", "Patient.ID"))

  # Determine the sample type based on the subtype and tumor_subtype
  df <- .determine_sample_type(df, tumor_subtype)

  # Filter the dataframe based on the sample type
  df <- df %>% dplyr::filter(type %in% sample_type)
  
  # Get the gene columns (excluding non-expression columns)
  non_expr_cols <- c("ID", "subtype", "dataset", "tissue", "Patient.ID", "type")
  gene_cols <- setdiff(colnames(df), non_expr_cols)
  
  # Convert expression columns to numeric
  df <- df %>% dplyr::mutate(dplyr::across(all_of(gene_cols), as.numeric))

  # Get the target gene column (first gene column after non-expression columns)
  target_col <- gene_cols[1]
  
  # Get the geneset columns (remaining gene columns)
  sig_cols <- gene_cols[-1]
  
  if (length(sig_cols) == 0) {
    warning("No geneset columns found")
    return(NULL)
  }

  # Use the utility function for correlation calculation
  plist <- .calculate_correlation(df, target_col = target_col, 
                                   sig_cols = sig_cols, cor_method = cor_method,
                                   adjust_method = adjust_method, conf_level = conf_level)
  
  return(plist)
}
