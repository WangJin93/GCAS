#' @title Perform correlation analysis
#' @description Perform correlation analysis of the expression data in multiple datasets.
#' @import dplyr psych
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param geneset_data The expression data of a genelist in multiple datasets, obtained by the get_expr_data() function.
#' @param tumor_subtype Tumor subtype used for correlation analysis, default is NULL.
#' @param sample_type Sample type used for correlation analysis, default all types: c("Tumor", "Normal").
#' @param cor_method Method for correlation analysis, default "pearson".
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
                              cor_method = "pearson") {

  # Return NULL if the input dataframe is NULL
  if (is.null(df)) {
    return(NULL)
  }

  # Merge the target gene expression data with the genelist expression data
  df <- merge(df, geneset_data,
              by = c("ID", "subtype", "dataset", "tissue", "Patient.ID"))

  # Determine the sample type based on the subtype and tumor_subtype
  if (is.null(tumor_subtype)) {
    df$type <- ifelse(df$subtype %in% c("Normal", "Adjacent"), "Normal", "Tumor")
  } else {
    tumor_subtype <- extract_subset(subtype,tumor_subtype)
    df <- df %>% dplyr::filter(subtype %in% c(tumor_subtype, "Normal", "Adjacent"))
    df$type <- ifelse(df$subtype %in% tumor_subtype, "Tumor", "Normal")
  }

  # Filter the dataframe based on the sample type
  df <- df %>% dplyr::filter(type %in% sample_type)
  df[6:(ncol(df)-1)] <- apply(df[6:(ncol(df)-1)], 2, as.numeric)

  # Split the dataframe by dataset
  sss <- split(df, df$dataset)
  dataset <- names(sss)
  nrow <- length(dataset)
  sig <- colnames(df)[7:(ncol(df)-1)]
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
      corr_result <- psych::corr.test(x = sss_can[, colnames(df)[6]],
                                      y = sss_can[, sig],
                                      method = cor_method)
      rvalue[i, ] <- corr_result[["r"]]
      pvalue[i, ] <- corr_result[["p"]]
    }
  }

  # Transpose correlation and p-value matrices
  rvalue_T <- t(rvalue)
  pvalue_T <- t(pvalue)

  # Create a list to hold the results
  plist <- list(r = rvalue_T, p = pvalue_T, sss = sss)

  return(plist)
}
