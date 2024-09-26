#' @title Correlation Analysis of Drug Sensitivity
#' @description Calculate the correlation between target gene expression and anti-tumor drug sensitivity in multiple datasets.
#' @import dplyr psych
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param cor_method Method for correlation analysis, default "pearson".
#' @param Target.pathway The signaling pathways of anti-tumor drug targets, default "Cell cycle". Use "drug_info" to get the detailed information of these drugs.
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
                          Target.pathway = c("Cell cycle")) {

  # Check if the dataframe is NULL, return NULL if so
  if (is.null(df)) {
    return(NULL)
  }

  # Filter out normal and adjacent subtypes and join with GCAS drug data
  df <- df %>%
    dplyr::filter(!subtype %in% c("Normal", "Adjacent")) %>%
    dplyr::inner_join(GCAS_drug, by = "ID")

  # Split the dataframe by dataset
  sss <- split(df, df$dataset)
  dataset <- names(sss)
  nrow <- length(dataset)

  # Extract target drugs from drug_info based on Target.pathway
  sig <- drug_info[drug_info$Target.pathway %in% Target.pathway, ][c("ID", "Name")]
  sig$drug <- paste(sig$Name, sig$ID, sep = "_")
  sig <- sig$drug
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

  # Transpose and remove rows with NA values
  rvalue_T <- t(na.omit(rvalue))
  pvalue_T <- t(na.omit(pvalue))

  # Create a list to hold the results
  plist <- list(r = rvalue_T, p = pvalue_T, sss = sss)

  return(plist)
}
