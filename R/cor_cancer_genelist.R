#' @title Perform correlation analysis
#' @description Perform correlation analysis of the mRNA/protein expression data in CPTAC database.
#' @import dplyr stringr
#' @param dataset Dataset name. Use 'dataset$Abbre' to get all datasets.
#' @param id1 Gene symbol, you can input one gene symbol.
#' @param id2 Gene symbols, you can input one or multiple symbols.
#' @param tumor_subtype Tumor subtype used for correlation analysis, default is NULL.
#' @param sample_type Sample type used for correlation analysis, default all types: c("Tumor", "Normal").
#' @param cor_method Method for correlation analysis, default "pearson".
#' @return A list containing the correlation results and the merged data.
#' @examples
#' \dontrun{
#' results <- cor_cancer_genelist(dataset = "GSE62113",
#'                                id1 = "STAT3",tumor_subtype = "LC",
#'                                id2 = c("TNS1", "TP53"),
#'                                sample_type = c("Tumor", "Normal"),
#'                                cor_method = "pearson")
#' }
#' @export
cor_cancer_genelist <- function(dataset = "GSE62113",
                                id1 = "STAT3",
                                id2 = c("TNS1", "TP53"),
                                tumor_subtype = NULL,
                                sample_type = c("Tumor", "Normal"),
                                cor_method = "pearson") {

  # Retrieve expression data for the first gene
  data1 <- get_expr_data(dataset, id1)

  # Check if data1 is NULL, return NULL if so
  if (is.null(data1)) {
    return(NULL)
  }

  # Retrieve expression data for the second gene(s)
  data2 <- get_expr_data(dataset, id2)

  # Merge the two datasets on common columns
  data_merge <- merge(data1, data2, by = c("ID", "subtype", "dataset", "tissue", "Patient.ID"))

  # Adjust sample type based on tumor_subtype parameter
  if (is.null(tumor_subtype)) {
    data_merge$type <- ifelse(data_merge$subtype %in% c("Normal", "Adjacent"), "Normal", "Tumor")
  } else {
    tumor_subtype <- extract_subset(subtype,tumor_subtype)
    data_merge <- data_merge %>% dplyr::filter(subtype %in% c(tumor_subtype, "Normal", "Adjacent"))
    data_merge$type <- ifelse(data_merge$subtype %in% tumor_subtype, "Tumor", "Normal")
  }

  # Filter data based on sample_type
  data_merge <- data_merge %>% dplyr::filter(.data$type %in% sample_type)

  # Convert relevant columns to numeric
  data_merge[6:(ncol(data_merge)-1)] <- apply(data_merge[6:(ncol(data_merge)-1)], 2, as.numeric)

  # Perform correlation test
  result <- corr.test(data_merge[, 6], data_merge[, 7:(ncol(data_merge)-1)], method = cor_method)

  # Extract correlation coefficients and p-values
  result <- cbind(t(as.data.frame(result[1])), t(as.data.frame(result[4])))
  result <- data.frame(
    "Symbol" = stringr::str_remove(row.names(result), "r."),
    "Correlation coefficient" = result[, 1],
    "P value" = result[, 2]
  )

  if (length(id2) == 1) {
    result[1, 1] <- id2
  }

  row.names(result) <- NULL

  # Create a list to hold the results and merged data
  p_df <- list("cor_result" = result, "cor_data" = data_merge)

  return(p_df)
}
