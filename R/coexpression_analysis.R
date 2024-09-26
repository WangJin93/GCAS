#' @title Calculate Gene Correlations
#' @description This function calculates the correlation between a given gene and all other genes
#' in the provided expression matrix. It also provides the corresponding p-values.
#' @import tibble
#' @param expression_matrix A numeric matrix where rows represent genes and columns represent samples.
#' @param gene A character string representing the gene for which the correlations will be calculated.
#' @param method A character string specifying the correlation method to be used. Default is "pearson".
#' @return A data frame containing gene names, correlation coefficients, and p-values.
#' @export
#' @examples
#' \dontrun{
#'   expression_matrix <- get_OSF_data(table = "GSE74706", action = "geo_data")
#'   results <- coexpression_analysis(expression_matrix, "RPN1")
#'   print(results)
#' }
coexpression_analysis <- function(expression_matrix, gene, method = "pearson") {
  # Remove duplicated rows based on the 'ID' column
  expression_matrix <- expression_matrix[!duplicated(expression_matrix$ID), ]

  # Ensure row names are properly set
  rownames(expression_matrix) <- NULL
  expression_matrix <- column_to_rownames(expression_matrix, "ID")

  # Check if the gene of interest exists in the expression matrix
  if (!(gene %in% rownames(expression_matrix))) {
    stop("The gene of interest does not exist in the expression matrix")
  }

  # Convert the R matrix to a format usable in C++
  expression_matrix_cpp <- as.matrix(expression_matrix)
  rownames_cpp <- rownames(expression_matrix)

  # Call a C++ function for correlation calculation (assuming cal_correlation is defined in C++)
  result <- cal_correlation(expression_matrix_cpp, rownames_cpp, gene, method) %>%
    na.omit()

  return(result)
}
