#' @title Merge sample_info data
#' @description
#' Get sample_info data and merge it with expression data.
#' @import dplyr
#' @param table Character. The name of the dataset table to retrieve sample information from. Default is "GSE19188".
#' @param data_input Data frame. Expression data obtained from get_expr_data() function.
#' @return Data frame. Merged data containing both expression data and sample information.
#' @examples
#' \dontrun{
#' data_input <- get_expr_data("GSE19188", "TP53")
#' results <- merge_clinic_data("GSE19188",data_input)
#' }
#' @export
merge_clinic_data <- function(table = "GSE19188", data_input) {

  # Retrieve sample information data
  sample_info <- get_OSF_data(table, "sample_info")

  # Merge expression data with sample information, filtering out 'Normal' subtypes
  merge_data <- data_input %>%
    dplyr::filter(subtype != "Normal") %>%
    merge(., sample_info, by = "ID")

  # Return the merged data
  return(merge_data)
}
