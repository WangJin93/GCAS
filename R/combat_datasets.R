#' @title Batch Correction using ComBat
#' @description This function performs batch correction on multiple datasets using the ComBat function from the sva package.
#' @param tables A character vector of table names to be processed.
#' @param tumor_subtype A character string specifying the tumor subtype to filter the datasets. If NULL, all subtypes are included.
#' @return A list containing the combined and batch-corrected data matrix and the sample information.
#' @import sva
#' @import dplyr
#' @import tibble
#' @export
#' @examples
#' \dontrun{
#'   tables <- c("GSE31210", "GSE74706")
#'   result <- combat_datasets(tables, tumor_subtype = "LC")
#'   combined_data <- result$combined_data
#'   sample_info <- result$sample_info
#' }
combat_datasets <- function(tables, tumor_subtype = NULL) {

  # Initialize a list to store individual datasets
  data_list <- list()

  # Read and process each dataset
  for (table in tables) {
    # Get the dataset using the get_OSF_data function
    data <- get_OSF_data(table = table, action = "geo_data")
    data <- data[!duplicated(data$ID), ]
    rownames(data) <- NULL
    # Set the first column as row names and remove it from data
    data <- data %>%
      column_to_rownames(var = "ID") %>%
      as.data.frame()

    # Append the dataset to the list
    data_list[[table]] <- data
  }

  # Find the intersection of all gene names
  common_genes <- Reduce(intersect, lapply(data_list, rownames))

  # Align each dataset to common genes
  for (table in tables) {
    data_list[[table]] <- data_list[[table]][common_genes, , drop = FALSE]
  }

  # Combine all datasets by columns
  combined_data <- do.call(cbind, data_list)
  colnames(combined_data) <- stringr::str_remove(colnames(combined_data), paste0(paste0(names(data_list), collapse = ".|"), "."))

  # Extract sample information
  tumor_subtype2 <- c(extract_subset(subtype, tumor_subtype), "Normal")

  if (is.null(tumor_subtype2)) {
    sample_info <- sample_subtype %>%
      filter(ID %in% colnames(combined_data)) %>%
      select(dataset, ID, subtype)
  } else {
    sample_info <- sample_subtype %>%
      filter(ID %in% colnames(combined_data)) %>%
      filter(subtype %in% tumor_subtype2) %>%
      select(dataset, ID, subtype)
  }

  dd <- sample_info %>% filter(subtype != "Normal") %>% .[,"dataset"]

  if (length(unique(dd)) != length(tables)) {
    cat(paste0("Warning! Only the following datasets contain the input subtype: ", paste(unique(dd), collapse = ","), "\n"))
  }

  if (length(unique(dd)) < 2) {
    cat(paste0("ERROR! The input dataset contains ", tumor_subtype, " data that is < 2, confirm your input!\n"))
    return(NULL)
  }

  sample_info$type <- ifelse(sample_info$subtype == "Normal", "Normal",
                             ifelse(sample_info$subtype %in% c("Adenoma", "polyp", "Cirrhosis", "HSIL", "CIN", "MGUS", "SMM"),
                                    "Premalignant", "Tumor"))
  sample_info <- sample_info %>% dplyr::filter(type != "Premalignant") %>% tibble::column_to_rownames("ID")

  combined_data <- combined_data[rownames(sample_info)]

  # Prepare the data for ComBat
  # Extract the batch variable (assuming dataset names are used as batches)
  mod <- model.matrix(~ as.factor(type), data = sample_info)

  final_data <- ComBat(dat = combined_data,
                       batch = sample_info$dataset,
                       mod = mod,
                       par.prior = TRUE)

  return(list(combined_data = final_data,
              sample_info = sample_info))
}

