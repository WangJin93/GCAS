#' @title Summarize Datasets by Tumor Subtype
#' @description This function summarizes the sample counts and pairing status for datasets based on a specified tumor subtype.
#' @param tumor_subtype A character string specifying the tumor subtype to filter the datasets. Default is "NSCLC".
#' @return A data frame summarizing the number of Normal, Tumor, and Premalignant samples, as well as the pairing status for each dataset.
#' @export
#' @examples
#' \dontrun{
#' dt_sum <- datasets_summary("NSCLC")
#' }
datasets_summary <- function(tumor_subtype = "NSCLC") {
  # Filter samples based on the tumor subtype
  dd <- sample_subtype[sample_subtype$subtype %in% extract_subset(subtype, tumor_subtype), ]
  dd_normal <- sample_subtype[(sample_subtype$dataset %in% dd$dataset) & (sample_subtype$subtype == "Normal"), ]
  dd <- rbind(dd, dd_normal)

  # Classify samples into Normal, Premalignant, and Tumor categories
  dd$type <- ifelse(dd$subtype == "Normal", "Normal",
                    ifelse(dd$subtype %in% c("Adenoma", "polyp", "Cirrhosis", "HSIL", "CIN", "MGUS", "SMM"),
                           "Premalignant", "Tumor"))

  # Load dplyr for data manipulation
  library(dplyr)

  # Summarize sample counts and pairing status for each dataset
  summary <- dd %>%
    group_by(tissue, dataset) %>%
    summarise(
      Normal_Adjacent = sum(type == "Normal"),
      Tumor = sum(type == "Tumor"),
      Premalignant = sum(type == "Premalignant"),
      Paired = ifelse(any(Patient.ID != ""), "T", "F")
    ) %>%
    ungroup()

  # Convert summary to data frame
  dataset_info <- summary %>% as.data.frame()
  return(dataset_info)
}

