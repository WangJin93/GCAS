#' @title Summary Statistics for Gene Expression Data
#' @description Compute summary statistics (mean, standard deviation, etc.) and perform hypothesis testing (t-test or Wilcoxon test) for gene expression data across different datasets.
#' @import dplyr
#' @param df A dataframe containing gene expression data, with columns: 'dataset', 'subtype', and the gene expression values.
#' @param tumor_subtype A character string specifying the tumor subtype to be analyzed. Default is NULL, which means all tumor subtypes will be included.
#' @param method A character string specifying the method for hypothesis testing. Options are "t.test" for t-test and "wilcox" for Wilcoxon test. Default is "t.test".
#' @return A dataframe with summary statistics and p-values for each dataset.
#' @examples
#' \dontrun{
#' df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
#' results <- data_summary(df, tumor_subtype = "LUAD")
#' }
#' @export
data_summary <- function(df, tumor_subtype = NULL, method = "t.test") {
  # Return NULL if the input dataframe is NULL
  if (is.null(df)) {
    return(NULL)
  }

  gene <- colnames(df)[3]
  colnames(df)[3] <- "value"

  # Assign type based on tumor_subtype
  if (is.null(tumor_subtype)) {
    df$type <- ifelse(df$subtype %in% c("Normal", "Adjacent"), "Normal", "Tumor")
  } else {
    tumor_subtype <- extract_subset(subtype,tumor_subtype)
    df <- df %>% filter(subtype %in% c(tumor_subtype, "Normal", "Adjacent"))
    df$type <- ifelse(df$subtype %in% tumor_subtype, "Tumor", "Normal")
  }

  df <- df %>% dplyr::filter(type %in% c("Normal", "Tumor"))
  df <- na.omit(df)

  # Compute mean, standard deviation, and sample size for each dataset and type
  summary_data <- df %>%
    group_by(dataset, type) %>%
    summarise(
      mean_value = mean(value),
      sd_value = sd(value),
      n = n(),
      .groups = "drop"
    )

  # Extract data for Normal and Tumor types
  normal_data <- summary_data %>%
    filter(type == "Normal") %>%
    dplyr::rename(mean_Normal = mean_value, sd_Normal = sd_value, n_Normal = n) %>%
    dplyr::select(-type) %>%
    mutate(se_Normal = sd_Normal / sqrt(n_Normal))

  tumor_data <- summary_data %>%
    filter(type == "Tumor") %>%
    dplyr::rename(mean_Tumor = mean_value, sd_Tumor = sd_value, n_Tumor = n) %>%
    dplyr::select(-type) %>%
    mutate(se_Tumor = sd_Tumor / sqrt(n_Tumor))

  # Merge Normal and Tumor data
  results <- normal_data %>%
    left_join(tumor_data, by = "dataset") %>%
    mutate(
      SE = sqrt(se_Normal^2 + se_Tumor^2),
      logFC = mean_Tumor - mean_Normal
    )

  # Compute p-values for each dataset
  p_values <- df %>%
    group_by(dataset) %>%
    summarise(
      p.value = if (sum(type == "Normal") > 0 & sum(type == "Tumor") > 0 & length(unique(type)) == 2) {
        if (method == "t.test") {
          tryCatch({
            t.test(value ~ type)$p.value
          }, error = function(e) NA_real_)
        } else if (method == "wilcox") {
          tryCatch({
            wilcox.test(value ~ type)$p.value
          }, error = function(e) NA_real_)
        } else {
          NA_real_
        }
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  # Merge p-values with results
  results <- results %>%
    left_join(p_values, by = "dataset") %>%
    mutate(gene = gene)

  return(results)
}
