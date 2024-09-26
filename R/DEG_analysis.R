#' @title Differential Expression Gene (DEG) Analysis
#' @description Perform differential expression gene analysis on a given dataset.
#' @param df A dataframe containing gene expression data with sample IDs as columns.
#' @param tumor_subtype A character vector specifying the tumor subtypes to be analyzed. Default is NULL, which means all tumor subtypes will be included.
#' @param ... Additional arguments passed to `lmFit`, `contrasts.fit`, and `eBayes`.
#' @return A dataframe with DEG analysis results, including log fold changes and p-values.
#' @examples
#' \dontrun{
#' df <- get_OSF_data(table = "GSE74706", action = "geo_data")
#' results <- DEGs_analysis(df, tumor_subtype = c("NSCLC"))
#' }
#' @export
DEGs_analysis <- function(df, tumor_subtype = NULL, ...) {
  # Return NULL if the input dataframe is NULL
  if (is.null(df)) {
    return(NULL)
  }

  # Filter sample information based on the provided sample IDs
  sample_info <- sample_subtype[sample_subtype$ID %in% colnames(df), ]

  # Assign type based on tumor_subtype
  if (is.null(tumor_subtype)) {
    sample_info$type2 <- ifelse(sample_info$subtype %in% c("Normal", "Adjacent"), "Normal", "Tumor")
  } else {
    tumor_subtype <- extract_subset(subtype,tumor_subtype)
    sample_info <- sample_info %>% filter(subtype %in% c(tumor_subtype, "Normal", "Adjacent"))
    sample_info$type2 <- ifelse(sample_info$subtype %in% tumor_subtype, "Tumor", "Normal")
  }

  # Subset df to include only relevant sample columns
  df <- df[, c("ID", sample_info$ID)]

  # Define group and design matrix for linear modeling
  group <- sample_info$type2
  design <- model.matrix(~0 + factor(group))
  colnames(design) <- levels(factor(group))

  # Create contrast matrix for Tumor vs Normal comparison
  contrast.matrix <- makeContrasts(contrasts = "Tumor - Normal", levels = design)

  # Fit linear model
  fit <- lmFit(df, design, ...)

  # Apply contrasts and compute eBayes statistics
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, ...)

  # Extract top differentially expressed genes
  DEG <- topTable(fit2, coef = 1, n = Inf, ...) %>% na.omit()

  # Rename the first column to "gene"
  colnames(DEG)[1] <- "gene"

  return(DEG)
}
