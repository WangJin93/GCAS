#' @title Plot logFC scatter
#' @description
#' This function generates a scatter of log fold changes (logFC) for a gene across different datasets.
#' It includes significance annotations based on p-values.
#' @import ggplot2
#' @import dplyr
#' @param results Data frame. The results data frame containing columns for dataset, gene, logFC, and p.value.
#' @param p.cut Numeric. The cutoff for adjusted p-value to determine significance. Default is 0.05.
#' @param logFC.cut Numeric. The cutoff for log fold change to determine significance. Default is 1.
#' @param colors A vector of color panel, default c("blue", "grey20", "red").
#' @return A ggplot object representing the scatter.
#' @examples
#' \dontrun{
#' df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
#' results <- data_summary(df, tumor_subtype = "LUAD")
#' scatter <- plot_logFC_scatter(results, logFC.cut = 0.5, colors = c("blue","grey20", "red"))
#' print(scatter)
#' }
#' @export
plot_logFC_scatter <- function(results, p.cut = 0.05,
                               logFC.cut = 1, colors = c("blue", "grey20", "red")) {

  # Remove rows with NA values
  results <- na.omit(results)
  results$p.value <- as.numeric(results$p.value)
  results$logFC <- as.numeric(results$logFC)

  # Determine change status
  results$change <- ifelse(results$p.value < p.cut,
                           ifelse(results$logFC < -logFC.cut, "Down",
                                  ifelse(results$logFC > logFC.cut, "Up", "No")),
                           "No")
  print(table(results$change))
  color_map <- c(
    "Down" = colors[1],
    "No" = colors[2],
    "Up" = colors[3]
  )

  # 为 unique_values 分配颜色
  colors <- color_map[unique(results$change)]
  # Create ggplot object
  p <- ggplot2::ggplot(data = results) +
    ggplot2::aes(x =logFC, y = -log10(p.value)) +
    ggplot2::geom_point(size = 2, shape = 16, alpha = 0.8) +
    ggplot2::theme_light(base_size = 20) +
    ggplot2::aes(color = change) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = -log10(p.cut), linetype = "dashed", color = "grey30", size = 0.8) +
    ggplot2::geom_vline(xintercept = c(logFC.cut, -logFC.cut), linetype = "dashed", color = "grey30", size = 0.8) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = bquote(''*Log[2]*' (Fold change)'), y = bquote(''*-Log[10]*' (P value)'))

  return(p)
}
