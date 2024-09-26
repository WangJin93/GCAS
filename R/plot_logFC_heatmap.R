#' @title Plot logFC Heatmap
#' @description
#' This function generates a heatmap of log fold changes (logFC) for a gene across different datasets.
#' It includes significance annotations based on p-values.
#' @import ggplot2
#' @import dplyr
#' @param results Data frame. The results data frame containing columns for dataset, gene, logFC, and p.value.
#' @param direction Ploting direction, horizontal or vertical.
#' @return A ggplot object representing the heatmap.
#' @examples
#' \dontrun{
#' df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
#' results <- data_summary(df, tumor_subtype = "LUAD")
#' heatmap <- plot_logFC_heatmap(results)
#' print(heatmap)
#' }
#' @export
plot_logFC_heatmap <- function(results,direction = "horizontal") {

  # Remove rows with NA values
  results <- na.omit(results)

  # Reorder factor levels of the dataset column to match the order of appearance
  results$dataset <- factor(results$dataset, levels = unique(results$dataset))

  # Create significance annotations based on p-values
  results$text <- dplyr::case_when(
    results$p.value < 0.001 ~ "***",
    results$p.value < 0.01 ~ "**",
    results$p.value < 0.05 ~ "*",
    TRUE ~ ""  # This handles p.value >= 0.05
  )
  results <- results[order(results$logFC),]
  results$dataset <- factor(results$dataset,levels = results$dataset)
  # Generate the heatmap using ggplot2
  if (direction == "horizontal"){
    heat <- ggplot2::ggplot(results, ggplot2::aes(dataset, gene)) +
    ggplot2::geom_tile(ggplot2::aes(fill = logFC), colour = "grey", linewidth = 1) +
    ggplot2::scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1.2, vjust = 1.2, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.margin = ggplot2::margin(l = 20)
    ) +
    ggplot2::geom_text(ggplot2::aes(label = text), col = "black", size = 3) +
    ggplot2::labs(
      fill = paste0("  *   p < 0.05", "  |  ", " **  p < 0.01", "  |  ", "*** p < 0.001      ")
    ) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::scale_y_discrete(position = "right")
  }
  if (direction == "vertical"){
    heat <- ggplot2::ggplot(results, ggplot2::aes( gene,dataset)) +
      ggplot2::geom_tile(ggplot2::aes(fill = logFC), colour = "grey", linewidth = 1) +
      ggplot2::scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12)
      ) +
      ggplot2::geom_text(ggplot2::aes(label = text), col = "black", size = 3) +
      ggplot2::labs(
        fill = paste0("\n","\n","\n","\n","\n","  *   p < 0.05","\n"," **  p < 0.01","\n","*** p < 0.001","\n\n","Log2FC")
      ) +
      ggplot2::scale_x_discrete(position = "bottom") +
      ggplot2::scale_y_discrete(position = "left")
  }

  # Return the generated heatmap
  return(heat)
}
