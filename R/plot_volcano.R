#' @title Volcano plot for DEGs
#' @description
#' Plotting volcano plot for DEGs between tumor and normal samples in CPTAC datasets.
#' @import dplyr ggplot2 ggrepel
#' @importFrom dplyr arrange
#' @importFrom ggrepel geom_text_repel
#' @param results DataFrame. The results from DEGs analysis containing columns 'adj.P.Val', 'P.Value', 'logFC', and 'gene'.
#' @param p.cut Numeric. The cutoff for adjusted p-value to determine significance. Default is 0.05.
#' @param logFC.cut Numeric. The cutoff for log fold change to determine significance. Default is 1.
#' @param show.top Logical. If TRUE, labels the top 5 up- and downregulated genes. Default is FALSE.
#' @param show.labels Character vector. Specific gene labels to show. Default is NULL.
#' @param colors A vector of color panel, default c("blue", "grey20", "red").
#' @return A ggplot2 object representing the volcano plot.
#' @examples
#' \dontrun{
#' df <- get_OSF_data(table = "GSE74706", action = "geo_data")
#' results <- DEGs_analysis(df)
#' plot_volcano(results)
#' }
#' @export
plot_volcano <- function(results, p.cut = 0.05, logFC.cut = 1,
                         show.top = FALSE, show.labels = NULL,
                         colors = c("blue", "grey20", "red")) {

  # Ensure numeric types
  results$adj.P.Val <- as.numeric(results$adj.P.Val)
  results$P.Value <- as.numeric(results$P.Value)
  results$logFC <- as.numeric(results$logFC)

  # Determine change status
  results$change <- ifelse(results$adj.P.Val < p.cut,
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

  # Arrange results by logFC
  results <- dplyr::arrange(results, logFC)

  # Label top genes if show.top is TRUE
  if (show.top) {
    top5 <- c(results[c(1:5), 1], results[c((nrow(results) - 4):nrow(results)), 1])
    results$label <- ifelse(results$gene %in% top5, results$gene, "")
  }

  # Label specific genes if show.labels is provided
  if (!is.null(show.labels)) {
    results$label <- ifelse(results$gene %in% show.labels, results$gene, "")
  }

  # Create ggplot object
  p <- ggplot2::ggplot(data = results) +
    ggplot2::aes(x = logFC, y = -log10(adj.P.Val)) +
    ggplot2::geom_point(size = 2, shape = 16, alpha = 0.8) +
    ggplot2::theme_light(base_size = 20) +
    ggplot2::aes(color = change) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = -log10(p.cut), linetype = "dashed", color = "grey30", size = 0.8) +
    ggplot2::geom_vline(xintercept = c(logFC.cut, -logFC.cut), linetype = "dashed", color = "grey30", size = 0.8) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = bquote(''*Log[2]*' (Fold change)'), y = bquote(''*-Log[10]*' (P.adj value)'))

  # Add text labels if required
  if (show.top | !is.null(show.labels)) {
    p <- p + ggrepel::geom_text_repel(
      data = results,
      ggplot2::aes(label = label),
      size = 5,
      color = "black",
      nudge_x = 0.2,
      nudge_y = 0.2,
      max.overlaps = 100000,
      box.padding = grid::unit(0.9, "lines"),
      point.padding = grid::unit(0.8, "lines"),
      show.legend = FALSE
    )
  }

  # Return the plot
  return(p)
}
