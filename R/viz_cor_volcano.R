#' @title Volcano plot for the correlation analysis
#' @description
#' Plotting volcano plot for the correlation analysis of a specific gene from the result  of cor_gcas_genelist().
#' @import dplyr ggplot2 ggrepel
#' @importFrom dplyr arrange
#' @importFrom ggrepel geom_text_repel
#' @param cor_result DataFrame. The results from cor_gcas_genelist analysis.
#' @param item Character. Gene symbol or drug name or immune cell.
#' @param p.cut Numeric. The cutoff for adjusted p-value to determine significance. Default is 0.05.
#' @param r.cut Numeric. The cutoff for log fold change to determine significance. Default is 0.3.
#' @param colors A vector of color panel, default c("blue", "grey20", "red").
#' @return A ggplot2 object representing the volcano plot.
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
#' df <- get_expr_data(genes = "TNS1",datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
#' cor_result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))
#' viz_cor_volcano(cor_result, "LILRB4", p.cut = 0.5, r.cut = 0.1,colors = c("blue" ,"grey20", "red"))
#' }
#' @export
viz_cor_volcano <- function(cor_result, item,
                            p.cut = 0.05, r.cut = 0.3,
                         colors = c("blue", "grey20", "red")) {

  # Ensure numeric types
  results <- data.frame(r = cor_result$r[item,],
                        p = cor_result$p[item,]) %>% na.omit()

  # Determine change status
  results$change <- ifelse(results$p < p.cut,
                           ifelse(results$r < -r.cut, "negative",
                                  ifelse(results$r > r.cut, "positive", "no")),
                           "no")
  print(table(results$change))

  # 定义颜色映射
  color_map <- c(
    "negative" = colors[1],
    "no" = colors[2],
    "positive" = colors[3]
  )

  # 为 unique_values 分配颜色
  colors <- color_map[unique(results$change)]

  # Arrange results by r
  results <- dplyr::arrange(results, r)

  # Label top genes if show.top is TRUE
  results$label <- ifelse(results$change != "no", rownames(results), "")

  # Create ggplot object
  p <- ggplot2::ggplot(data = results) +
    ggplot2::aes(x =  r, y = -log10(p)) +
    ggplot2::geom_point(size = 2, shape = 16, alpha = 0.8) +
    ggplot2::theme_light(base_size = 20) +
    ggplot2::aes(color = change) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = -log10(p.cut), linetype = "dashed", color = "grey30", size = 0.8) +
    ggplot2::geom_vline(xintercept = c(r.cut, -r.cut), linetype = "dashed", color = "grey30", size = 0.8) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = bquote(''*Log[2]*' (Correlation coefficient)'), y = bquote(''*-Log[10]*' (P value)'))

  # Add text labels if required
  p <- p + ggrepel::geom_text_repel(
    data = results,
    ggplot2::aes(label = label),
    size = 5,
    color = "black",
    nudge_x = 0.2,
    nudge_y = 0.2,
    max.overlaps = 6,
    box.padding = grid::unit(0.9, "lines"),
    point.padding = grid::unit(0.8, "lines"),
    segment.size = 0.5,  # Adjust segment size to make connecting lines clearer
    show.legend = FALSE
  )

  # Return the plot
  return(p)
}
