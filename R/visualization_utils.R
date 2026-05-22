#' @title GCAS ggplot2 theme
#' @description A consistent theme for all GCAS visualizations.
#' @import ggplot2
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @export
#' @examples
#' \dontrun{
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   gcas_theme()
#' }
gcas_theme <- function(base_size = 16, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = base_size),
      axis.text.y = element_text(size = base_size),
      axis.title.x = element_text(size = base_size),
      axis.title.y = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.8),
      legend.title = element_text(size = base_size),
      strip.text = element_text(size = base_size),
      plot.title = element_text(size = base_size * 1.2, hjust = 0.5)
    )
}

#' @title GCAS color palette
#' @description Default color palette for GCAS visualizations.
#' @param n Number of colors needed.
#' @param palette Type of palette: "main", "diverging", or "sequential".
#' @return A vector of colors.
#' @export
#' @examples
#' \dontrun{
#' ggplot(mtcars, aes(mpg, wt, color = factor(cyl))) +
#'   geom_point(size = 3) +
#'   scale_color_manual(values = gcas_palette(3))
#' }
gcas_palette <- function(n = 3, palette = "main") {
  palettes <- list(
    main = c("#00AFBB", "#FC4E07", "#E7B800", "#2E9FDF", "#72B7B2"),
    diverging = c("#5C5DAF", "white", "#EA2E2D"),
    sequential = colorRampPalette(c("#E7B800", "#FC4E07"))(n)
  )
  
  if (!palette %in% names(palettes)) {
    warning(paste("Palette", palette, "not found. Using 'main' palette."))
    palette <- "main"
  }
  
  if (palette == "main" && n > length(palettes$main)) {
    warning(paste("Requested", n, "colors but main palette only has", 
                  length(palettes$main), ". Recycling colors."))
    return(rep(palettes$main, length.out = n))
  }
  
  if (palette == "diverging") {
    return(palettes$diverging)
  }
  
  if (palette == "sequential") {
    return(palettes$sequential)
  }
  
  return(palettes$main[1:n])
}

#' @title Add significance annotations to plot
#' @description Helper function to add significance labels to ggplot.
#' @import ggplot2
#' @param p A ggplot object.
#' @param pv A data frame containing p-values and significance labels.
#' @param show.label Logical, whether to show significance labels.
#' @param show.value Logical, whether to show p-values.
#' @param y.position Y position for annotations.
#' @return A ggplot object with annotations.
#' @keywords internal
.add_significance <- function(p, pv, show.label = TRUE, show.value = FALSE, 
                              y.position = NULL) {
  
  if (is.null(pv) || nrow(pv) == 0) {
    return(p)
  }
  
  if (is.null(y.position)) {
    # Extract y range from plot data
    plot_data <- ggplot2::ggplot_build(p)$data[[1]]
    y_range <- range(plot_data$y, na.rm = TRUE)
    y.position <- y_range[2] + diff(y_range) * 0.1
  }
  
  if (show.label && "p.signif" %in% colnames(pv)) {
    p <- p + 
      geom_text(
        aes(label = .data$p.signif),
        data = pv,
        y = y.position,
        size = 8,
        inherit.aes = FALSE
      )
  }
  
  if (show.value && "p" %in% colnames(pv)) {
    p <- p + 
      geom_text(
        aes(label = ifelse(.data$p < 0.001, "P < 0.001", 
                          paste0("P = ", signif(.data$p, 3)))),
        data = pv,
        y = y.position,
        size = 6,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}

#' @title Add sample count labels to plot
#' @description Helper function to add sample count labels to ggplot.
#' @import ggplot2
#' @param p A ggplot object.
#' @param count_data A data frame containing sample counts.
#' @param y.position Y position for labels.
#' @param color_values Color values for labels.
#' @return A ggplot object with sample count labels.
#' @keywords internal
.add_sample_counts <- function(p, count_data, y.position = NULL, color_values = NULL) {
  
  if (is.null(count_data) || nrow(count_data) == 0) {
    return(p)
  }
  
  if (is.null(y.position)) {
    plot_data <- ggplot2::ggplot_build(p)$data[[1]]
    y_range <- range(plot_data$y, na.rm = TRUE)
    y.position <- y_range[1] - diff(y_range) * 0.2
  }
  
  p <- p +
    geom_text(
      data = count_data,
      aes(label = .data$n, color = .data$type),
      y = y.position,
      size = 6,
      inherit.aes = FALSE
    )
  
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values)
  }
  
  return(p)
}