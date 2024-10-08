% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_logFC_heatmap.R
\name{plot_logFC_heatmap}
\alias{plot_logFC_heatmap}
\title{Plot logFC Heatmap}
\usage{
plot_logFC_heatmap(results, direction = "horizontal")
}
\arguments{
\item{results}{Data frame. The results data frame containing columns for dataset, gene, logFC, and p.value.}

\item{direction}{Ploting direction, horizontal or vertical.}
}
\value{
A ggplot object representing the heatmap.
}
\description{
This function generates a heatmap of log fold changes (logFC) for a gene across different datasets.
It includes significance annotations based on p-values.
}
\examples{
\dontrun{
df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
heatmap <- plot_logFC_heatmap(results)
print(heatmap)
}
}
