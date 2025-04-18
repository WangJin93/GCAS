% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_DEGs_list.R
\name{get_DEGs_list}
\alias{get_DEGs_list}
\title{Extract Differentially Expressed Genes (DEGs) Lists}
\usage{
get_DEGs_list(DEGs_lists, logFC_cut = 1, p_cut = 0.05)
}
\arguments{
\item{DEGs_lists}{A list of dataframes containing DEG analysis results.}

\item{logFC_cut}{A numeric value specifying the log fold change cutoff for significant DEGs. Default is 1.}

\item{p_cut}{A numeric value specifying the p-value cutoff for significant DEGs. Default is 0.05.}
}
\value{
A list containing two lists: one for upregulated genes and one for downregulated genes across the provided datasets.
}
\description{
Extract significantly upregulated and downregulated genes from multiple DEG analysis results.
}
\examples{
\dontrun{
df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
results1 <- DEGs_analysis(df1)
df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
results2 <- DEGs_analysis(df2)
DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
results <- get_DEGs_list(DEGs_lists)
}
}
