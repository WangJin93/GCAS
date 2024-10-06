\name{dataset_info}
\alias{dataset_info}
\docType{data}
\title{
Dataset informations
}
\description{
Summary of the general informations of the GEO datasets in this tool.
}
\usage{data("dataset_info")}
\format{
  A data frame with 288 observations on the following 6 variables.
  \describe{
    \item{\code{type}}{a character vector}
    \item{\code{dataset}}{a character vector}
    \item{\code{Normal_Adjacent}}{a numeric vector}
    \item{\code{Tumor}}{a numeric vector}
    \item{\code{Premalignant}}{a numeric vector}
    \item{\code{Paired}}{a character vector}
  }
}
\examples{
data(dataset_info)
str(dataset_info)
}
